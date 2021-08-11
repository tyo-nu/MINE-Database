from itertools import product
from pathlib import Path
from typing import List, Set

import numpy as np
import rdkit.rdBase as rkrb
import rdkit.RDLogger as rkl
import torch
import torch.utils.data
from keras.models import model_from_json
from rdkit.Chem import AllChem, MolFromSmiles, MolToSmiles, RemoveHs
from torch import nn
from torch.nn import functional as F

from minedatabase.filters.base_filter import Filter
from minedatabase.pickaxe import Pickaxe


logger = rkl.logger()
logger.setLevel(rkl.ERROR)
rkrb.DisableLog("rdApp.error")

device = torch.device("cpu")

CHARSET = [
    " ",
    "#",
    "(",
    ")",
    "+",
    "-",
    "/",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "=",
    "@",
    "B",
    "C",
    "F",
    "H",
    "I",
    "N",
    "O",
    "P",
    "S",
    "[",
    "\\",
    "]",
    "c",
    "l",
    "n",
    "o",
    "r",
    "s",
]

cwd = Path(__file__).parent


class OneHotFeaturizer(object):
    def __init__(self, charset=CHARSET, padlength=120):
        self.charset = CHARSET
        self.pad_length = padlength

    def featurize(self, smiles):
        return np.array([self.one_hot_encode(smi) for smi in smiles])

    def one_hot_array(self, i):
        return [int(x) for x in [ix == i for ix in range(len(self.charset))]]

    def one_hot_index(self, c):
        return self.charset.index(c)

    def pad_smi(self, smi):
        return smi.ljust(self.pad_length)

    def one_hot_encode(self, smi):
        return np.array(
            [self.one_hot_array(self.one_hot_index(x)) for x in self.pad_smi(smi)]
        )

    def one_hot_decode(self, z):
        z1 = []
        for i in range(len(z)):
            s = ""
            for j in range(len(z[i])):
                oh = np.argmax(z[i][j])
                s += self.charset[oh]
            z1.append([s.strip()])
        return z1

    def decode_smiles_from_index(self, vec):
        return "".join(map(lambda x: CHARSET[x], vec)).strip()


class ReactionFeasibilityFilter(Filter):
    """A filter that filters compounds based on their reaction feasibility.

    Uses the following analysis Biotechnol J. 2021 May;16(5):e2000605. doi: 10.1002/biot.202000605

    Parameters
    ----------
    use_unpredicted : Bool
        Whether or not use to unpredicted reactions, by default False
    use_wildcards : Bool
        Whether or not to allow wildcard reactions through, by default True
    generation_list : list
        Generations to apply filter -- empty list filters all, by default empty list
    last_generation_only : bool
        Whether or not to only filter last generation, by default False

    Attributes
    ----------
    use_unpredicted : bool
        Whether or not use to unpredicted reactions, by default False
    use_wildcards : Bool
        Whether or not to allow wildcard reactions through, by default True
    generation_list : list
        Generations to apply filter -- empty list filters all, by default empty list
    last_generation_only : bool
        Whether or not to only filter last generation, by default False
    """

    def __init__(
        self,
        use_unpredicted: bool = False,
        use_wildcards: bool = True,
        generation_list: List[int] = [],
        last_generation_only: bool = False,
    ) -> None:
        # Import modules that are only required for the RFC
        self._filter_name = "Reaction Feasibility Filter"
        self.use_unpredicted = use_unpredicted
        self.generation_list = generation_list
        self.last_generation_only = last_generation_only

        # Threshold for feasibility
        self._feas_threshold = 0.32

    @property
    def filter_name(self) -> str:
        return self._filter_name

    def _pre_print(self) -> None:
        """Print before filtering."""
        print(f"Filter out reactions for feasibility.")

    def _post_print(
        self, pickaxe: Pickaxe, n_total: int, n_filtered: int, time_sample: float
    ) -> None:
        """Print after filtering."""
        print(
            (
                f"{n_filtered} of {n_total} "
                "compounds selected after feasibility filtering in time_sample."
            )
        )

    def _get_inputs(self, rxn_list, pickaxe):
        """rxn_list will be pickaxe eventually"""
        # Get reactions information
        def get_cpd_smiles(cpd_id):
            return pickaxe.compounds[cpd_id]["SMILES"]

        reactions_info = {}
        for rxn_id in rxn_list:
            rxn = pickaxe.reactions[rxn_id]
            reactants = [
                get_cpd_smiles(v[1]) for v in rxn["Reactants"] if v[1].startswith("C")
            ]
            products = [
                get_cpd_smiles(v[1]) for v in rxn["Products"] if v[1].startswith("C")
            ]

            pairs = product(reactants, products)
            reactions_info[rxn["_id"]] = list(pairs)

        # Process this information
        input_info = {}
        input_fails = {}
        for rxn_id, reaction_pairs in reactions_info.items():
            if not reaction_pairs:
                continue

            for i, (reactant_smiles, product_smiles) in enumerate(reaction_pairs):
                if len(reactant_smiles) <= 120:
                    if len(product_smiles) <= 120:
                        mol1 = MolFromSmiles(reactant_smiles)
                        mol2 = MolFromSmiles(product_smiles)
                        mol1 = RemoveHs(mol1)
                        mol2 = RemoveHs(mol2)
                        reactant_smiles = MolToSmiles(mol1)
                        product_smiles = MolToSmiles(mol2)

                        # TODO what does this fix? from original code
                        if "M" in reactant_smiles or "M" in product_smiles:
                            input_fails[rxn_id + "_" + str(i)] = None
                        else:
                            input_info[rxn_id + "_" + str(i)] = [
                                reactant_smiles,
                                product_smiles,
                            ]
                    else:
                        input_fails[rxn_id + "_" + str(i)] = None
                else:
                    input_fails[rxn_id + "_" + str(i)] = None

        return input_info, input_fails

    def _choose_items_to_filter(self, pickaxe: Pickaxe, processes: int = 1) -> Set[str]:
        """
        Check the compounds against the MW constraints and return
        compounds to filter.
        """

        def get_cpd_smiles(cpd_id):
            return pickaxe.compounds[cpd_id]["SMILES"]

        cpds_remove_set = set()
        rxns_remove_set = set()

        # TODO put these statements together
        # No reactions to filter for
        if len(pickaxe.reactions) == 0:
            print("No reactions to calculate feasibilty for.")
            return cpds_remove_set, rxns_remove_set

        if self.last_generation_only and pickaxe.generation != self.generation:
            print("Not filtering for this generation using feasibility.")
            return cpds_remove_set, rxns_remove_set

        if self.generation_list and (self.generation - 1) not in self.generation_list:
            print("Not filtering for this generation using feasibility.")
            return cpds_remove_set, rxns_remove_set

        reactions_to_check = []
        for cpd in pickaxe.compounds.values():
            # Compounds are in generation and correct type
            if cpd["Generation"] == pickaxe.generation and cpd["Type"] not in [
                "Coreactant",
                "Target Compound",
            ]:
                reactions_to_check.extend(cpd["Product_of"])

        reactions_to_check = set(reactions_to_check)

        rfc_input, rfc_fails = self._get_inputs(reactions_to_check, pickaxe)

        feasibility_dict = {}
        rfc_predictions = _get_feasibility(rfc_input, self._feas_threshold)

        for rxn_id, feas_scores in rfc_predictions.items():
            rxn_id = rxn_id.split("_")[0]
            # Set to true if any are feasible, set to false if doesn't
            if feas_scores[2] == "feasible":
                feasibility_dict[rxn_id] = True
            elif rxn_id not in feasibility_dict:
                feasibility_dict[rxn_id] = False

        # Assign values for unpredicted
        for rxn_id in reactions_to_check:
            rxn = pickaxe.reactions[rxn_id]
            use_rxn = False
            if rxn_id not in feasibility_dict:
                use_rxn = self.use_unpredicted

                # Check for asterisks
                reactants = [
                    get_cpd_smiles(v[1])
                    for v in rxn["Reactants"]
                    if v[1].startswith("C")
                ]
                products = [
                    get_cpd_smiles(v[1])
                    for v in rxn["Products"]
                    if v[1].startswith("C")
                ]

                compounds = reactants.extend(products)
                if any([("*" in cpd) for cpd in compounds]) and self.use_wildcards:
                    use_rxn = True

            feasibility_dict[rxn_id] = use_rxn

        rxns_remove_set = set(
            rxn_id for rxn_id, keep in feasibility_dict.items() if not keep
        )

        return cpds_remove_set, rxns_remove_set


def _get_feasibility(input_info, feas_threshold=0.32):
    # inputs
    # Load Model

    weight_file = cwd.parent / "data/feasibility/final_model.h5"
    vae_model_file = cwd.parent / "data/feasibility/vae_model.pth"
    model_file = cwd.parent / "data/feasibility/final_model.json"

    json_file = open(model_file, "r")
    loaded_model_json = json_file.read()
    json_file.close()

    loaded_model = model_from_json(loaded_model_json)
    loaded_model.load_weights(weight_file)

    vae_model = MolecularVAE()
    vae_model.load_state_dict(torch.load(vae_model_file, map_location=device))
    vae_model.to(device)
    vae_model.eval()

    feature_info = _calculate_features(input_info, vae_model)
    results = _predict_reaction_feasibility(feature_info, loaded_model, feas_threshold)

    return results


def _predict_reaction_feasibility(feature_info, model, feas_threshold):
    results = {}
    for each_label in feature_info:
        if all([type(a) == type(None) for a in feature_info[each_label]]):
            continue

        X1 = feature_info[each_label][0]
        X2 = feature_info[each_label][1]

        val_list = []
        for i in range(10):
            pred_result = model.predict([X1, X2])
            val = pred_result[0][0]
            val_list.append(val)
        mean_val = np.mean(val_list)
        std_val = np.std(val_list)
        final_val = mean_val - (std_val * 0.5)
        if final_val >= feas_threshold:
            feasibility = "feasible"
        else:
            feasibility = "infeasible"

        results[each_label] = [mean_val, std_val, feasibility]

    return results


def _calculate_features(input_info, vae_model):
    r = 2
    feature_info = {}
    for each_label in input_info:

        reactant_smi = input_info[each_label][0]
        product_smi = input_info[each_label][1]

        if "9" in reactant_smi or "9" in product_smi:
            # Does do large rings
            feature_info[each_label] = [None, None]
            continue

        z1 = _calc_z(vae_model, reactant_smi)
        z1 = z1[0]

        z2 = _calc_z(vae_model, product_smi)
        z2 = z2[0]

        reactant_mol = MolFromSmiles(reactant_smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(reactant_mol, radius=r, nBits=1024)
        bits = fp.ToBitString()
        reactant_feature = []
        for f in z1:
            reactant_feature.append(float(f))
        for f in bits:
            reactant_feature.append(int(f))

        product_mol = MolFromSmiles(product_smi)
        fp = AllChem.GetMorganFingerprintAsBitVect(product_mol, radius=r, nBits=1024)
        bits = fp.ToBitString()
        product_feature = []
        for f in z2:
            product_feature.append(float(f))
        for f in bits:
            product_feature.append(int(f))

        reactant_feature = [reactant_feature]
        product_feature = [product_feature]
        reactant_feature = np.asarray(reactant_feature)
        product_feature = np.asarray(product_feature)

        feature_info[each_label] = [reactant_feature, product_feature]
    return feature_info


def _calc_z(model, smi):
    start = smi
    start = start.ljust(120)
    oh = OneHotFeaturizer()
    start_vec = torch.from_numpy(oh.featurize([start]).astype(np.float32)).to(device)
    z = model(start_vec)[1].cpu().detach().numpy()
    return z


class MolecularVAE(nn.Module):
    def __init__(self):
        super(MolecularVAE, self).__init__()

        self.conv1d1 = nn.Conv1d(120, 9, kernel_size=9)
        self.conv1d2 = nn.Conv1d(9, 9, kernel_size=9)
        self.conv1d3 = nn.Conv1d(9, 10, kernel_size=11)
        self.fc0 = nn.Linear(90, 435)
        self.fc11 = nn.Linear(435, 292)
        self.fc12 = nn.Linear(435, 292)

        self.fc2 = nn.Linear(292, 292)
        self.gru = nn.GRU(292, 501, 3, batch_first=True)
        self.fc3 = nn.Linear(501, 35)

    def encode(self, x):
        h = F.relu(self.conv1d1(x))
        h = F.relu(self.conv1d2(h))
        h = F.relu(self.conv1d3(h))
        h = h.view(h.size(0), -1)
        h = F.selu(self.fc0(h))
        return self.fc11(h), self.fc12(h)

    def reparametrize(self, mu, logvar):
        if self.training:
            std = torch.exp(0.5 * logvar)
            eps = 1e-2 * torch.randn_like(std)
            w = eps.mul(std).add_(mu)
            return w
        else:
            return mu

    def decode(self, z):
        z = F.selu(self.fc2(z))
        z = z.view(z.size(0), 1, z.size(-1)).repeat(1, 120, 1)
        out, h = self.gru(z)
        out_reshape = out.contiguous().view(-1, out.size(-1))
        y0 = F.softmax(self.fc3(out_reshape), dim=1)
        y = y0.contiguous().view(out.size(0), -1, y0.size(-1))
        return y

    def forward(self, x):
        mu, logvar = self.encode(x)
        z = self.reparametrize(mu, logvar)
        return self.decode(z), mu, logvar

    def sample(self, x, decode_attempts=1000, noise_norm=1.0):
        target_atoms = ["C", "N", "O", "S", "F", "Cl", "Br", "H"]

        oh = OneHotFeaturizer()
        mu, logvar = self.encode(x)
        z = self.reparametrize(mu, logvar)
        z2 = z.cpu().detach().numpy()
        smiles_results = {}
        for i in range(decode_attempts):
            Z = self.perturb_z(z2, noise_norm, False)
            z3 = torch.from_numpy(Z.astype(np.float32)).to("cuda")
            X = self.decode(z3)
            X = X.cpu().detach().numpy()
            y = np.argmax(X, axis=2)
            smiles = oh.decode_smiles_from_index(y[0])
            mol = MolFromSmiles(smiles)
            if mol != None:
                smiles_results[smiles] = mol

        for i in range(decode_attempts):
            Z = self.perturb_z(z2, noise_norm, True)
            z3 = torch.from_numpy(Z.astype(np.float32)).to("cuda")
            X = self.decode(z3)
            X = X.cpu().detach().numpy()
            y = np.argmax(X, axis=2)
            smiles = oh.decode_smiles_from_index(y[0])
            mol = MolFromSmiles(smiles)
            if mol != None:
                smiles_results[smiles] = mol

        for smiles in smiles_results:
            flag = True
            each_tmp_mol = smiles_results[smiles]
            for each_atom in each_tmp_mol.GetAtoms():
                if each_atom.GetSymbol() not in target_atoms:
                    flag = False
            if flag:
                print(smiles)
        return

    def perturb_z(self, z, noise_norm, constant_norm=False):
        if noise_norm > 0.0:
            noise_vec = np.random.normal(0, 1, size=z.shape)
            noise_vec = noise_vec / np.linalg.norm(noise_vec)
            if constant_norm:
                return z + (noise_norm * noise_vec)
            else:
                noise_amp = np.random.uniform(0, noise_norm, size=(z.shape[0], 1))
                return z + (noise_amp * noise_vec)
        else:
            return z


if __name__ == "__main__":
    pass
