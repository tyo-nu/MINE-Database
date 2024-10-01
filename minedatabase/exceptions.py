"""Submodule providing exceptions used in the mine database package."""
from rdkit.Chem import Mol
from rdkit.Chem import MolToSmiles


class DisconnectedMoleculeException(ValueError):
    """Exception raised when a molecule is disconnected."""
    def __init__(self, molecule: Mol):
        molecule_smile = MolToSmiles(molecule)
        super().__init__((
            f"The provided molecule with SMILES '{molecule_smile}' is "
            "disconnected. This is most common when compounds are salts. "
            "If you want to include disconnected molecules, you can set "
            "the `fragmented_mols` parameter to True."
        ))