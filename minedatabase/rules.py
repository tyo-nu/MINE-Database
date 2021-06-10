"""Generate rules to use in pickaxe runs."""
from bisect import bisect_right
from io import StringIO
from os import remove
from pathlib import Path
from typing import List, Tuple

import pandas as pd


pwd = Path(__file__).parent

pattern_dictionary = {
    "aromatic": r":\[|\]:",
    "aromatic_oxygen": r"^\[#6:\d+\]:|\[#6:\d+\]:|\[#6:\d+\]\d:|\[#6:\d+\]:\d",
    "carbonyl": r"=\[#8:\d\]|\[#8:\d\]=",
    "nitrogen": r"\[#7:\d+\]",
    "oxygen": r"\[#8:\d+\]",
    "fluorine": r"\[#9:\d+\]",
    "phosphorus": r"\[#15:\d+\]",
    "sulfur": r"\[#16:\d+\]",
    "chlorine": r"\[#17:\d+\]",
    "bromine": r"\[#35:\d+\]",
    "iodine": r"\[#53:\d+\]",
    "halogen": r"\[#(9|17|35|53):\d+\]",
}


def metacyc_generalized(
    n_rules: int = None,
    fraction_coverage: float = None,
    anaerobic: float = False,
    include_containing: List[str] = None,
    exclude_containing: List[str] = None,
    **kwargs,
) -> Tuple[StringIO, StringIO, str]:
    """Generate generalize metacyc rule subsets.

    Generate subsets of the metacyc generalized reaction opreators by specifying the
    number of rules of the fraction coverage of metacyc desired. Rules are chosen
    in the order of rules that map the most reactions to least. For fractional coverage
    the lowest number of rules that give a coverage less than or equal to the specified
    coverage is given.

    Specific groups can be specified to be excluded or used as well to prune to specific
    rules. This is a two step process:
        1) Select rules by include_containing. If none specified, use all rules.
        2) Remove rules by excluded_containing.

    Parameters
    ----------
    n_rules : int, optional
        Number of rules to use. If excluded rules result in less than specified
        number then all rules are taken, by default None.
    fraction_coverage : float, optional
        The fraction of coverage desired. This may be impossible to reach depending
        on which rules are excluded, by default None.
    anaerobic: float, optional
        Whether to remove oxygen requiring reactions.
    include_containing: List[str], optional
        A list containing features to include. Valid features are:
            - aromatic
            - aromatic_oxygen
            - carbonyl
            - halogen
            - nitrogen
            - oxygen
            - phosphorus
            - sulfur
            - fluorine
            - chlorine
            - bromine
            - iodine
        sending None gives all groups, by default None.
    exclude_containing: List[str], optional
        A list containing features to exclude.
            - aromatic
            - aromatic_oxygen
            - carbonyl
            - halogen
            - nitrogen
            - oxygen
            - phosphorus
            - sulfur
            - fluorine
            - chlorine
            - bromine
            - iodine
        By default None.

    Returns
    -------
    Tuple[StringIO, StringIO, str]
        A tuple containing two streams (reaction rules and cofactor) and
        the rule name.
    """
    # Whether or not to return counts
    return_counts = kwargs.get("return_counts", False)
    return_all = kwargs.get("return_all", False)
    # Metacyc Rules
    reaction_mapping = pwd / Path("data/metacyc_rules/metacyc_mapped.tsv")
    rules = pwd / Path("data/metacyc_rules/metacyc_generalized_rules.tsv")
    coreactants = pwd / Path("data/metacyc_rules/metacyc_coreactants.tsv")

    # Get dataframe to help select rules
    rule_df = pd.read_csv(reaction_mapping, delimiter="\t", usecols=["rule"]).rename(
        columns={"rule": "Name"}
    )

    # Generate CDF for determining fraction coverage
    rule_df.Name = rule_df.Name.map(lambda s: s.split("_")[0])
    rule_counts = rule_df.value_counts().rename_axis("Name").reset_index(name="counts")

    # Attach rules to this DF
    rule_df = pd.read_csv(
        rules,
        delimiter="\t"
        # usecols=["Name", "Reactants", "SMARTS", "Prod"]
    )
    removed_df = pd.DataFrame()

    # Filter out any reactions based on filtering
    name_append = ""
    removed_ids = set()
    if anaerobic:
        removed_ids.update(
            rule_df[rule_df["Reactants"].str.contains("^O2|;O2|O2;")].index.to_list()
        )
        # rule_df = rule_df[~rule_df["Reactants"].str.contains("^O2|;O2|O2;")]
        name_append += "_anaerobic"

    if include_containing:
        patterns_to_include = r"|".join(
            pattern_dictionary[feature] for feature in include_containing
        )
        removed_ids.update(
            rule_df[
                ~rule_df["SMARTS"].str.contains(patterns_to_include)
            ].index.to_list()
        )
        # rule_df = rule_df[rule_df["SMARTS"].str.contains(patterns_to_include)]
        name_append += "_with_inclusion"

    if exclude_containing:
        patterns_to_exclude = r"|".join(
            pattern_dictionary[feature] for feature in exclude_containing
        )
        removed_ids.update(
            rule_df[rule_df["SMARTS"].str.contains(patterns_to_exclude)].index.to_list()
        )
        # rule_df = rule_df[~rule_df["SMARTS"].str.contains(patterns_to_exclude)]
        name_append += "_with_exclusion"

    # Reorder dataframe
    rule_ids = rule_df.index.to_list()
    keep_ids = [i for i in rule_ids if i not in removed_ids]
    removed_ids = [i for i in rule_ids if i in removed_ids]
    new_ids = keep_ids + removed_ids

    # Merge with rule counts, saving old ids
    rule_df["index"] = rule_df.index
    rule_df = pd.merge(rule_counts, rule_df, on="Name").set_index("index")

    # TODO why is this happening?
    # Calculate ids to keep, some rules are dropped due to not having reactions
    new_ids = [i for i in new_ids if i in rule_df.index]
    rule_df = rule_df.loc[new_ids]

    # Calculate CDF and drop unwanted rules
    rule_df["cdf"] = rule_df["counts"].cumsum() / rule_df["counts"].sum()

    if return_all:
        rule_df.loc[rule_df.index.isin(keep_ids), "type"] = "keep"
        rule_df.loc[~rule_df.index.isin(keep_ids), "type"] = "drop"
    else:
        rule_df = rule_df[rule_df.index.isin(keep_ids)]
        rule_df = rule_df.reset_index(drop=False)

    # Get regex pattern to search for rules to filter dataframe
    if n_rules:
        # Default to n_rules if all are specified
        # pat = "|".join(list(rule_df.iloc[0:n_rules]["rule"]))
        rule_name = f"Metacyc_generalized_{n_rules}_rules"
    elif fraction_coverage:
        # Calculate CDF
        n_rules = bisect_right(rule_df["cdf"].values, fraction_coverage) + 1
        # n_rules = rule_df["cdf"].sub(fraction_coverage).abs().idxmin() + 1
        # pat = "|".join(list(rule_df.iloc[0:n_rules]["rule"]))
        rule_name = (
            f"Metacyc_generalized_{fraction_coverage}_fraction_coverage".replace(
                ".", ","
            )
        )
    else:
        # pat = "|".join(list(rule_counts["rule"]))
        n_rules = len(rule_df)
        rule_name = "Metacyc_generalized"

    rule_name += name_append

    # get stream to new rules
    new_rules = rule_df.iloc[0:n_rules]
    if not return_counts:
        new_rules = new_rules.drop(columns=["counts"])
        if "cdf" in new_rules.columns:
            new_rules = new_rules.drop(columns=["cdf"])

    rule_stream = StringIO()
    new_rules.to_csv(rule_stream, sep="\t", index=False)
    rule_stream = StringIO(rule_stream.getvalue())

    return rule_stream, coreactants, rule_name


def metacyc_intermediate(
    n_rules: int = None,
    fraction_coverage: float = None,
    anaerobic: float = False,
    include_containing: List[str] = None,
    exclude_containing: List[str] = None,
) -> Tuple[StringIO, StringIO, str]:
    """Generate intermediate metacyc rule subsets.

    Generate subsets of the metacyc intermediate reaction opreators by specifying the
    number of rules of the fraction coverage of metacyc desired. Coverage and number of
    rules are taken from the generalized operators and their intermediate operators are
    chosen.

    Parameters
    ----------
    n_rules : int, optional
        Number of rules to use, by default None.
    fraction_coverage : float, optional
        The fraction of coverage desired, by default None.
    anaerobic: float, optional
        Whether to remove oxygen requiring reactions.
    include_containing: List[str], optional
        A list containing features to include. Valid features are:
            - aromatic
            - aromatic_oxygen
            - carbonyl
            - halogen
            - nitrogen
            - oxygen
            - phosphorus
            - sulfur
            - fluorine
            - chlorine
            - bromine
            - iodine
        sending None gives all groups, by default None.
    exclude_containing: List[str], optional
        A list containing features to exclude. Valid features are:
            - aromatic
            - aromatic_oxygen
            - carbonyl
            - halogen
            - nitrogen
            - oxygen
            - phosphorus
            - sulfur
            - fluorine
            - chlorine
            - bromine
            - iodine
        by default None

    Returns
    -------
    Tuple[StringIO, StringIO, str]
        A tuple containing two streams that contain the reaction rule information and
        the rule name.
    """
    # Metacyc Rules
    rules = pwd / Path("data/metacyc_rules/metacyc_intermediate_rules.tsv")
    # uniprot_rules = pwd / Path(
    #     "data/metacyc_rules/metacyc_intermediate_rules_uniprot.tsv"
    # )
    coreactants = pwd / Path("data/metacyc_rules/metacyc_coreactants.tsv")

    # Get intermediate dataframe and sort by reactions mapped
    rule_df = pd.read_csv(rules, delimiter="\t")
    # rule_df["counts"] = rule_df.Comments.map(lambda s: len(s.split(";")))
    total_rxns = rule_df.counts.sum()

    # Filter DF to only desired compounds using general rules pattern matching
    general_rule_stream, _, _ = metacyc_generalized(
        n_rules=None,
        fraction_coverage=1,
        anaerobic=anaerobic,
        include_containing=include_containing,
        exclude_containing=exclude_containing,
        return_counts=True,
    )
    general_rule_df = pd.read_csv(
        general_rule_stream,
        delimiter="\t",
        usecols=["Name", "Reactants", "SMARTS", "Products", "counts", "Comments"],
    )
    valid_rules_pattern = "|".join(rule for rule in general_rule_df["Name"].values)
    rule_df = rule_df[rule_df["Name"].str.contains(valid_rules_pattern)]

    # Some generalized have no intermediate, bring those in
    # missing_rules = set(general_rule_df.Name) - set(
    #     v.split("_")[0] for v in rule_df["Name"].values
    # )
    # if missing_rules:
    #     missing_df = general_rule_df[
    #         general_rule_df["Name"].str.contains("|".join(missing_rules))
    #     ]
    #     rule_df = rule_df.append(missing_df)

    # Generate
    # CDF for determining fraction coverage
    rule_df = rule_df.sort_values(by="counts", ascending=False)
    rule_df["cdf"] = rule_df["counts"].cumsum() / total_rxns

    # saved file to do this instead
    # # Change Comments to be uniprot
    # uniprot_df = pd.read_csv(
    #     uniprot_rules, delimiter="\t", usecols=["Name", "Comments"]
    # )
    # uniprot_df = uniprot_df.rename(columns={"Comments": "Uniprot"})
    # # rule_df = rule_df.drop(columns=["Comments"])
    # rule_df = pd.merge(rule_df, uniprot_df, on="Name", how="left")

    # Filter out any reactions based on filtering
    name_append = ""
    if anaerobic:
        name_append += "_anaerobic"

    if include_containing:
        name_append += "_with_inclusion"

    if exclude_containing:
        name_append += "_with_exclusion"

    # Get regex pattern to search for rules to filter dataframe
    if n_rules:
        # Default to n_rules if all are specified
        rule_name = f"Metacyc_intermediate_{n_rules}_rules"
    elif fraction_coverage:
        # Calculate CDF
        n_rules = bisect_right(rule_df["cdf"].values, fraction_coverage) + 1
        rule_name = (
            f"Metacyc_intermediate_{fraction_coverage}_fraction_coverage".replace(
                ".", ","
            )
        )
    else:
        # pat = "|".join(list(rule_counts["rule"]))
        n_rules = len(rule_df)
        rule_name = "Metacyc_intermediate"

    # get stream to new rules
    rule_name += name_append

    # get stream to new rules
    new_rules = rule_df.iloc[0:n_rules]
    if "cdf" in new_rules.columns:
        new_rules = new_rules.drop(columns=["cdf"])
    if "counts" in new_rules.columns:
        new_rules = new_rules.drop(columns=["counts"])
    rule_stream = StringIO()
    new_rules.to_csv(rule_stream, sep="\t", index=False)
    rule_stream = StringIO(rule_stream.getvalue())

    return rule_stream, coreactants, rule_name


def BNICE() -> Tuple[Path, Path, str]:
    """Generate BNICE rules.

    Generate the original BNICE rules that were use before the improved MetaCyc rules
    were generated.

    Returns
    -------
    Tuple[Path, Path, str]
        The path to the rules and coreactants and the rule name.
    """
    rules = pwd / Path("data/original_rules/EnzymaticReactionRules.tsv")
    coreactants = pwd / Path("data/original_rules/EnzymaticCoreactants.tsv")
    rule_name = "BNICE"

    return rules, coreactants, rule_name


import seaborn as sns


sns.scatterplot()
