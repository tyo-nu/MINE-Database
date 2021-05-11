"""Generate rules to use in pickaxe runs."""
from io import StringIO
from pathlib import Path
from typing import Tuple, List

from bisect import bisect_right

import pandas as pd


pwd = Path(__file__).parent

pattern_ignore_dictionary = {
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
    ignore_containing: List[str] = None,
    **kwargs,
) -> Tuple[StringIO, StringIO, str]:
    """Generate generalize metacyc rule subsets.

    Generate subsets of the metacyc generalized reaction opreators by specifying the
    number of rules of the fraction coverage of metacyc desired. Rules are chosen
    in the order of rules that map the most reactions to least. For fractional coverage
    the lowest number of rules that give a coverage less than or equal to the specified
    coverage is given.

    Specific rules can be ignored as well to prune the reaction list to specific
    reactions.

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
    ignore_containing: List[str], optional
        A list containing features to ignore. Valid features are:
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

    Returns
    -------
    Tuple[StringIO, StringIO, str]
        A tuple containing two streams that contain the reaction rule information and
        the rule name.
    """
    # Whether or not to return counts
    return_counts = kwargs.get("return_counts", False)
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
    rule_df = pd.merge(rule_counts, rule_df, on="Name")

    # Filter out any reactions based on filtering
    name_append = ""
    if anaerobic:
        rule_df = rule_df[~rule_df["Reactants"].str.contains("^O2|;O2|O2;")]
        name_append += "_anaerobic"

    if ignore_containing:
        patterns_to_ignore = r"|".join(
            pattern_ignore_dictionary[feature] for feature in ignore_containing
        )
        rule_df = rule_df[~rule_df["SMARTS"].str.contains(patterns_to_ignore)]
        name_append += "_ignore_containing"

    # Reindex DF
    rule_df = rule_df.reset_index(drop=True)

    # Get regex pattern to search for rules to filter dataframe
    if n_rules:
        # Default to n_rules if all are specified
        # pat = "|".join(list(rule_df.iloc[0:n_rules]["rule"]))
        rule_name = f"Metacyc_generalized_{n_rules}_rules"
    elif fraction_coverage:
        # Calculate CDF
        rule_df["cdf"] = rule_df["counts"].cumsum() / sum(rule_counts["counts"])
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

    stream = StringIO()
    new_rules.to_csv(stream, sep="\t", index=False)
    stream = StringIO(stream.getvalue())

    return stream, coreactants, rule_name


def metacyc_intermediate(
    n_rules: int = None,
    fraction_coverage: float = None,
    anaerobic: float = False,
    ignore_containing: List[str] = None,
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
    ignore_containing: List[str], optional
        A list containing features to ignore. Valid features are:
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
    uniprot_rules = pwd / Path(
        "data/metacyc_rules/metacyc_intermediate_rules_uniprot.tsv"
    )
    coreactants = pwd / Path("data/metacyc_rules/metacyc_coreactants.tsv")

    # Get intermediate dataframe and sort by reactions mapped
    rule_df = pd.read_csv(rules, delimiter="\t")
    rule_df["counts"] = rule_df.Comments.map(lambda s: len(s.split(";")))
    total_rxns = rule_df.counts.sum()

    # Filter DF to only desired compounds using general rules pattern matching
    general_rule_stream, _, _ = metacyc_generalized(
        n_rules=None,
        fraction_coverage=1,
        anaerobic=anaerobic,
        ignore_containing=ignore_containing,
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
    missing_rules = set(general_rule_df["Name"].values).difference(
        set(v.split("_")[0] for v in rule_df["Name"].values)
    )
    missing_df = general_rule_df[
        general_rule_df["Name"].str.contains("|".join(missing_rules))
    ]
    rule_df = rule_df.append(missing_df)

    # Generate
    # CDF for determining fraction coverage
    rule_df = rule_df.sort_values(by="counts", ascending=False)

    # Change Comments to be uniprot
    uniprot_df = pd.read_csv(
        uniprot_rules, delimiter="\t", usecols=["Name", "Comments"]
    )
    uniprot_df = uniprot_df.rename(columns={"Comments": "Uniprot"})
    # rule_df = rule_df.drop(columns=["Comments"])
    rule_df = pd.merge(rule_df, uniprot_df, on="Name", how="left")

    # Filter out any reactions based on filtering
    name_append = ""
    if anaerobic:
        name_append += "_anaerobic"

    if ignore_containing:
        name_append += "_ignore_containing"

    # Get regex pattern to search for rules to filter dataframe
    if n_rules:
        # Default to n_rules if all are specified
        rule_name = f"Metacyc_intermediate_{n_rules}_rules"
    elif fraction_coverage:
        # Calculate CDF
        rule_df["cdf"] = rule_df["counts"].cumsum() / total_rxns
        n_rules = bisect_right(rule_df["cdf"].values, fraction_coverage) + 1
        rule_name = (
            f"Metacyc_intermediate_{fraction_coverage}_fraction_coverage".replace(
                ".", ","
            )
        )
    else:
        # pat = "|".join(list(rule_counts["rule"]))
        n_rules = len(rule_df)
        rule_name = "Metacyc_generalized"

    # get stream to new rules
    rule_name += name_append

    # get stream to new rules
    new_rules = rule_df.iloc[0:n_rules]
    new_rules = new_rules.drop(columns=["counts", "cdf"])
    stream = StringIO()
    new_rules.to_csv(stream, sep="\t", index=False)
    stream = StringIO(stream.getvalue())

    return stream, coreactants, rule_name


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
