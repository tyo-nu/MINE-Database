"""Generate rules to use in pickaxe runs."""
from io import StringIO
from pathlib import Path
from typing import Tuple

import pandas as pd


pwd = Path(__file__).parent


def metacyc_generalized(
    n_rules: int = None, fraction_coverage: float = None
) -> Tuple[StringIO, StringIO, str]:
    """Generate generalize metacyc rule subsets.

    Generate subsets of the metacyc generalized reaction opreators by specifying the
    number of rules of the fraction coverage of metacyc desired. Rules are chosen
    in the order of rules that map the most reactions to least. Ties will be broken by
    lowest rule number.

    Parameters
    ----------
    n_rules : int, optional
        Number of rules to use, by default None.
    fraction_coverage : float, optional
        The fraction of coverage desired, by default None.

    Returns
    -------
    Tuple[StringIO, StringIO, str]
        A tuple containing two streams that contain the reaction rule information and
        the rule name.
    """
    # Metacyc Rules
    reaction_mapping = pwd / Path("data/metacyc_rules/metacyc_mapped.tsv")
    rules = pwd / Path("data/metacyc_rules/metacyc_generalized_rules.tsv")
    coreactants = pwd / Path("data/metacyc_rules/metacyc_coreactants.tsv")

    # Get dataframe to help select rules
    headers = ["rule", "source", "SMARTS", "Map"]
    rule_df = pd.read_csv(
        reaction_mapping, delimiter="\t", names=headers, usecols=["rule"]
    )

    # Generate CDF for determining fraction coverage
    rule_df.rule = rule_df.rule.map(lambda s: s.split("_")[0])
    rule_counts = rule_df.value_counts().rename_axis("rule").reset_index(name="counts")
    rule_counts["rule_id"] = rule_counts.index
    rule_counts["cdf"] = rule_counts["counts"].cumsum() / sum(rule_counts["counts"])

    # Get regex pattern to search for rules to filter dataframe
    if n_rules:
        # Default to n_rules if all are specified
        pat = "|".join(list(rule_counts.iloc[0:n_rules]["rule"]))
        rule_name = f"Metacyc_generalized_{n_rules}_rules"
    elif fraction_coverage:
        n_rules = rule_counts["cdf"].sub(fraction_coverage).abs().idxmin() + 1
        pat = "|".join(list(rule_counts.iloc[0:n_rules]["rule"]))
        rule_name = (
            f"Metacyc_generalized_{fraction_coverage}_fraction_coverage".replace(
                ".", ","
            )
        )
    else:
        pat = "|".join(list(rule_counts["rule"]))
        rule_name = "Metacyc_generalized"

    # get stream to new rules
    full_rule_df = pd.read_csv(rules, delimiter="\t")
    new_rules = full_rule_df[full_rule_df["Name"].str.contains(pat)]
    stream = StringIO()
    new_rules.to_csv(stream, sep="\t")
    stream = StringIO(stream.getvalue())

    return stream, coreactants, rule_name


def metacyc_intermediate(
    n_rules: int = None, fraction_coverage: float = None
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

    Returns
    -------
    Tuple[StringIO, StringIO, str]
        A tuple containing two streams that contain the reaction rule information and
        the rule name.
    """
    # Metacyc Rules
    reaction_mapping = pwd / Path("data/metacyc_rules/metacyc_mapped.tsv")
    rules = pwd / Path("data/metacyc_rules/metacyc_intermediate_rules.tsv")
    coreactants = pwd / Path("data/metacyc_rules/metacyc_coreactants.tsv")

    # Get dataframe to help select rules
    headers = ["rule", "source", "SMARTS", "Map"]
    rule_df = pd.read_csv(
        reaction_mapping, delimiter="\t", names=headers, usecols=["rule"]
    )

    # Generate CDF for determining fraction coverage
    rule_df.rule = rule_df.rule.map(lambda s: s.split("_")[0])
    rule_counts = rule_df.value_counts().rename_axis("rule").reset_index(name="counts")
    rule_counts["rule_id"] = rule_counts.index
    rule_counts["cdf"] = rule_counts["counts"].cumsum() / sum(rule_counts["counts"])

    # Get regex pattern to search for rules to filter dataframe
    if n_rules:
        # Default to n_rules if all are specified
        pat = "|".join(list(rule_counts.iloc[0:n_rules]["rule"]))
        rule_name = f"Metacyc_intermediate_{n_rules}_rules"
    elif fraction_coverage:
        n_rules = rule_counts["cdf"].sub(fraction_coverage).abs().idxmin() + 1
        pat = "|".join(list(rule_counts.iloc[0:n_rules]["rule"]))
        rule_name = (
            f"Metacyc_intermediate_{fraction_coverage}_fraction_coverage".replace(
                ".", ","
            )
        )
    else:
        pat = "|".join(list(rule_counts["rule"]))
        rule_name = "Metacyc_intermediate"

    # get stream to new rules
    full_rule_df = pd.read_csv(rules, delimiter="\t")
    new_rules = full_rule_df[full_rule_df["Name"].str.contains(pat)]
    stream = StringIO()
    new_rules.to_csv(stream, sep="\t")
    stream = StringIO(stream.getvalue())

    return stream, coreactants, rule_name


def metacyc_intermediate_uniprot(
    n_rules: int = None, fraction_coverage: float = None
) -> Tuple[StringIO, StringIO, str]:
    """Generate intermediate metacyc rule subsets with uniprot mapping.

    Generate a subset of the metacyc intermediate reaction operators with uniprot
    coverage by specifying the number of rules of the fraction coverage of metacyc
    desired. Coverage and number of rules are taken from the generalized operators and
    their intermediate operators are chosen.

    Only rules with uniprot mappings are returned.

    Parameters
    ----------
    n_rules : int, optional
        Number of rules to use, by default None.
    fraction_coverage : float, optional
        The fraction of coverage desired, by default None.

    Returns
    -------
    Tuple[StringIO, StringIO, str]
        A tuple containing two streams that contain the reaction rule information and
        the rule name.
    """
    # Metacyc Rules
    reaction_mapping = pwd / Path("data/metacyc_rules/metacyc_mapped.tsv")
    rules = pwd / Path("data/metacyc_rules/metacyc_intermediate_rules_uniprot.tsv")
    coreactants = pwd / Path("data/metacyc_rules/metacyc_coreactants.tsv")

    # Get dataframe to help select rules
    headers = ["rule", "source", "SMARTS", "Map"]
    rule_df = pd.read_csv(
        reaction_mapping, delimiter="\t", names=headers, usecols=["rule"]
    )

    # Generate CDF for determining fraction coverage
    rule_df.rule = rule_df.rule.map(lambda s: s.split("_")[0])
    # Only keep those with uniprot mapping
    full_rule_df = pd.read_csv(rules, delimiter="\t")
    full_rule_set = set(full_rule_df.Name.map(lambda s: s.split("_")[0]))
    rule_df = rule_df[rule_df["rule"].isin(full_rule_set)]
    # Calc CDF
    rule_counts = rule_df.value_counts().rename_axis("rule").reset_index(name="counts")
    rule_counts["rule_id"] = rule_counts.index
    rule_counts["cdf"] = rule_counts["counts"].cumsum() / sum(rule_counts["counts"])

    # Get regex pattern to search for rules to filter dataframe
    if n_rules:
        # Default to n_rules if all are specified
        pat = "|".join(list(rule_counts.iloc[0:n_rules]["rule"]))
        rule_name = f"Metacyc_intermediate_uniprot_{n_rules}_rules"
    elif fraction_coverage:
        n_rules = rule_counts["cdf"].sub(fraction_coverage).abs().idxmin() + 1
        pat = "|".join(list(rule_counts.iloc[0:n_rules]["rule"]))
        rule_name = (
            f"Metacyc_intermediate_uniprot_{fraction_coverage}"
            "_fraction_coverage".replace(".", ",")
        )
    else:
        pat = "|".join(list(rule_counts["rule"]))
        rule_name = "Metacyc_intermediate_uniprot"

    # get stream to new rules
    new_rules = full_rule_df[full_rule_df["Name"].str.contains(pat)]
    stream = StringIO()
    new_rules.to_csv(stream, sep="\t")
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
