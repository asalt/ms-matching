# digest.py
from itertools import combinations
from pyteomics import parser
from pyteomics import mass
from ms_matching.io import read_fasta
from ms_matching import match

import numpy as np

from typing import Dict, Tuple, List

from ms_matching.db import (
    open_db,
    store_mod_definitions,
    create_schema,
    load_mod_definitions,
)


default_tag_map = {
    15.995: "ox",
    79.966: "ph",
    42.011: "ac",
}
default_aa_mass = dict(mass.std_aa_mass)


def mz_to_bin(mz_array, bin_width):  # TODO let it take ppm
    # TODO add check if not np.array
    return np.round((mz_array - bin_width / 2) / bin_width).astype(int)


def build_aa_mass_from_variable_mods(
    config, tag_map=default_tag_map, aa_mass=default_aa_mass
):
    """
    Parses MSFragger-style variable_mod entries into a pyteomics-compatible aa_mass dict.

    Args:
        config: dict from parse_msfragger_config()
        tag_map: optional map from rounded delta mass to prefix tag (e.g. 15.995 → 'ox')

    Returns:
        aa_mass: modified dict for pyteomics
        mod_aliases: mapping from base AA to modified tag name
    """
    # aa_mass = dict(mass.std_aa_mass)
    mod_aliases = {}

    keys = [x for x in config.keys() if x.startswith("variable_mod_")]
    for key in keys:
        value = config[key]
        parts = value.split()
        if len(parts) < 2:
            continue
        delta_mass = float(parts[0])
        targets = parts[1]
        rounded_mass = round(delta_mass, 3)

        tag_prefix = tag_map.get(rounded_mass, f"mod{int(abs(delta_mass) * 1000)}")
        for aa in targets:
            if aa in mass.std_aa_mass:
                tag = f"{tag_prefix}{aa}"
                aa_mass[tag] = mass.std_aa_mass[aa] + delta_mass
                mod_aliases.setdefault(aa, []).append(tag)
            else:
                pass  # could log/handle N-term/C-term cases separately later
    return aa_mass, mod_aliases


# def build_aa_mass_from_mod_def_dict(mod_def_dict: Dict[str, Dict[str, float]]) -> Tuple[Dict[str, float], Dict[str, str]]:
#     """
#     Rebuilds a pyteomics-compatible aa_mass dictionary and mod_aliases
#     from mod_definitions loaded from SQLite.
#
#     Args:
#         mod_def_dict: e.g. {
#             'ox': {'mass': 15.9949, 'targets': 'M'},
#             'ph': {'mass': 79.9663, 'targets': 'STY'}
#         }
#
#     Returns:
#         aa_mass: e.g. { 'oxM': 147.0354, ... }
#         mod_aliases: e.g. { 'M': 'oxM', 'S': 'phS' }
#     """
#     aa_mass = dict(mass.std_aa_mass)
#     mod_aliases = {}
#
#     for tag, entry in mod_def_dict.items():
#         delta = entry['mass']
#         for aa in entry['targets']:
#             mod_tag = f"{tag}{aa}"
#             aa_mass[mod_tag] = aa_mass[aa] + delta
#             mod_aliases.setdefault(aa, []).append(mod_tag)
#
#     return aa_mass, mod_aliases


def build_aa_mass_from_mod_def_dict(
    mod_def_dict: Dict[str, Dict[str, float]]
) -> Tuple[Dict[str, float], Dict[str, List[str]]]:
    from pyteomics import mass

    SPECIAL_TARGETS = {
        "n": "N-term",
        "^": "N-term",
        "[": "N-term",
        "c": "C-term",
        "]": "C-term",
    }

    aa_mass = dict(mass.std_aa_mass)
    mod_aliases = {}

    for tag, entry in mod_def_dict.items():
        delta = entry["mass"]
        for raw_aa in entry["targets"]:
            aa = SPECIAL_TARGETS.get(raw_aa, raw_aa)

            if aa in mass.std_aa_mass:
                mod_tag = f"{tag}{aa}"
                aa_mass[mod_tag] = aa_mass[aa] + delta
                mod_aliases.setdefault(aa, []).append(mod_tag)

            elif aa in ("N-term", "C-term"):
                # You may decide not to add to aa_mass, just alias
                mod_aliases.setdefault(aa, []).append(tag)
            else:
                print(f"⚠️ Warning: Unknown mod target '{raw_aa}' → skipped")

    return aa_mass, mod_aliases


def digest_protein(
    protein_sequence: str,
    enzyme: str = "trypsin/P",
    missed_cleavages: int = 2,
    min_length: int = 7,
    max_length: int = 50,
) -> list:
    """
    Digest a protein sequence into peptides using the specified enzyme.

    Parameters
    ----------
    protein_sequence : str
        The protein sequence to digest.
    enzyme : str, optional
        The enzyme to use for digestion. Default is 'trypsin/P'.
    missed_cleavages : int, optional
        The number of missed cleavages allowed. Default is 2.
    min_length : int, optional
        The minimum length of the peptides. Default is 7.
    max_length : int, optional
        The maximum length of the peptides. Default is 50.

    Returns
    -------
    list
        A list of digested peptides.
    """
    # Create a digest object
    if enzyme == "trypsin/P":
        enzyme = "[KR]"

    peptides = parser.cleave(
        protein_sequence,
        rule=enzyme,
        missed_cleavages=missed_cleavages,
        min_length=min_length,
        max_length=max_length,
    )

    # assert len(peptides) > 0, "No peptides were generated from the protein sequence."

    return peptides


# def generate_mod_variants(
#     peptide: str, mod_aliases: dict, max_mods_per_peptide: int = 3
# ):
#     """
#     Generate all variants of a peptide with up to N modifications.
#
#     Args:
#         peptide: e.g. "STYMQL"
#         mod_aliases: e.g. {'S': ['phS'], 'T': ['phT'], 'M': ['oxM']}
#         max_mods_per_peptide: maximum allowed modifications of each type
#
#     Returns:
#         List of parsed peptide sequences suitable for fast_mass2
#     """
#     variants = set()
#
#     # Get all modifiable positions and the mod tag to apply
#     mod_sites = [
#         (i, mod_aliases[aa][0]) for i, aa in enumerate(peptide) if aa in mod_aliases
#     ]
#
#     for n in range(1, max_mods_per_peptide + 1):
#         for mod_pos_combo in combinations(mod_sites, n):
#             mod_seq = list(peptide)
#             for i, mod_tag in mod_pos_combo:
#                 mod_seq[i] = mod_tag
#             variants.add(tuple(mod_seq))  # Use tuple to deduplicate
#     return [parser.parse("".join(seq)) for seq in variants]


def apply_terminal_mods(
    parsed_peptide: List[str],
    mod_aliases: Dict[str, List[str]],
    selected_mods: Dict[str, str] = None,
) -> List[str]:
    """
    Apply terminal modifications to a parsed peptide sequence.

    Args:
        parsed_peptide: e.g. ['M', 'A', 'S']
        mod_aliases: e.g. { 'N-term': ['ac'], 'C-term': ['modC'] }
        selected_mods: optional dict like {'N-term': 'ac', 'C-term': 'modC'}

    Returns:
        Modified peptide: e.g. ['ac', 'M', 'A', 'S', 'modC']
    """
    peptide = list(parsed_peptide)  # make a copy

    if selected_mods:
        if "N-term" in selected_mods:
            peptide.insert(0, selected_mods["N-term"])
        if "C-term" in selected_mods:
            peptide.append(selected_mods["C-term"])

    else:
        # Auto-apply the first available terminal mod if defined
        if "N-term" in mod_aliases:
            peptide.insert(0, mod_aliases["N-term"][0])
        if "C-term" in mod_aliases:
            peptide.append(mod_aliases["C-term"][0])

    return peptide


def generate_mod_variants(
    peptide: str,
    mod_aliases: Dict[str, List[str]],
    max_mods_per_peptide: int = 3,
    include_terminal_mods: bool = True,
) -> List[List[str]]:
    """
    Generate all variants of a peptide with modifications.

    Args:
        peptide: base sequence, e.g. "MAST"
        mod_aliases: dict like {'M': ['oxM'], 'S': ['phS'], 'N-term': ['ac']}
        max_mods_per_peptide: max internal mods
        include_terminal_mods: whether to try N-term/C-term variants

    Returns:
        List of parsed peptide sequences with modifications
    """
    base = list(peptide)
    variants = set()

    modifiable_positions = [
        (i, mods)
        for i, aa in enumerate(base)
        if aa in mod_aliases
        for mods in mod_aliases[aa]
    ]

    for n in range(0, max_mods_per_peptide + 1):
        for pos_combo in combinations(modifiable_positions, n):
            modded = list(base)
            for i, mod_tag in pos_combo:
                modded[i] = mod_tag

            parsed = modded

            if include_terminal_mods:
                # try combinations of N-term and C-term
                terminal_variants = [
                    apply_terminal_mods(parsed, mod_aliases, {}),
                    (
                        apply_terminal_mods(
                            parsed, mod_aliases, {"N-term": mod_aliases["N-term"][0]}
                        )
                        if "N-term" in mod_aliases
                        else parsed
                    ),
                    (
                        apply_terminal_mods(
                            parsed, mod_aliases, {"C-term": mod_aliases["C-term"][0]}
                        )
                        if "C-term" in mod_aliases
                        else parsed
                    ),
                    (
                        apply_terminal_mods(
                            parsed,
                            mod_aliases,
                            {
                                "N-term": mod_aliases["N-term"][0],
                                "C-term": mod_aliases["C-term"][0],
                            },
                        )
                        if "N-term" in mod_aliases and "C-term" in mod_aliases
                        else parsed
                    ),
                ]
                for p in terminal_variants:
                    variants.add(tuple(p))  # hashable
            else:
                variants.add(tuple(parsed))

    return [parser.parse("".join(seq)) for seq in variants]


def generate_peptides(protein="AAAQQQQSMETK", **kwargs):
    # under construction

    proteins = kwargs.get("proteins", ["AAAQQQQSMETK", "AAAQQQCPALK"])
    max_modis = kwargs.get("max_modis", 3)
    config = kwargs.get(
        "config",
        {
            "variable_mod_01": "15.9949 M 3",
            # 'variable_mod_02': '42.0106 [^ 1',
            "variable_mod_03": "79.96633 STY 3",
        },
    )
    aa_mass, mod_aliases = build_aa_mass_from_variable_mods(
        config, tag_map=default_tag_map
    )
    peptides = list()  # we will not do this later
    # for header, aa_seq in proteins:
    digest = digest_protein(protein)  # TODO add more args here
    for peptide in digest:
        peptide_modis = generate_mod_variants(
            peptide, mod_aliases, max_mods_per_peptide=max_modis
        )
        peptides.append(peptide)
        peptides += peptide_modis  # this part is wrong, probably
    return peptides, aa_mass


def generate_theoretical_fragments(parsed_peptide, aa_mass, bin_width=0.5, conn=None):

    frags = match.make_fragments(parsed_peptide, aa_mass=aa_mass)  # generator
    # TODO keep track of b/y
    frags = mz_to_bin(np.array(list(frags)), bin_width)
    return frags


def prepare_matching_environment(conn):
    mod_defs = load_mod_definitions(conn)
    aa_mass, mod_aliases = build_aa_mass_from_mod_def_dict(mod_defs)
    return aa_mass, mod_aliases
