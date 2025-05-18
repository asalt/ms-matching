# digest.py
from itertools import combinations, product
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

from pyteomics import mass
PROTON_MASS = mass.nist_mass.get("H+", {0: (1.00727646677, 1.0), 1: (1.00727646677, 1.0)})[0][0]


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


# def build_aa_mass_from_mod_def_dict(
#     mod_def_dict: Dict[str, Dict[str, float]]
# ) -> Tuple[Dict[str, float], Dict[str, List[str]]]:
#     from pyteomics import mass
# 
#     SPECIAL_TARGETS = {
#         "n": "peptide-N-term",
#         #"^": "N-term",
#         "[^": "protein-N-term",
#         "c": "peptide-C-term",
#         "]": "protein-C-term",
#     }
# 
#     aa_mass = dict(mass.std_aa_mass)
#     mod_aliases = {}
# 
#     for tag, entry in mod_def_dict.items():
#         delta = entry["mass"]
#         for raw_aa in entry["targets"]:
#             aa = SPECIAL_TARGETS.get(raw_aa, raw_aa)
# 
#             if aa in mass.std_aa_mass:
#                 mod_tag = f"{tag}{aa}"
#                 aa_mass[mod_tag] = aa_mass[aa] + delta
#                 mod_aliases.setdefault(aa, []).append(mod_tag)
# 
#             elif aa in ("N-term", "C-term"):
#                 # You may decide not to add to aa_mass, just alias
#                 mod_aliases.setdefault(aa, []).append(tag)
#             else:
#                 print(f"⚠️ Warning: Unknown mod target '{raw_aa}' → skipped")
# 
#     return aa_mass, mod_aliases

# def build_aa_mass_from_mod_def_dict(mod_def_dict):
#     from pyteomics import mass
# 
#     aa_mass = dict(mass.std_aa_mass)
#     mod_aliases = {}
# 
#     for tag, entry in mod_def_dict.items():
#         delta = entry["mass"]
#         residues = entry.get("residues", [])
#         position = entry.get("position")
# 
#         for aa in residues:
#             mod_tag = f"{tag}{aa}"
#             aa_mass[mod_tag] = aa_mass[aa] + delta
#             mod_aliases.setdefault(aa, []).append(mod_tag)
# 
#         # 🔧 Add this block for terminal mods
#         if position:
#             mod_aliases.setdefault(position, []).append(tag)
# 
#     return aa_mass, mod_aliases

def build_aa_mass_from_mod_def_dict_orig(mod_def_dict: Dict,
        aa_mass:Dict[str, float] | None = None,):
    from pyteomics import mass as _mass
    aa_mass = aa_mass or dict(_mass.std_aa_mass)
    mod_aliases = {}

    for tag, entry in mod_def_dict.items():
        delta = entry["mass"]
        residues = entry.get("residues", [])
        position = entry.get("position")

        # 1. Residue-specific mods
        for aa in residues:
            mod_tag = f"{tag}{aa}"
            aa_mass[mod_tag] = aa_mass[aa] + delta
            mod_aliases.setdefault(aa, []).append(mod_tag)

        # 2. Positional mods — generate modAA entries for ALL amino acids
        if position:
            for aa in aa_mass:
                mod_tag = f"{tag}{aa}"
                aa_mass[mod_tag] = aa_mass[aa] + delta
            mod_aliases.setdefault(position, []).append(tag)

    return aa_mass, mod_aliases

def build_aa_mass_from_mod_def_dict(
    mod_def_dict: Dict[str, Dict],
    aa_mass: Dict[str, float] | None = None,
) -> Tuple[Dict[str, float], Dict[str, list]]:
    """
    • Adds every (mod-tag + residue) mass to `aa_mass`
    • Builds `mod_aliases` so that:
        ─ residue → [modAA, …]
        ─ 'protein-N-term' / 'peptide-N-term' / … → [terminal-tag, …]
    The function never mutates a dictionary it is iterating over.
    """
    from pyteomics import mass as _mass

    aa_mass = dict(aa_mass or _mass.std_aa_mass)      # make a private copy
    canonical_aas = [aa for aa in aa_mass if len(aa) == 1 and aa.isupper()]

    mod_aliases: Dict[str, list] = {}

    for tag, entry in mod_def_dict.items():
        delta     = entry["mass"]
        residues  = entry.get("residues", [])
        position  = entry.get("position")              # None | 'protein-N-term' …

        # ── 1. Residue-specific adjustments ────────────────────────────

        # ── 2. Positional mods (N- / C-term) ───────────────────────────
        if position:
            #import pdb; pdb.set_trace()
            # use *snapshot* of canonical residues → no mutation during iter
            for aa in canonical_aas:
                mod_tag = f"{tag}"#{aa}"
                # skip if we already created it in the residue loop
                if mod_tag not in aa_mass:
                    aa_mass[mod_tag] = delta # aa_mass[aa] + delta
            mod_aliases.setdefault(position, []).append(tag)
        else:
            for aa in residues:
                mod_tag = f"{tag}{aa}"
                aa_mass[mod_tag] = aa_mass[aa] + delta
                mod_aliases.setdefault(aa, []).append(mod_tag)

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


def apply_terminal_mods( # check if this is used
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



def generate_mod_variants_orig(
    peptide: str,
    mod_aliases: Dict[str, List[str]],
    max_mods_per_peptide: int = 3,
    include_terminal_mods: bool = True,
) -> List[List[str]]:
    base = list(peptide)
    variants = set()

    # Step 1: Internal mods
    modifiable = [
        (i, aa, mod_aliases[aa])
        for i, aa in enumerate(base)
        if aa in mod_aliases
    ]

    for n in range(max_mods_per_peptide + 1):
        for positions in combinations(modifiable, n):
            mod_choices = [
                [(i, mod) for mod in mods]
                for (i, aa, mods) in positions
            ]
            for mods_applied in product(*mod_choices):
                modded = base.copy()
                for i, mod_tag in mods_applied:
                    modded[i] = mod_tag
                variants.add(tuple(modded))

    # Step 2: Terminal mods
    if include_terminal_mods:
        extended = set()
        for seq in variants or {tuple(base)}:
            first = seq[0]
            last = seq[-1]
            first_aa = first[-1]  # Assumes final char is AA
            last_aa = last[-1]

            nterm_mods = mod_aliases.get("protein-N-term", [])
            cterm_mods = mod_aliases.get("protein-C-term", [])
            residue_mods = mod_aliases.get(first_aa, [])

            # Generate all N-term combinations
            for nmod in nterm_mods:
                # N-term only
                extended.add((f"{nmod}{first}", *seq[1:]))

                # If base AA not already modified, try stacked mods
                if first == first_aa:
                    for rmod in residue_mods:
                        stacked = f"{nmod}{rmod}{first_aa}"
                        extended.add((stacked, *seq[1:]))

            # C-term only
            for cmod in cterm_mods:
                extended.add((*seq[:-1], f"{cmod}{last}"))

            # Both N and C-term
            for nmod in nterm_mods:
                for cmod in cterm_mods:
                    ntag = f"{nmod}{first}"
                    ctag = f"{cmod}{last}"
                    extended.add((ntag, *seq[1:-1], ctag))

                    if first == first_aa:
                        for rmod in residue_mods:
                            stacked = f"{nmod}{rmod}{first_aa}"
                            extended.add((stacked, *seq[1:-1], ctag))

            # Always keep original
            extended.add(seq)

        variants = extended

    ret = [list(seq) for seq in variants]
    # import pdb; pdb.set_trace()
    return ret


def _stack(new_mod: str, token: str) -> str: # not using this
    """
    Prepend `new_mod` to `token`.
    If `token` is already modified (i.e. its length > 1), insert the ProForma
    ‘&’ separator so Pyteomics can parse the stacked modifications.

    >>> _stack("mod42010", "S")
    'mod42010S'
    >>> _stack("mod42010-", "mod79966S")
    'mod42010-mod79966S'
    """
    return f"{new_mod}{'&' if len(token) > 1 else ''}{token}"


def generate_mod_variants(
    peptide: str,
    mod_aliases: Dict[str, List[str]],
    max_mods_per_peptide: int = 3,
    include_terminal_mods: bool = True,
) -> List[List[str]]:
    """

    """
    base: List[str] = list(peptide)
    variants: Set[Tuple[str, ...]] = set()

    # ------------------------------------------------------------------
    # 1. internal residue mods
    # ------------------------------------------------------------------

    variants.add(tuple(peptide))
    modifiable = [
        (i, aa, mod_aliases[aa])
        for i, aa in enumerate(base)
        if aa in mod_aliases
    ]

    for n in range(max_mods_per_peptide + 1):
        for positions in combinations(modifiable, n):
            mod_choices = [
                [(i, mod) for mod in mods]  # one choice per residue
                for (i, aa, mods) in positions
            ]
            for mods_applied in product(*mod_choices):
                modded = base.copy()
                for i, mod_tag in mods_applied:
                    modded[i] = mod_tag
                variants.add(tuple(modded))

    # ------------------------------------------------------------------
    # 2. N- and C-terminal mods (with optional stacking)
    # ------------------------------------------------------------------
    if include_terminal_mods:
        extended: Set[Tuple[str, ...]] = set()
        for seq in variants or {tuple(base)}:
            first, last = seq[0], seq[-1]
            first_aa, last_aa = first[-1], last[-1]

            nterm_mods = mod_aliases.get("protein-N-term", [])
            cterm_mods = mod_aliases.get("protein-C-term", [])
            residue_mods = mod_aliases.get(first_aa, [])

            # N-term only
            # import pdb; pdb.set_trace()
            for nmod in nterm_mods:
                extended.add((nmod, *seq))
                #extended.add((_stack(nmod, first), *seq[1:]))

                # stacked (N-term + residue) if base AA isn’t already modified
                if first == first_aa:
                    # import pdb; pdb.set_trace()
                    for rmod in residue_mods:
                        stacked = (f"{nmod}", f"{rmod}")
                        extended.add((*stacked, *seq[1:]))

            # C-term only
            
            if last == last_aa: # not sure this is right
                for cmod in cterm_mods:
                    extended.add((*seq, cmod))
                    #extended.add((*seq[:-1], _stack(cmod, last)))

            # Both N- and C-term.  # this is impossible, unless the sequence only has 1 aa
            # import pdb; pdb.set_trace()
            # for nmod in nterm_mods:
            #     for cmod in cterm_mods:
            #         #ntag = _stack(nmod, first)
            #         #ctag = _stack(cmod, last)
            #         extended.add((nmod, *seq, cmod))

            #         if first == first_aa:
            #             import pdb; pdb.set_trace()
            #             for rmod in residue_mods:
            #                 stacked = f"{nmod}&{rmod}{first_aa}" # check this
            #             stacked = (f"{nmod}", f"{first_aa}{rmod}")
            #                 extended.add((ntag, *seq, ctag))

            # # keep the original un-modified sequence
            # extended.add(seq)

        variants = variants | extended

    return [list(seq) for seq in variants]




# def generate_mod_variants(
#     peptide: str,
#     mod_aliases: Dict[str, List[str]],
#     max_mods_per_peptide: int = 3,
#     include_terminal_mods: bool = True,
# ) -> List[List[str]]:
#     """
#     Generate all variants of a peptide with modifications.
# 
#     Args:
#         peptide: base sequence, e.g. "MAST"
#         mod_aliases: dict like {'M': ['oxM'], 'S': ['phS'], 'N-term': ['ac']}
#         max_mods_per_peptide: max internal mods
#         include_terminal_mods: whether to try N-term/C-term variants
# 
#     Returns:
#         List of parsed peptide sequences with modifications
#     """
#     base = list(peptide)
#     variants = set()
# 
#     modifiable_positions = [
#         (i, mods)
#         for i, aa in enumerate(base)
#         if aa in mod_aliases
#         for mods in mod_aliases[aa]
#     ]
# 
#     for n in range(0, max_mods_per_peptide + 1):
#         for pos_combo in combinations(modifiable_positions, n):
#             modded = list(base)
#             for i, mod_tag in pos_combo:
#                 modded[i] = mod_tag
# 
#             parsed = modded
# 
#             if include_terminal_mods:
#                 # try combinations of N-term and C-term
#                 terminal_variants = [
#                     apply_terminal_mods(parsed, mod_aliases, {}),
#                     (
#                         apply_terminal_mods(
#                             parsed, mod_aliases, {"N-term": mod_aliases["N-term"][0]}
#                         )
#                         if "N-term" in mod_aliases
#                         else parsed
#                     ),
#                     (
#                         apply_terminal_mods(
#                             parsed, mod_aliases, {"C-term": mod_aliases["C-term"][0]}
#                         )
#                         if "C-term" in mod_aliases
#                         else parsed
#                     ),
#                     (
#                         apply_terminal_mods(
#                             parsed,
#                             mod_aliases,
#                             {
#                                 "N-term": mod_aliases["N-term"][0],
#                                 "C-term": mod_aliases["C-term"][0],
#                             },
#                         )
#                         if "N-term" in mod_aliases and "C-term" in mod_aliases
#                         else parsed
#                     ),
#                 ]
#                 for p in terminal_variants:
#                     variants.add(tuple(p))  # hashable
#             else:
#                 variants.add(tuple(parsed))
# 
#     return [parser.parse("".join(seq)) for seq in variants]


# def generate_peptides(protein="AAAQQQQSMETK", **kwargs):
#     # under construction
# 
#     proteins = kwargs.get("proteins", ["AAAQQQQSMETK", "AAAQQQCPALK"])
#     max_modis = kwargs.get("max_modis", 3)
#     config = kwargs.get(
#         "config",
#         {
#             "variable_mod_01": "15.9949 M 3",
#             # 'variable_mod_02': '42.0106 [^ 1',
#             "variable_mod_03": "79.96633 STY 3",
#         },
#     )
#     import pdb; pdb.set_trace()
#     aa_mass, mod_aliases = build_aa_mass_from_variable_mods(
#         config, tag_map=default_tag_map
#     )
#     peptides = list()  # we will not do this later
#     # for header, aa_seq in proteins:
#     digest = digest_protein(protein)  # TODO add more args here
#     for peptide in digest:
#         peptide_modis = generate_mod_variants(
#             peptide, mod_aliases, max_mods_per_peptide=max_modis
#         )
#         peptides.append(peptide)
#         peptides += peptide_modis  # this part is wrong, probably
#     return peptides, aa_mass

# def generate_peptides(protein="AAAQQQQSMETK", config=None, **kwargs):
#     from ms_matching.read_config import extract_mod_definitions
#     from ms_matching.digest import build_aa_mass_from_mod_def_dict
# 
#     if config is None:
#         config = kwargs.get("config", {
#             "variable_mod_01": "15.9949 M",
#             "variable_mod_03": "79.96633 STY"
#         })
# 
#     mod_defs = extract_mod_definitions(config)
#     aa_mass, mod_aliases = build_aa_mass_from_mod_def_dict(mod_defs)
# 
#     max_modis = kwargs.get("max_modis", 3)
#     peptides = []
# 
#     for peptide in digest_protein(protein):
#         peptides.append(peptide)
#         mod_variants = generate_mod_variants(peptide, mod_aliases, max_mods_per_peptide=max_modis)
#         peptides += mod_variants
# 
#     return peptides, aa_mass


def generate_theoretical_fragments(parsed_peptide, aa_mass, bin_width=0.5, maxcharge=2, **kwargs):

    frags = make_fragments(parsed_peptide, aa_mass=aa_mass, maxcharge=maxcharge, **kwargs)  # generator
    # TODO keep track of b/y
    frags = mz_to_bin(np.array(list(frags)), bin_width)
    return frags


def prepare_matching_environment(conn):
    mod_defs = load_mod_definitions(conn)
    aa_mass, mod_aliases = build_aa_mass_from_mod_def_dict(mod_defs)
    return aa_mass, mod_aliases



# if __name__ == "__main__":
#     peptides, aa_mass = generate_peptides(protein="MELTINGQQQQK", max_modis=2)
#     for p in peptides:
#         print("Parsed:", p)
#         print("Mass:", mass.fast_mass(p, aa_mass=aa_mass))
#         print("-" * 40)



def calculate_precursor_mass(peptide, aa_mass, charge):
    from pyteomics.mass import fast_mass
    mass = fast_mass(peptide, aa_mass=aa_mass) + (charge * PROTON_MASS)
    mz = mass / charge
    return mz


# we can also use pyteomics.mass.fast_mass to perform all of these calculations for us
def make_fragments(
    peptide, ion_types=("b", "y"), maxcharge=1, aa_mass=mass.std_aa_mass
):  # from https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    from pyteomics import parser, mass
    # import pdb; pdb.set_trace()
    for i in range(1, len(peptide)):
        for ion_type in ion_types:
            for charge in range(1, maxcharge + 1):
                if ion_type[0] in "abc":
                    yield mass.fast_mass2(
                        peptide[:i], ion_type=ion_type, charge=charge, aa_mass=aa_mass
                    )
                else:
                    yield mass.fast_mass2(
                        peptide[i:], ion_type=ion_type, charge=charge, aa_mass=aa_mass
                    )

