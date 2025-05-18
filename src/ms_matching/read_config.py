import os
import re
import warnings
from typing import Dict, Tuple, List



def parse_msfragger_config(config_path):
    """
    Parses a simple MSFragger-style key=value config file into a dict.
    Ignores comments, blank lines, and parses numeric values when possible.
    """
    if not os.path.exists(config_path):
        raise ValueError(f"{config_path} does not exist")
    config = {}
    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.split("#", 1)[0].strip()

            # Try to convert to number
            if value.replace(".", "", 1).replace("-", "", 1).isdigit():
                value = float(value) if "." in value else int(value)
            elif "," in value:
                # Try parsing as list of numbers
                items = [v.strip() for v in value.split(",")]
                try:
                    value = [float(v) if "." in v else int(v) for v in items]
                except ValueError:
                    value = items  # fallback to string list
            config[key] = value
    return config




def extract_mod_definitions(config: Dict[str, str],
                            tag_map: Dict[float, str] = None,
                            aa_mass: Dict[str, float] = None) -> Dict[str, Dict]:

    from pyteomics import mass as _mass
    tag_map = tag_map or {}
    aa_mass = aa_mass or _mass.std_aa_mass
    ALL_AA: List[str] = sorted([aa for aa in aa_mass if len(aa) == 1 and aa.isupper()])


    # ------------------------------------------------------------------ #
    #  Tables                                                            #
    # ------------------------------------------------------------------ #
    # token ➜ (position, dash-direction)  -- dash is applied once later
    TERMINAL = {
        "[^": ("protein-N-term",  "right"),   # dash after tag
        "n":  ("peptide-N-term",  "right"),
        "]":  ("protein-C-term",  "left"),    # dash before tag
        "c":  ("peptide-C-term",  "left"),
    }

    # residue tokens are just single uppercase letters already in aa_mass

    TOKEN_REGEX = r"\[\^|n|c|\]|[A-Z]"

    # ------------------------------------------------------------------ #
    #  Inner parser                                                      #
    # ------------------------------------------------------------------ #
    def parse_entry(item: Tuple[str, str]):
        key, value = item
        if not key.startswith("variable_mod_"):
            return None

        parts = value.split()          # split on any whitespace
        if len(parts) < 2:
            warnings.warn(f"Ignoring malformed mod line: {key} = {value}")
            return None

        # ---------- mass & base tag ----------
        try:
            delta_mass = float(parts[0])
        except ValueError:
            warnings.warn(f"Invalid mass value in: {key} = {value}")
            return None

        rounded = round(delta_mass, 3)
        tag_base = tag_map.get(rounded, f"mod{int(abs(delta_mass*1000))}")

        # ---------- token scan ----------
        target_str = parts[1]
        tokens = re.findall(TOKEN_REGEX, target_str)

        residues: List[str] = []
        termini: set[str] = set()          # {'left'} or {'right'} or empty
        position: str | None = None

        for tok in tokens:
            # TERMINAL token?
            if tok in TERMINAL:
                pos, dash_dir = TERMINAL[tok]

                # conflict?  (e.g. both N- and C-term in the same string)
                if dash_dir in termini:
                    continue        # same token twice – harmless, ignore
                if termini:
                    raise ValueError(
                        f"Both N- and C-terminal tokens in '{target_str}'"
                    )

                termini.add(dash_dir)
                position = pos

                if tok == "[^":          # special “all residues” token
                    residues.extend(ALL_AA)

                continue    # token handled – go to next

            # Residue token?
            if tok in aa_mass:
                residues.append(tok)
                continue

            # Unknown token
            warnings.warn(
                f"Unknown mod target token '{tok}' in: {key} = {value}"
            )

        # ---------- post-processing ----------
        # unique residues, preserve order
        seen: set[str] = set()
        residues = [aa for aa in residues if not (aa in seen or seen.add(aa))]

        # add dash just once according to the collected direction, if any
        tag_final = {
            "left":  lambda t: f"-{t}",
            "right": lambda t: f"{t}-",
        }.get(next(iter(termini), None), lambda t: t)(tag_base)

        return tag_final, {
            "mass": delta_mass,
            "residues": residues,
            "position": position,
            "targets": target_str,   # keep raw string for downstream use
        }

    # ------------------------------------------------------------------ #
    #  Parse all lines & return dict                                     #
    # ------------------------------------------------------------------ #
    return {
        tag: info
        for parsed in map(parse_entry, config.items())
        if parsed is not None
        for tag, info in [parsed]
    }


    # def parse_entry(item: Tuple[str, str]):
    #     key, value = item
    #     if not key.startswith("variable_mod_"):
    #         return None

    #     parts = re.split(r"\s+", value.strip())
    #     if len(parts) < 2:
    #         warnings.warn(f"Ignoring malformed mod line: {key} = {value}")
    #         return None

    #     # ── mass ───────────────────────────────────────────────────────────
    #     try:
    #         delta_mass = float(parts[0])
    #     except ValueError:
    #         warnings.warn(f"Invalid mass value in mod line: {key} = {value}")
    #         return None
    #     rounded_mass = round(delta_mass, 3)
    #     tag = tag_map.get(rounded_mass, f"mod{int(abs(delta_mass * 1000))}")

    #     # ── target parsing ─────────────────────────────────────────────────
    #     target_str = parts[1]
    #     terminal_map = {
    #         '[^': "protein-N-term",
    #         'n':  "peptide-N-term",
    #         ']':  "protein-C-term",
    #         'c':  "peptide-C-term",
    #     }

    #     position = None
    #     residues: List[str] = []

    #     tag_modi_func = lambda x: x
    #     tokens = re.findall(r'\[\^|n|c|\]|[A-Z]', target_str)
    #     for tok in tokens:
    #         if tok in terminal_map:
    #             # set or overwrite the terminus position
    #             position = terminal_map[tok]

    #             # adjust tag for N- or C-terminal markers
    #             if "N-term" in position:
    #                 tag_modi_func = lambda x: x + "-"
    #             else:               # C-term
    #                 tag_modi_func = lambda x: "-" + x 

    #             # “[ˆ” means *every* residue
    #             if tok == '[^':
    #                 residues.extend(ALL_AA)

    #         elif tok in aa_mass:        # a canonical amino-acid code
    #             residues.append(tok)

    #         else:
    #             warnings.warn(f"Unknown mod target token '{tok}' "
    #                           f"in line: {key} = {value}")

    #     # deduplicate while preserving original order
    #     seen = set()
    #     residues = [aa for aa in residues if (aa not in seen and not seen.add(aa))]

    #     tag_final = tag_modi_func(tag)
    #     return tag_final, {
    #         "mass": delta_mass,
    #         "residues": residues,
    #         "position": position,
    #         "targets": target_str,
    #     }

    # # ----------------------------------------------------------------------
    return {
        tag: info
        for parsed in map(parse_entry, config.items())
        if parsed is not None
        for tag, info in [parsed]
    }

