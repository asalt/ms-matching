# driver.py
# incomplete
import itertools
import math
from ms_matching.read_config import parse_msfragger_config, extract_mod_definitions
from ms_matching.db import open_db, create_schema, store_mod_definitions, load_mod_definitions
from ms_matching.digest import (
    digest_protein, generate_mod_variants,
    build_aa_mass_from_mod_def_dict,
    generate_theoretical_fragments,
    calculate_precursor_mass
)
from ms_matching.io import read_fasta
from ms_matching.match import do_match
from ms_matching.db import store_peptide_fragments
import warnings

def initialize_database(config_dict: dict, db_path: str):
    mod_defs = extract_mod_definitions(config_dict)
    # import pdb; pdb.set_trace()

    with open_db(db_path) as conn:
        create_schema(conn)
        store_mod_definitions(mod_defs, conn)

def index_peptides_from_fasta(fasta_path: str, db_path: str, bin_width: float = 0.5, precursor_charge_range = (2,6), pH=2, ion_types=("b", "y")):
    from pyteomics import electrochem
    with open_db(db_path) as conn:
        mod_defs = load_mod_definitions(conn)
        aa_mass, mod_aliases = build_aa_mass_from_mod_def_dict(mod_defs)


    for header, sequence in read_fasta(fasta_path):
        for peptide in digest_protein(sequence):
            max_charge = math.ceil(electrochem.charge(peptide, pH=pH)) 
            charge_states = [int(x) for x in range(precursor_charge_range[0], max(max_charge,2)+1)]

            parsed_peptides = generate_mod_variants(peptide, mod_aliases)

            for parsed_peptide, charge in itertools.product(parsed_peptides, charge_states):
                precursor_mass = calculate_precursor_mass(parsed_peptide, aa_mass, charge=charge)
                fragments = generate_theoretical_fragments(
                    parsed_peptide=parsed_peptide,
                    aa_mass=aa_mass,
                    bin_width=bin_width,
                    maxcharge=max(charge-1, 1),
                    ion_types=ion_types,
                    #conn=conn
                    # hdf5_file will likely need to be passed here too
                )
                yield {"fragments" : fragments,
                        "peptide" : parsed_peptide,
                        "aa_mass" : aa_mass,
                        "mod_aliases" : mod_aliases,
                        "bin_width": bin_width,
                        "charge": charge,
                        }


def match_spectrum(spectrum_data: dict, precursor_mz: float, charge: int, db_path: str, bin_width: float):
    from ms_matching.db import query_candidates_by_mz, load_fragment_bins
    with open_db(db_path) as conn:
        candidates = query_candidates_by_mz(conn, precursor_mz, charge, tolerance_da=1.0)
        matches = []
        for peptide, mods, hdf5_path in candidates:
            theo_bins, frag_types = load_fragment_bins(hdf5_path)  # TODO: refactor load_fragment_bins to accept conn/hdf5
            score = do_match(peptide, spectrum_data, theo_bins=theo_bins)
            matches.append((peptide, mods, score))
        return matches


def run(config_path:str = None, db_path: str = None, fasta_path: str=None, config_dict=None, ion_types=None, precursor_charge_range=None):

    if config_dict is None:
        config_dict = parse_msfragger_config(config_path)
        initialize_database(config_dict, db_path)

    if ion_types is None:
        ion_types = config_dict.get("fragment_ion_series", ("b", "y"))
    if precursor_charge_range is None:
        precursor_charge_range_str = config_dict.get("precursor charge", "2 6")
        try:
            precursor_charge_range = [ int(x) for x in precursor_charge_range_str.split(" ") ]
        except ValueError:
            warnings.warn(f"invalid precursor charge {precursor_charge_range_str}, falling back to default")
            precursor_charge_range = [2, 6]

    fragment_data = index_peptides_from_fasta(fasta_path, db_path, precursor_charge_range=precursor_charge_range)
    yield from index_peptides_from_fasta(fasta_path, db_path, precursor_charge_range=precursor_charge_range)


def run_and_store(config_path: str, db_path: str, fasta_path: str):

    initialize_database(config_dict, db_path)

    config_dict = parse_msfragger_config(config_path)
    ion_types = config_dict.get("fragment_ion_series", ("b", "y"))
    precursor_charge_range_str = config_dict.get("precursor charge", "2 6")
    try:
        precursor_charge_range = [ int(x) for x in precursor_charge_range_str.split(" ") ]
    except ValueError:
        warnings.warn(f"invalid precursor charge {precursor_charge_range_str}, falling back to default")
        precursor_charge_range = [2, 6]

    fragment_data = run(db_path=db_path, fasta_path=fasta_path, config_dict=config_dict, ion_types=ion_types, precursor_charge_range=precursor_charge_range)

    with open_db(db_path) as conn:
        for frag in fragment_data:
            store_peptide_fragments(
                peptide = ''.join(parsed_peptide), #: str,
                # mods: str,
                # charge: int,
                # mass: float,
                # binned_fragments: List[int],
                # frag_types: List[str],
                #ordinals: List[int],
                hdf5_file=None, #,: h5py.File,
                sqlite_conn=conn,
            )



def store(fragpep_ix_gen): # does nothing at the moment
    for res in fragpep_ix_gen:
        pass

        # store_peptide_fragments(
        #     peptide = ''.join(parsed_peptide), #: str,
        #     mods: str,
        #     charge: int,
        #     mass: float,
        #     binned_fragments: List[int],
        #     frag_types: List[str],
        #     ordinals: List[int],
        #     hdf5_file: h5py.File,
        #     sqlite_conn: sqlite3.Connection,
        # )
