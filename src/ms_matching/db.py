# db.py
import sqlite3
from contextlib import contextmanager
import h5py
import numpy as np
from typing import List, Tuple, Optional, Dict
import json


@contextmanager
def open_db(db_path: str, read_only=False):
    """
    Context manager that opens a SQLite connection and ensures it is committed and closed.

    Usage:
        with open_db("my.db") as conn:
            conn.execute("...")
    """
    mode = "ro" if read_only else "rwc"
    uri = f"file:{db_path}?mode={mode}"
    conn = sqlite3.connect(uri, uri=True)
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise
    finally:
        conn.close()


def create_schema(sqlite_conn: sqlite3.Connection) -> None:
    sqlite_conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS peptide_index (
            id INTEGER PRIMARY KEY,
            peptide TEXT NOT NULL,
            mods TEXT,
            charge INTEGER NOT NULL,
            mass REAL NOT NULL,
            hdf5_path TEXT NOT NULL
        );
        CREATE INDEX IF NOT EXISTS idx_mass_charge ON peptide_index (mass, charge);
        CREATE INDEX IF NOT EXISTS idx_peptide_charge ON peptide_index (peptide, charge);
        CREATE INDEX IF NOT EXISTS idx_mods ON peptide_index (mods);
        """
    )

    sqlite_conn.executescript(
        """CREATE TABLE IF NOT EXISTS mod_definitions (
            tag TEXT PRIMARY KEY,
            mass REAL NOT NULL,
            targets TEXT NOT NULL,
            residues TEXT,
            position TEXT
        );
    """
    )


# def store_mod_definitions(
#     mod_def_dict: Dict[str, Dict[str, float]], sqlite_conn: sqlite3.Connection
# ) -> None:
#     """
#     Store mod definitions into the mod_definitions table.
# 
#     mod_def_dict: e.g. {
#         'ox': {'mass': 15.9949, 'targets': 'M'},
#         'ph': {'mass': 79.9663, 'targets': 'STY'}
#     }
#     """
#     rows = [
#         (tag, modinfo["mass"], modinfo["targets"])
#         for tag, modinfo in mod_def_dict.items()
#     ]
# 
#     # add logging here
#     sqlite_conn.executemany(
#         """
#         INSERT OR REPLACE INTO mod_definitions (tag, mass, targets)
#         VALUES (?, ?, ?)
#     """,
#         rows,
#     )

def store_mod_definitions(
    mod_def_dict: Dict[str, Dict[str, float]], sqlite_conn: sqlite3.Connection
) -> None:
    rows = [
        (
            tag,
            modinfo["mass"],
            modinfo.get("targets", ""),
            json.dumps(modinfo.get("residues", [])),
            modinfo.get("position"),
        )
        for tag, modinfo in mod_def_dict.items()
    ]

    sqlite_conn.executemany(
        """
        INSERT OR REPLACE INTO mod_definitions (tag, mass, targets, residues, position)
        VALUES (?, ?, ?, ?, ?)
        """,
        rows,
    )


# def load_mod_definitions(
#     sqlite_conn: sqlite3.Connection,
# ) -> Dict[str, Dict[str, float]]:
#     """
#     Load mod definitions from the database into a dictionary.
# 
#     Returns: e.g. {
#         'ox': {'mass': 15.9949, 'targets': 'M'},
#         'ph': {'mass': 79.9663, 'targets': 'STY'}
#     }
#     """
#     cur = sqlite_conn.execute("SELECT tag, mass, targets FROM mod_definitions")
#     return {
#         tag: {"mass": mass, "targets": targets} for tag, mass, targets in cur.fetchall()
#     }

def load_mod_definitions(sqlite_conn: sqlite3.Connection) -> Dict[str, Dict[str, float]]:
    cur = sqlite_conn.execute(
        "SELECT tag, mass, targets, residues, position FROM mod_definitions"
    )
    mod_defs = {}
    for tag, mass, targets, residues_json, position in cur.fetchall():
        mod_defs[tag] = {
            "mass": mass,
            "targets": targets,
            "residues": json.loads(residues_json) if residues_json else [],
            "position": position,
        }
    return mod_defs

def store_peptide_fragments(
    peptide: str,
    mods: str,
    charge: int,
    mass: float,
    binned_fragments: List[int],
    frag_types: List[str],
    ordinals: List[int],
    hdf5_file: h5py.File,
    sqlite_conn: sqlite3.Connection,
):
    group_name = f"/frags/{peptide}/z{charge}/{mods or 'none'}"

    grp = hdf5_file.require_group(group_name)
    grp.create_dataset("bins", data=np.array(binned_fragments, dtype=np.int32))
    grp.create_dataset("frag_types", data=np.array(frag_types, dtype="S1"))
    grp.create_dataset("ordinals", data=np.array(ordinals, dtype=np.int8))

    cur = sqlite_conn.cursor()
    cur.execute(
        """
        INSERT INTO peptide_index (peptide, mods, charge, mass, hdf5_path)
        VALUES (?, ?, ?, ?, ?)
    """,
        (peptide, mods, charge, mass, group_name),
    )
    sqlite_conn.commit()


def query_candidates_by_mz(
    sqlite_conn: sqlite3.Connection,
    precursor_mz: float,
    charge: int,
    tolerance_da: float = 1.0,
) -> List[Tuple[str, str, str]]:
    # Convert m/z to neutral mass
    mass = (precursor_mz * charge) - (charge * 1.007276)  # proton mass

    cur = sqlite_conn.cursor()
    cur.execute(
        """
        SELECT peptide, mods, hdf5_path FROM peptide_index
        WHERE charge = ? AND mass BETWEEN ? AND ?
    """,
        (charge, mass - tolerance_da, mass + tolerance_da),
    )

    return cur.fetchall()


def load_fragment_bins(
    hdf5_file: h5py.File, hdf5_path: str
) -> Tuple[np.ndarray, np.ndarray]:
    grp = hdf5_file[hdf5_path]
    bins = grp["bins"][...]
    frag_types = grp["frag_types"][...].astype(str)
    return bins, frag_types
