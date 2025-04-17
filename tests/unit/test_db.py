# test_db.py

from ms_matching.db import (
    open_db,
    create_schema,
    store_mod_definitions,
    load_mod_definitions,
)
from ms_matching.read_config import parse_msfragger_config, extract_mod_definitions


def test_store_and_load_mod_defs(tmp_path):
    db_path = tmp_path / "test.db"
    mod_defs = {
        "ox": {"mass": 15.9949, "targets": "M"},
        "ph": {"mass": 79.9663, "targets": "STY"},
    }

    with open_db(db_path) as conn:
        create_schema(conn)
        store_mod_definitions(mod_defs, conn)
        loaded = load_mod_definitions(conn)

    assert loaded == mod_defs
