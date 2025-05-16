# tests/unit/test_read_config.py
import pytest
from ms_matching.read_config import extract_mod_definitions

def test_single_variable_mod():
    config = {"variable_mod_01": "15.9949 M"}
    mods = extract_mod_definitions(config)
    assert "mod15994" in mods
    assert mods["mod15994"]["mass"] == 15.9949
    assert mods["mod15994"]["targets"] == "M"

def test_multiple_variable_mods():
    config = {
        "variable_mod_01": "15.9949 M",
        "variable_mod_02": "42.0106 K"
    }
    mods = extract_mod_definitions(config)
    assert len(mods) == 2
    assert "mod15994" in mods
    assert "mod42010" in mods

def test_custom_tag_map():
    config = {"variable_mod_01": "15.9949 M"}
    tag_map = {15.995: "Oxidation"}
    mods = extract_mod_definitions(config, tag_map)
    assert "Oxidation" in mods

def test_ignore_invalid_lines():
    config = {
        "variable_mod_01": "15.9949 M",
        "some_other_key": "foo bar",
        "variable_mod_02": "not a valid line"
    }
    with pytest.warns(UserWarning, match="Invalid mass value"):
        mods = extract_mod_definitions(config)
    assert "mod15994" in mods
    assert len(mods) == 1

