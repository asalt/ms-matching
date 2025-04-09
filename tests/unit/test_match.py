
from ms_matching.match import do_match


def test_positive_match():
    mz_array = [148.08,
    978.51,
    277.12,
    849.47,
    391.16,
    735.43,
    504.25,
    622.34,
    607.25,
    519.33,
    735.35,
    391.24,
    848.43,
    278.15,
    979.47,
    147.11]
    intensity_array = [1] * len(mz_array)
    raw_mz = {
        "params": {"charge": [2]},
        "m/z array": mz_array,
        "intensity array": intensity_array,
    }
    peptide_sequence = "FENLCKIMK"
    error_tol = 0.5
    error_tol_units = "Da"
    result = do_match(
        peptide_sequence=peptide_sequence,
        raw_mz=raw_mz,
        error_tol=error_tol,
        error_tol_units=error_tol_units,
    )
    assert result is not None
    assert result == 1 * len(mz_array)


def test_negative_match():
    mz_array = [148.08,
    978.51,
    277.12,
    849.47,
    391.16,
    735.43,
    504.25,
    622.34,
    607.25,
    519.33,
    735.35,
    391.24,
    848.43,
    278.15,
    979.47,
    147.11][::-1]
    intensity_array = [0] * len(mz_array)
    raw_mz = {
        "params": {"charge": [2]},
        "m/z array": mz_array,
        "intensity array": intensity_array,
    }
    peptide_sequence = "FENLCKIMK"
    error_tol = 0.5
    error_tol_units = "Da"
    result = do_match(
        peptide_sequence=peptide_sequence,
        raw_mz=raw_mz,
        error_tol=error_tol,
        error_tol_units=error_tol_units,
        verbose=False,
    )
    assert result is not None
    assert result == 0