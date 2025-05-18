from pyteomics import parser, mass


# we can also use pyteomics.mass.fast_mass to perform all of these calculations for us
def make_fragments(
    peptide, types=("b", "y"), maxcharge=1, aa_mass=mass.std_aa_mass
):  # from https://pyteomics.readthedocs.io/en/latest/examples/example_msms.html
    """
    The function generates all possible m/z for fragments of types
    `types` and of charges from 1 to `maxharge`.
    """
    # import pdb; pdb.set_trace()
    for i in range(1, len(peptide)):
        for ion_type in types:
            for charge in range(1, maxcharge + 1):
                if ion_type[0] in "abc":
                    yield mass.fast_mass2(
                        peptide[:i], ion_type=ion_type, charge=charge, aa_mass=aa_mass
                    )
                else:
                    yield mass.fast_mass2(
                        peptide[i:], ion_type=ion_type, charge=charge, aa_mass=aa_mass
                    )


def do_match(
    peptide_sequence,
    raw_mz,
    error_tol=20,
    error_tol_units="ppm",
    verbose=True,
    maxcharge=None,
    **kwargs,
):
    """
    Match the peptide sequence to the raw m/z data.
    Parameters
    peptide_sequence: str
        The peptide sequence to match.
    raw_mz: dict
        The raw m/z data to match against.
        from pyteomics.ms2.read or similar, with keys:
            - "params": parameters for the raw data
            - "m/z array": array of m/z values
            - "intensity array": array of intensity values
    error_tol: float
        The error tolerance for the match.
    error_tol_units: str
        The units for the error tolerance. Can be 'ppm' or 'Da'.
    """
    if error_tol_units not in ("ppm", "Da"):
        raise ValueError(f"error_tol_units must be one of ('ppm', 'Da')")
    if error_tol_units == "ppm":
        error_tol = (
            error_tol * 1e-6 * raw_mz["m/z array"].mean()
        )  # convert ppm to Da, approximating the mean m/z as the average observed value

    params = raw_mz["params"]
    parent_ion_charge = int(params["charge"][0])
    if maxcharge is None:
        maxcharge = parent_ion_charge - 1
    mz_array = raw_mz["m/z array"]
    intensity_array = raw_mz["intensity array"]

    b_fragments = list(
        make_fragments(peptide_sequence, types="b", maxcharge=maxcharge)
    )
    y_fragments = list(
        make_fragments(peptide_sequence, types="y", maxcharge=maxcharge)
    )

    # how does this work?
    # could be rewritten to be more clear
    b_matches = [
        x for x in b_fragments if any(abs(x - y) < error_tol for y in mz_array)
    ]
    y_matches = [
        x for x in y_fragments if any(abs(x - y) < error_tol for y in mz_array)
    ]
    # how do we get the intensities for the matched fragments?

    # solution
    assert len(mz_array) == len(
        intensity_array
    ), "m/z array and intensity array must be the same length"
    # we can assume this will always be true, but we can check it anyway

    b_fragment_intensities = []
    for b_fragment in b_fragments:
        for mz, intensity in zip(mz_array, intensity_array):
            if abs(b_fragment - mz) < error_tol:
                if verbose:
                    print(f"Match found for b fragment {b_fragment} at m/z {mz}")
                b_fragment_intensities.append(intensity)
                break  # break because we have found the match

    y_fragment_intensities = []
    for y_fragment in y_fragments:
        for mz, intensity in zip(mz_array, intensity_array):
            if abs(y_fragment - mz) < error_tol:
                if verbose:
                    print(
                        f"Match found for y fragment {y_fragment} at m/z {mz} with intensity {intensity}"
                    )
                b_fragment_intensities.append(intensity)
                break  # break because we have found the match

    total_intensity = sum(b_fragment_intensities) + sum(y_fragment_intensities)
    n_matches = len(b_fragment_intensities) + len(y_fragment_intensities)
    if verbose:
        print(f"Number of matches: {n_matches}")
        print(f"Total intensity: {total_intensity}")

    # ideally this would return something more informative; a proper match score
    return total_intensity
