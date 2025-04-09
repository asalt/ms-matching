from pyteomics import parser

# example


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

    assert len(peptides) > 0, "No peptides were generated from the protein sequence."


    return peptides



#