

def read_fasta(filepath):
    """
    Generator that reads a FASTA file and yields (header, sequence) tuples line-by-line.

    Args:
        filepath (str): Path to FASTA file

    Yields:
        Tuple[str, str]: (header, sequence)
    """
    header = None
    seq_lines = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, ''.join(seq_lines)
                header = line[1:].strip()  # remove '>'
                seq_lines = []
            else:
                seq_lines.append(line)

        if header:
            yield header, ''.join(seq_lines)
