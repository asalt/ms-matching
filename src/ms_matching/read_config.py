def parse_msfragger_config(config_path):
    """
    Parses a simple MSFragger-style key=value config file into a dict.
    Ignores comments, blank lines, and parses numeric values when possible.
    """
    config = {}
    with open(config_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if '=' not in line:
                continue
            key, value = line.split('=', 1)
            key = key.strip()
            value = value.split('#', 1)[0].strip()

            # Try to convert to number
            if value.replace('.', '', 1).replace('-', '', 1).isdigit():
                value = float(value) if '.' in value else int(value)
            elif ',' in value:
                # Try parsing as list of numbers
                items = [v.strip() for v in value.split(',')]
                try:
                    value = [float(v) if '.' in v else int(v) for v in items]
                except ValueError:
                    value = items  # fallback to string list
            config[key] = value
    return config




