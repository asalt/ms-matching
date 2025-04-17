import re


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


def extract_mod_definitions(config: dict, tag_map: dict = None) -> dict:
    tag_map = tag_map or {}

    def parse_entry(item):
        key, value = item
        if not key.startswith("variable_mod_"):
            return None
        parts = re.split(r"\s+", value.strip())
        if len(parts) < 2:
            return None

        delta_mass = float(parts[0])
        targets = parts[1]
        rounded_mass = round(delta_mass, 3)
        tag = tag_map.get(rounded_mass, f"mod{int(abs(delta_mass * 1000))}")
        return tag, {"mass": delta_mass, "targets": targets}

    return {
        tag: mod_info
        for tag, mod_info in map(parse_entry, config.items())
        if tag is not None
    }
