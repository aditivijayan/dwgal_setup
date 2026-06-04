from pathlib import Path 

def read_galic_params(param_file):
    values = {}

    with open(param_file) as f:
        for line in f:
            # Strip comments
            line = line.split('%', 1)[0].strip()

            if not line:
                continue

            parts = line.split()
            if len(parts) < 2:
                continue

            key = parts[0]
            val = parts[1]

            values[key] = val

    # Extract what you need
    output_dir = values.get("OutputDir")
    unit_mass  = values.get("UnitMass_in_g")

    if output_dir is None:
        raise KeyError("OutputDir not found in param file")
    if unit_mass is None:
        raise KeyError("UnitMass_in_g not found in param file")

    return Path(output_dir), float(unit_mass)
