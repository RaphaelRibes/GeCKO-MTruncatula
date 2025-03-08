import yaml
import sys

def load_yaml(file_path):
    try:
        with open(file_path, "r") as file:
            return yaml.safe_load(file)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.", file=sys.stderr)
        sys.exit(1)
    except yaml.YAMLError as e:
        print(f"Error: Unable to read the YAML file '{file_path}'.\n{e}", file=sys.stderr)
        sys.exit(1)

if len(sys.argv) != 3:
    print("Usage: python3 parse_yaml.py <path_to_yaml_file> <path_to_model_yaml>", file=sys.stderr)
    sys.exit(1)

yaml_file = sys.argv[1]
model_file = sys.argv[2]

config = load_yaml(yaml_file)
model_config = load_yaml(model_file)

if not isinstance(config, dict) or not isinstance(model_config, dict):
    print("Both YAML files must be lists of simple key-value pairs.", file=sys.stderr)
    sys.exit(1)

for expected_var in model_config.keys():
    if expected_var not in config:
        print(f"\nERROR: The expected variable '{expected_var}' was not found in your config file ({yaml_file}). Please make sure to include it.", file=sys.stderr)
        print(f"\nList of expected variables:\n{', '.join(model_config.keys())}", file=sys.stderr)
        print("\nAs a reminder:", file=sys.stderr)
        print("All expected variables must be specified in YAML format, with one variable per row followed by a ':' and a space before its assigned value. For example: VARIABLE_NAME: value.", file=sys.stderr)
        print("To leave a variable empty, use quotes: VARIABLE_NAME: \"\" ", file=sys.stderr)
        #print("\nExiting.", file=sys.stderr)
        sys.exit(1)


for key, value in config.items():
    print(f"{key}\t{value}")
