import json
import yaml
import sys

def convert_json_to_yaml(json_file, yaml_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    with open(yaml_file, 'w') as f:
        yaml.dump(data, f)

# Usage example
if len(sys.argv) != 3:
    print("Usage: python json_to_yaml.py <input_file> <output_file>")
    sys.exit(1)

json_file = sys.argv[1]
yaml_file = sys.argv[2]
convert_json_to_yaml(json_file, yaml_file)
