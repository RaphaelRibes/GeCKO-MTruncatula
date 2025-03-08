#requires python 3.8+
import importlib.metadata
import argparse
import sys

def get_library_version(lib):
    try:
        return importlib.metadata.version(lib)
    except importlib.metadata.PackageNotFoundError:
        return "Not installed"

parser = argparse.ArgumentParser(description="Check Python library versions")
parser.add_argument("-i", "--input", required=True, help="Input file with list of libraries")
parser.add_argument("-o", "--output", required=True, help="Output file for library versions")
args = parser.parse_args()


with open(args.input, "r") as f:
    python_libraries = [line.strip() for line in f if line.strip()]

output_file = args.output


versions = {lib: get_library_version(lib) for lib in python_libraries}

with open(output_file, "w") as f:
    f.write(f"Python\t{sys.version.split()[0]}\n")
    for lib, version in versions.items():
        f.write(f"{lib}\t{version}\n")

