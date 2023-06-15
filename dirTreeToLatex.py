import os
import sys

def print_dirtree(dir_path, level=1):
    for entry in sorted(os.listdir(dir_path)):
        entry_path = os.path.join(dir_path, entry)

        # Escape underscores in the filename
        escaped_entry = entry.replace('_', '\\_')

        # Print the current entry
        print(f".{level} {escaped_entry}.")

        # If the entry is a directory, recursively print its content
        if os.path.isdir(entry_path):
            print_dirtree(entry_path, level+1)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python dirtree.py <directory_path>")
        sys.exit(1)

    dir_path = sys.argv[1]

    if not os.path.isdir(dir_path):
        print(f"Error: {dir_path} is not a directory.")
        sys.exit(1)

    print("\\dirtree{%")
    print_dirtree(dir_path)
    print("}")
