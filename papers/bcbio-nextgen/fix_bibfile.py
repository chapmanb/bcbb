import sys

in_file, out_file = sys.argv[1:]

with open(in_file) as in_handle:
    with open(out_file, "w") as out_handle:
        for line in in_handle:
            if line.strip().startswith("language"):
                line = None
            elif line.strip().startswith("@misc"):
                base, cite = line.split("{")
                cite_parts = [x.strip() for x in cite.split("_") if x.strip()]
                line = base + "{" + cite_parts[0] + ",\n"
            if line:
                out_handle.write(line)
