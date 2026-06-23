#!/usr/bin/env python3
import re

xmfa = "4genomes.xmfa"
outfile = "links.txt"

links = []
block = []

with open(xmfa) as f:
    for line in f:

        # alignment block header
        if line.startswith(">"):
            m = re.search(r">\s*(\d+):(\d+)-(\d+)", line)
            if m:
                genome = "g" + m.group(1)
                start = int(m.group(2))
                end = int(m.group(3))
                block.append((genome, start, end))

        # end of alignment block
        elif line.startswith("="):
            if len(block) > 1:
                ref = block[0]
                for other in block[1:]:
                    links.append((ref, other))
            block = []

# write Circos links
with open(outfile, "w") as out:
    for ref, other in links:
        out.write(
            f"{ref[0]} {ref[1]} {ref[2]} "
            f"{other[0]} {other[1]} {other[2]}\n"
        )

print("Links written:", len(links))