#!/usr/bin/env python3
import sys
import re
import csv

BASES = ("A", "C", "G", "T")

def clean_bases(s: str) -> str:
    s = re.sub(r"\^.", "", s)          # read start + MAPQ char
    s = s.replace("$", "")             # read end

    # remove indels: +3ACT / -2ag etc.
    i = 0
    out = []
    while i < len(s):
        if s[i] in "+-":
            i += 1
            num = []
            while i < len(s) and s[i].isdigit():
                num.append(s[i])
                i += 1
            if num:
                i += int("".join(num))
        else:
            out.append(s[i])
            i += 1
    return "".join(out)

def count_acgt(ref: str, bases: str):
    counts = {b: 0 for b in BASES}
    bases = clean_bases(bases)

    for b in bases:
        if b in ".,":          # reference base on forward/reverse strand
            counts[ref] += 1
        else:
            b = b.upper()
            if b in counts:
                counts[b] += 1

    return counts

writer = csv.writer(sys.stdout)
writer.writerow(["CHR", "POS", "REF", "ALT", "COUNT_REF", "COUNT_ALT"])

for line in sys.stdin:
    if not line.strip():
        continue

    fields = line.rstrip("\n").split("\t")
    if len(fields) < 5:
        continue

    chrom, pos, ref, depth, bases = fields[:5]
    ref = ref.upper()

    if ref not in BASES:
        continue

    counts = count_acgt(ref, bases)
    count_ref = counts[ref]

    for alt in BASES:
        if alt == ref:
            continue

        writer.writerow([
            chrom,
            pos,
            ref,
            alt,
            count_ref,
            counts[alt],
        ])
