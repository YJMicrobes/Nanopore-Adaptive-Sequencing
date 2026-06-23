#!/usr/bin/env python3
"""
Merge Kraken2 report files into a taxonomy abundance matrix.

Input:
    7_kraken2_contigs/*/kraken2_report.txt

Output:
    kraken_summary_table.csv

Author: YJ
"""


import pandas as pd
import glob
import os


# -------------------------------------------------------
# Input directory
# -------------------------------------------------------

input_dir = "7_kraken2_contigs"

report_files = glob.glob(
    os.path.join(input_dir, "**", "kraken2_report.txt"),
    recursive=True
)


if not report_files:
    raise FileNotFoundError(
        "No Kraken2 report files found. Check input directory."
    )


print(f"Found {len(report_files)} Kraken2 reports")



# -------------------------------------------------------
# Parse Kraken2 reports
# -------------------------------------------------------

dfs = []


for report in sorted(report_files):

    # sample name = parent folder name
    sample = os.path.basename(os.path.dirname(report))

    print(f"Processing: {sample}")


    taxa = []


    with open(report) as fh:

        for line in fh:

            parts = line.rstrip().split("\t")


            if len(parts) < 6:
                continue


            # Kraken2 report format:
            # %  clade_reads  direct_reads  rank  taxid  name

            percent = parts[0].strip()
            name = parts[5].strip()


            taxa.append(
                [
                    name,
                    float(percent)
                ]
            )


    if len(taxa) == 0:
        print(f"Warning: {sample} contains no taxa")
        continue


    df = pd.DataFrame(
        taxa,
        columns=["taxon", sample]
    )


    df = df.set_index("taxon")


    dfs.append(df)



# -------------------------------------------------------
# Merge samples
# -------------------------------------------------------

summary = pd.concat(
    dfs,
    axis=1
).fillna(0)



# -------------------------------------------------------
# Save output
# -------------------------------------------------------

summary.to_csv(
    "kraken_summary_table.csv"
)


print("\nFinished!")
print(
    f"Output: kraken_summary_table.csv "
    f"({summary.shape[0]} taxa x {summary.shape[1]} samples)"
)