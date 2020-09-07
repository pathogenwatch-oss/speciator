import sys
from os.path import basename

import pandas as pd

original_csv = sys.argv[1]
csv_name = basename(original_csv).replace('.tsv', '')

# df = pd.read_csv(original_csv, sep="\t")
# .species.tsv should be read all as strings.
df = pd.read_csv(original_csv, sep="\t", dtype=str)
df.to_parquet(f"{csv_name}.pqt")
