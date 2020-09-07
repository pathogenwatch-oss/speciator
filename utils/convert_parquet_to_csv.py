from os.path import basename

import pandas as pd
import sys

original_pqt = sys.argv[1]
pqt_name = basename(original_pqt).replace('.pqt', '')

# df = pd.read_csv(original_csv, sep="\t")
# .species.tsv should be read all as strings.
df = pd.read_parquet(original_pqt)
df.to_csv(f"{pqt_name}.csv")
