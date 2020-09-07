import pandas

# types = {'taxonId': object, 'superkingdomId': str, 'genusId': str, 'speciesId': str}
codes = pandas.read_csv('raw_data/full_taxon_codes.csv', header=0, dtype=str, index_col=0)
names = pandas.read_csv('raw_data/full_taxon_names.csv', header=0, dtype=str, index_col=0)

combined = codes.join(names, on='taxonId', lsuffix='_code', rsuffix='_name')
combined.to_parquet('data/taxon_info.pqt', index=True)
