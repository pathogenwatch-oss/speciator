import csv
import re
import sys

import pandas as pd


# header = ',taxid,accession,GCF_accession,GCF_accession_without_version,filename,length,num_contigs,info,refseq_organism_name,refseq_full_organism_name,bacsort_organism_name,curated_organism_name,species_taxid,bioproject,biosample,refseq_category,infraspecific_name,assembly_level,asm_name,submitter,ftp_path,superkingdom_code,genus_code,species_code,superkingdom_name,genus_name,species_name'
# header_fields = header.split(',')


def extract_organism_name(raw: str) -> str:
    if raw == 'Salmonella_enterica/0001.fna.gz':
        return 'Salmonella enterica subsp. enterica serovar Typhi'
    subspecies_matcher = re.compile(r'\ssubsp\s\w+$')
    raw = raw.split('/')[0].replace('_', ' ').replace('.fna', '')
    return subspecies_matcher.sub('', raw)


def update_shigella_record(df, csv_row):
    accession = csv_row['kleborate']
    species_name = csv_row['species']
    species_code = csv_row['taxid']
    df.loc[df['accession'] == accession, 'genus_name'] = 'Shigella'
    df.loc[df['accession'] == accession, 'genus_code'] = '620'
    df.loc[df['accession'] == accession, 'species_name'] = species_name
    df.loc[df['accession'] == accession, 'refseq_organism_name'] = species_name
    df.loc[df['accession'] == accession, 'refseq_full_organism_name'] = species_name
    df.loc[df['accession'] == accession, 'bacsort_organism_name'] = species_name
    df.loc[df['accession'] == accession, 'curated_organism_name'] = species_name
    df.loc[df['accession'] == accession, 'species_code'] = species_code
    df.loc[df['accession'] == accession, 'taxid'] = species_code
    df.loc[df['accession'] == accession, 'species_taxid'] = species_code
    return df


# Script start
typhi_organism_name = 'Salmonella enterica subsp. enterica serovar Typhi'
mash_info = {'curated_organism_name': [], 'accession': [], 'length': []}

with open('kleborate_library/mash_info.txt', 'r') as fh:
    for line_number in range(0, 10):
        fh.readline()
    for line in fh.readlines():
        line = line.strip()
        if line == '':
            continue
        data = line.split()
        # Skip comment extension lines for now.
        if data[0] != '1000':
            continue
        mash_info['curated_organism_name'].append(extract_organism_name(data[2]))
        mash_info['accession'].append(data[2])
        mash_info['length'].append(int(data[1]))

mash_df = pd.DataFrame.from_dict(mash_info)
print(f'{mash_df.shape[0]} rows in original mash_info.', file=sys.stderr)

# Add taxon IDs
taxon_info_df = pd.read_parquet('data/taxon_info.pqt')
merged = mash_df.merge(taxon_info_df, left_on='curated_organism_name', right_on='species_name', how='inner')
merged = merged.drop_duplicates()
print(f'{merged.shape[0]} rows in merged df.', file=sys.stderr)

# Add missing columns
merged['taxid'] = merged['species_code']
merged['species_taxid'] = merged['species_code']
merged['GCF_accession'] = merged['accession']
merged['GCF_accession_without_version'] = merged['accession']
merged['filename'] = merged['accession']
merged['num_contigs'] = 1
merged['info'] = ''
merged['refseq_organism_name'] = merged['curated_organism_name']
merged['refseq_full_organism_name'] = merged['curated_organism_name']
merged['bacsort_organism_name'] = merged['curated_organism_name']
merged = merged.assign(bioproject='', biosample='', refseq_category='', infraspecific_name='', assembly_level='',
                       asm_name='', submitter='', ftp_path='')

# Add Typhi representative
typhi_dict = {
    'curated_organism_name': [typhi_organism_name],
    'accession': ['Salmonella_enterica/0001.fna.gz'],
    'GCF_accession': ['Salmonella_enterica/0001.fna.gz'],
    'GCF_accession_without_version': ['Salmonella_enterica/0001.fna.gz'],
    'filename': ['Salmonella_enterica/0001.fna.gz'],
    'num_contigs': [1],
    'info': '',
    'refseq_organism_name': typhi_organism_name,
    'refseq_full_organism_name': typhi_organism_name,
    'bacsort_organism_name': typhi_organism_name,
    'length': [5004298],
    'taxid': ['90370'],
    'species_taxid': ['28901'],
    'species_code': ['28901'],
    'species_name': ['Salmonella enterica'],
    'genus_code': ['590'],
    'genus_name': ['Salmonella'],
    'superkingdom_code': ['2'],
    'superkingdom_name': ['Bacteria'],
}

typhi_df = pd.DataFrame.from_dict(typhi_dict)

sys.stderr.write(f'Before: {merged.shape[0]}')
merged = merged.append(typhi_df, ignore_index=True)
sys.stderr.write(f' after: {merged.shape[0]}\n')

# Add Shigella fixes
shigella_dict = list()
with open('kleborate_library/shigella_links.csv') as sl_fh:
    shigella_updates = csv.DictReader(sl_fh)
    for shigella_update in shigella_updates:
        merged = update_shigella_record(merged, shigella_update)

# merged.loc[merged['curated_organism_name'] == typhi_organism_name, 'taxid'] = '90370'
# merged.loc[merged['curated_organism_name'] == typhi_organism_name, 'species_taxid'] = '28901'
# merged.loc[merged['curated_organism_name'] == typhi_organism_name, 'species_code'] = '28901'
# merged.loc[merged['curated_organism_name'] == typhi_organism_name, 'species_name'] = 'Salmonella enterica'
# merged.loc[merged['curated_organism_name'] == typhi_organism_name, 'genus_code'] = '590'
# merged.loc[merged['curated_organism_name'] == typhi_organism_name, 'genus_name'] = 'Salmonella'
# merged.loc[merged['curated_organism_name'] == typhi_organism_name, 'superkingdom_code'] = '2'
# merged.loc[merged['curated_organism_name'] == typhi_organism_name, 'superkingdom_name'] = 'Bacteria'

# merged = merged.drop(
#     columns=['species_code', 'genus_name', 'genus_code', 'superkingdom_code', 'superkingdom_name']).rename(
#     columns={'species_name': 'organism_name'})

# merged = merged[
#     ['accession', 'GCF_accession', 'GCF_accession_without_version', 'filename', 'length', 'num_contigs', 'info',
#      'refseq_organism_name', 'refseq_full_organism_name', 'bacsort_organism_name', 'curated_organism_name', 'taxid',
#      'species_taxid', 'bioproject', 'biosample', 'refseq_category', 'infraspecific_name', 'assembly_level',
#      'asm_name', 'submitter', 'ftp_path']]

# merged['bioproject'] = ''
# merged['biosample'] = ''
# merged['refseq_category'] = ''

merged.to_parquet('kleborate_library/Kleborate.k21s1000.species.pqt')
merged.to_csv('kleborate_library/Kleborate.k21s1000.species.csv')
