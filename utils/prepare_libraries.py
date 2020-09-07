import subprocess
import sys

import pandas


# Requires the .species.pqt & taxon_info.pqt
# Have a directory with all mash sketches in.
# Writes to a given output directory
def write_metadata(library_name: str, full_df, output_path: str):
    # drop(
    #     columns=['superkingdom_code', 'superkingdom_name', 'species_code', 'species_name', 'genus_code',
    #              'genus_name']).
    full_df.rename_axis('taxid').reset_index().to_parquet(f'{output_path}/{library_name}.species.pqt')


def write_library(genus_name: str, strain_df, mash_files_path: str, output_path: str):
    name = f"{genus_name.replace(' ', '_')}.k21s1000"
    write_metadata(name, strain_df, output_path)
    list_filename = f'{output_path}/{name}.txt'
    library_filename = f'{output_path}/{name}'
    strain_names = extract_mash_names(strain_df)
    with open(list_filename, 'w') as list_fh:
        list_fh.write('\n'.join([f'{mash_files_path}/{strainId}.fna.gz.msh' for strainId in strain_names]))
    result = subprocess.run(['mash', 'paste', '-l', library_filename, list_filename], check=True, capture_output=True)
    if result.returncode != 0:
        exit(result.stderr)


def extract_mash_names(strain_table) -> set:
    return {name.replace('.fna.gz', '') for name in strain_table['filename'].array}


data_dir = sys.argv[1]
mash_dir = sys.argv[2]
output_dir = sys.argv[3]
max_sample_size = 1000
do_not_sample = {'NoGenus'}
merge_sets = {
    "Escherichia": {'Salmonella', 'Shigella', 'Citrobacter'},
    'Macrococcus': {'Micrococcus'},
    'Bacteroides': {'Parabacteroides'}
    # 'Klebsiella': {'Raoultella'}
}
skip_genus = set().union(*merge_sets.values())

taxon_data = pandas.read_parquet('data/taxon_info.pqt').rename_axis('taxid').fillna('NoGenus')
all_strains = pandas.read_parquet('data/all_complete_refseq.k21s1000.species.pqt').astype({'taxid': int}).set_index(
    'taxid')

all_strains = all_strains[all_strains['curated_organism_name'] != 'Excluded']

# Split by species
viruses = taxon_data[taxon_data['superkingdom_name'] == 'Viruses']
fungi = taxon_data[taxon_data['superkingdom_name'] == 'Eukaryota']
bacteria = taxon_data[taxon_data['superkingdom_name'] == 'Bacteria']

merge_suffixes = ('_strains', '_taxons')
bacteria_strains = all_strains.merge(bacteria, left_on='taxid', right_on='taxid', how='inner', suffixes=merge_suffixes)
virus_strains = all_strains.merge(viruses, left_on='taxid', right_on='taxid', how='inner', suffixes=merge_suffixes)
fungi_strains = all_strains.merge(fungi, left_on='taxid', right_on='taxid', how='inner', suffixes=merge_suffixes)

# Write the virus & fungi libraries
write_library('Viruses', virus_strains, mash_dir, output_dir)
write_library('Eukaryota', fungi_strains, mash_dir, output_dir)

# Get the bacteria genus codes
genus_names = pandas.unique(bacteria['genus_name'])
genus_reps = dict()
# genus_reps = {genus_code: list() for genus_code in genus_names}

for genus_name in genus_names:
    if genus_name in skip_genus:
        print(f"Skipping {genus_name} for now.")
        continue
    # print(genus_name, file=sys.stderr)
    genus_strains = bacteria_strains[bacteria_strains['genus_name'] == genus_name]
    if genus_name in merge_sets.keys():
        for merge_genus in merge_sets[genus_name]:
            genus_set = bacteria_strains[bacteria_strains['genus_name'] == merge_genus]
            # print(f"Merging {merge_genus} into {genus_name} (extra seqs: {genus_set.shape[0]}")
            genus_strains = genus_strains.append(genus_set)
    # if genus_name == 'Escherichia':
    #     genus_strains = genus_strains.append(bacteria_strains[bacteria_strains['genus_name'] == 'Salmonella']).append(
    #         bacteria_strains[bacteria_strains['genus_name'] == 'Shigella']).append(
    #         bacteria_strains[bacteria_strains['genus_name'] == 'Citrobacter'])
    # elif genus_name == 'Macrococcus':
    #     genus_strains = genus_strains.append(bacteria_strains[bacteria_strains['genus_name'] == 'Micrococcus'])
    # elif genus_name == 'Bacteroides':
    #     genus_strains = genus_strains.append(bacteria_strains[bacteria_strains['genus_name'] == 'Parabacteroides'])
    if genus_strains.shape[0] == 0:
        print(f'No reps for {genus_name}', file=sys.stderr)
        continue
        # Write library for all members
    write_library(genus_name, genus_strains, mash_dir, output_dir)
    species_codes = pandas.unique(genus_strains['species_code'])
    species_reps = list()
    # Get a list of species and pick a representative from each one
    for species_code in species_codes:
        # print(species_code, file=sys.stderr)
        species_reps.append(genus_strains[genus_strains['species_code'] == species_code].sample(1))
    species_selection = pandas.concat(species_reps)
    if genus_name in do_not_sample:
        genus_reps[genus_name] = species_selection
    else:
        sample_size = species_selection.shape[0] if species_selection.shape[0] <= max_sample_size else max_sample_size
        genus_reps[genus_name] = species_selection.sample(sample_size)

# Exclude the NoGenus category from the filter library
del genus_reps['NoGenus']

# Create the merged libraries
merged_strains = pandas.concat(genus_reps.values()).append(virus_strains).append(fungi_strains)
write_library('filter', merged_strains, mash_dir, output_dir)
