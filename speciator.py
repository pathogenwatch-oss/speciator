import json
import os
import sys

import pandas

from bactinspector import commands


class AttributeDict(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__


def default_result() -> dict:
    return {
        'taxId': '32644',
        'speciesId': '32644',
        'speciesName': 'Unidentified',
        'genusId': '9999999',
        "genusName": 'Unclassified',
        "superkingdomId": '12908',
        "superkingdomName": 'Unclassified Sequences',
        "referenceId": 'None found',
        "mashDistance": 1.0,
        "pValue": 0,
        "matchingHashes": '0/1000',
        "confidence": 'None',
        "source": 'None'
    }


def build_pw_result(result_df, species_datafile, lib_name):
    if result_df.shape[0] == 0:
        return default_result()

    result = result_df.iloc[0]
    species_info = pandas.read_parquet(species_datafile)

    if result['species'] == 'No significant matches':
        return default_result()
    else:
        if species_info[species_info['species_name'] == result['species']].shape[0] > 0:
            species_md = species_info[species_info['species_name'] == result['species']].iloc[0]
        else:
            species_md = species_info[species_info['species_code'] == result['species_taxid']].iloc[0]
        if species_md['species_code'] != result['species_taxid']:
            # This is a re-written reference
            strain_id = species_md['species_code']
        else:
            strain_id = result['strain_taxid']
        species_md.fillna(value='', inplace=True)
        # print(f"Strain ID: {strain_id}; superkingdom_code: {species_md['superkingdom_code']}", file=sys.stderr)
        return {
            'taxId': strain_id,
            'speciesId': species_md['species_code'],
            'speciesName': species_md['species_name'],
            'genusId': species_md['genus_code'] if species_md['genus_code'] != '' else '9999999',
            "genusName": species_md['genus_name'] if species_md['genus_name'] != '' else 'Unclassified',
            "superkingdomId": species_md['superkingdom_code'],
            "superkingdomName": species_md['superkingdom_name'],
            "referenceId": result['top_hit'],
            "mashDistance": float(result['top_hit_distance']),
            "pValue": float(result['top_hit_p_value']),
            "matchingHashes": result['top_hit_shared_hashes'],
            "confidence": result['result'],
            "source": lib_name
        }


def run_mass_search(num_best_matches=20, taxon_info='bactinspector/data/taxon_info.pqt'):
    # Initial filter
    # print(f'Running filter search', file=sys.stderr)
    filter_result = build_pw_result(commands.run_check_species(
        AttributeDict({'distance_cutoff': 0.15,
                       'fasta_file_pattern': fasta_file,
                       'input_dir': input_dir,
                       'output_dir': '/tmp/',
                       'local_mash_and_info_file_prefix': f'{libraries_path}/filter.k21s1000',
                       'stdout_summary': True,
                       'num_best_matches': num_best_matches,
                       'parallel_processes': 1,
                       'mash_path': '',
                       'allowed_variance': 0.1,
                       'allowed_variance_rarer_species': 0.5
                       })), taxon_info, 'filter')

    genus = filter_result['genusName'].replace(' ', '_') if filter_result['genusName'] != 'Unclassified' else 'NoGenus'

    # print(f"Genus: {genus}", file=sys.stderr)
    collected_groups = {"Escherichia": {'Escherichia', 'Salmonella', 'Shigella', 'Citrobacter'},
                        'Macrococcus': {'Macrococcus', 'Micrococcus'},
                        'Bacteroides': {'Bacteroides', 'Parabacteroides'},
                        'Klebsiella': {'Klebsiella', 'Raoultella'}}

    if filter_result['superkingdomName'] != 'Bacteria' and filter_result[
        'superkingdomName'] != 'Unclassified Sequences':
        library, threshold = filter_result['superkingdomName'], 0.075
    elif genus in collected_groups['Escherichia']:
        library, threshold = 'Escherichia', 0.05
    elif genus in collected_groups['Klebsiella']:
        library, threshold = 'Klebsiella', 0.05
    elif genus in collected_groups['Macrococcus']:
        library, threshold = 'Macrococcus', 0.05
    elif genus in collected_groups['Bacteroides']:
        library, threshold = 'Bacteroides', 0.05
    elif genus == 'NoGenus':
        library, threshold = genus, 0.075
    else:
        library, threshold = genus, 0.05

    # print(f'Library={library}, threshold={threshold}, num_best_matches={num_best_matches}', file=sys.stderr)

    library_file = f'{libraries_path}/{library}.k21s1000'

    results_df = commands.run_check_species(
        AttributeDict({'distance_cutoff': threshold,
                       'fasta_file_pattern': fasta_file,
                       'input_dir': input_dir,
                       'output_dir': '/tmp/',
                       'local_mash_and_info_file_prefix': library_file,
                       'stdout_summary': True,
                       'num_best_matches': num_best_matches,
                       'parallel_processes': 1,
                       'mash_path': '',
                       'allowed_variance': 0.1,
                       'allowed_variance_rarer_species': 0.5
                       }))

    # Check if we need to try the hit n hope.
    if results_df.iloc[0]['species'] == 'No significant matches':
        library = 'Fallback'
        results_df = commands.run_check_species(
            AttributeDict({'distance_cutoff': threshold,
                           'fasta_file_pattern': fasta_file,
                           'input_dir': input_dir,
                           'output_dir': '/tmp/',
                           'local_mash_and_info_file_prefix': f'{libraries_path}/NoGenus.k21s1000',
                           'stdout_summary': True,
                           'num_best_matches': num_best_matches,
                           'parallel_processes': 1,
                           'mash_path': '',
                           'allowed_variance': 0.1,
                           'allowed_variance_rarer_species': 0.5
                           })
        )
    return results_df, library


fasta_path = sys.argv[1]
libraries_path = sys.argv[2]
taxon_info = sys.argv[3]

# print(f'Running {fasta_path} against {libraries_path}\n', file=sys.stderr)
input_dir = os.path.dirname(fasta_path)
fasta_file = os.path.basename(fasta_path)

# Curated library search
# print(f'Running curated library search', file=sys.stderr)
species_assignment = commands.run_check_species(
    AttributeDict({'distance_cutoff': 0.04,
                   'fasta_file_pattern': fasta_file,
                   'input_dir': input_dir,
                   'output_dir': '/tmp/',
                   'local_mash_and_info_file_prefix': f'{libraries_path}/Kleborate.k21s1000',
                   'stdout_summary': True,
                   'num_best_matches': 1,
                   'parallel_processes': 1,
                   'mash_path': '',
                   'allowed_variance': 0.1,
                   'allowed_variance_rarer_species': 0.5
                   }))
tag = 'Curated'
if species_assignment.iloc[0]['species'] == 'No significant matches' or species_assignment.iloc[0][
      'species'] in ['Escherichia coli', 'Shigella sonnei', 'Shigella flexneri']:
    species_assignment, tag = run_mass_search(taxon_info=taxon_info)

json_result = build_pw_result(species_assignment, taxon_info, tag)
print(json.dumps(json_result), file=sys.stdout)
