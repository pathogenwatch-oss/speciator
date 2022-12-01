import os
import sys
from io import StringIO

import pandas as pd

from bactinspector.utility_functions import add_new_file_extension, get_base_name, run_command


def run_mash_sketch(file, filetype, output_dir=None, mash_path=''):
    """
    run mash sketch on a fasta file and return the path to the resulting sketch file
    """
    if output_dir:
        sketch_file = os.path.join(output_dir, '{0}.msh'.format(get_base_name(file)))
    else:
        sketch_file = add_new_file_extension(file, 'msh')

    if not os.path.exists(sketch_file) or os.path.getsize(sketch_file) == 0:
        # sys.stderr.write('Sketching {0}\n'.format(get_base_name(file)))
        if filetype == 'fasta':
            command_and_arguments = [os.path.join(mash_path, 'mash'), 'sketch', file, '-o', sketch_file]
        else:
            command_and_arguments = [os.path.join(mash_path, 'mash'), 'sketch', '-m', '3', file, '-o', sketch_file]
        ret_code, out, err = run_command(command_and_arguments)
        if ret_code != 0:
            sys.stderr.write('Error whilst performing mash sketch: {0}\n'.format(err))
            sys.exit(ret_code)
    return sketch_file


def get_best_mash_matches(sample_sketch, ref_seq_sketch, mash_path='', number_of_best_matches=10,
                          distance_threshold=0.5):
    """
    run mash dist sample sketch file vs the ref_seq sketches and return the best matches
    """
    match_file = add_new_file_extension(sample_sketch, 'best_matches.txt')
    # TODO: Add unhappy paths
    if not os.path.exists(match_file) or os.path.getsize(match_file) == 0:
        mash_dists = execute_mashing(mash_path, ref_seq_sketch, sample_sketch, distance_threshold)
        # merge with refseq matches for potential filtering
        # mash_dists = mash_dists.merge(refseq_species_info, left_on='subject', right_on='filename', how='right')
        mash_dists = mash_dists.filter(['query', 'subject', 'distance', 'p-value', 'shared-hashes'])
        # sort by distance and output the subjects (match in refseq) 
        matches = mash_dists.sort_values('distance', ascending=True).head(number_of_best_matches)
        matches = matches.rename(columns={'subject': 'filename'}).filter(
            items=['filename', 'distance', 'p-value', 'shared-hashes'])

    return get_base_name(sample_sketch), matches


def execute_mashing(mash_path, ref_seq_sketch, sample_sketch, distance_threshold):
    # sys.stderr.write(f'Getting best match for {get_base_name(sample_sketch)}\n')
    # time1 = time.process_time()
    command_and_arguments = [os.path.join(mash_path, 'mash'), 'dist', '-d', str(distance_threshold), sample_sketch,
                             ref_seq_sketch]

    # print(command_and_arguments, file=sys.stderr)
    # result = subprocess.run(command_and_arguments, stdout=subprocess.PIPE, stderr=True, text=True)
    ret_code, out, err = run_command(command_and_arguments, text=True)
    if ret_code != 0:
        print('Error whilst performing mash dist: {0}'.format(err))
        sys.exit(ret_code)

    distances_fh = StringIO(out)
    mash_dists = pd.read_csv(distances_fh, sep="\t", names=['query', 'subject', 'distance', 'p-value', 'shared-hashes'])
    # if result.returncode != 0:
    #     print('Error whilst performing mash dist: {0}'.format(result.stderr))
    #     sys.exit(result.returncode)
    # time2 = time.process_time()
    # print(f'Mash time {time2 - time1}', file=sys.stderr)
    return mash_dists


def get_species_match_details(matches, refseq_species_info):
    """
    use pandas to merge best matches  with ref species info and return the merged data frame
    """

    best_match_species_df = matches.merge(
        refseq_species_info,
        on=['filename']
    )
    return best_match_species_df


def get_most_frequent_species_match(matches, refseq_species_info, distance_cutoff=0.05):
    """
    use pandas to merge best match file with ref species info and report the most frequent species
    return species and count
    """
    best_match_species_df = get_species_match_details(matches, refseq_species_info)
    # filter for close matches
    # best_match_species_df = best_match_species_df.loc[best_match_species_df['distance'] <= distance_cutoff]
    if len(best_match_species_df) == 0:
        return 'No significant matches', None, None, None, None, None, None, None, None
    else:
        # get most frequent species and count
        species_name_counts = best_match_species_df['curated_organism_name'].value_counts()
        if 1 == species_name_counts.size or species_name_counts[0] != species_name_counts[1]:
            most_frequent_species_name = species_name_counts.index[0]
            most_frequent_species_count = species_name_counts[0]
        else:
            equal_often_names = species_name_counts[species_name_counts == species_name_counts[0]].index.tolist()
            best_match = pd.concat(
                [best_match_species_df.loc[best_match_species_df['curated_organism_name'] == test_name].sort_values(
                    'distance').head(1) for test_name in equal_often_names]).sort_values('distance').iloc[0, :]
            most_frequent_species_name = best_match['curated_organism_name']
            most_frequent_species_count = species_name_counts[0]

        # get top hit of the most frequent species as measured by distance
        taxid_counts = best_match_species_df['species_code'].value_counts()
        # print(best_match_species_df.to_csv(sys.stderr))
        if 1 == taxid_counts.size or taxid_counts[0] != taxid_counts[1]:
            most_frequent_species_taxid = taxid_counts.index[0]
            top_hit = best_match_species_df.loc[
                          best_match_species_df['species_code'] == most_frequent_species_taxid].sort_values(
                'distance').iloc[0, :]
        else:
            equal_often_taxids = taxid_counts[taxid_counts == taxid_counts[0]].index.tolist()
            top_hit = pd.concat(
                [best_match_species_df.loc[best_match_species_df['species_code'] == test_taxid].sort_values(
                    'distance').head(1) for test_taxid in equal_often_taxids]).sort_values('distance').iloc[0, :]

        return (
            most_frequent_species_name,
            most_frequent_species_count,
            str(top_hit['species_code']),
            str(top_hit.taxid),
            top_hit.accession,
            len(best_match_species_df),
            top_hit['distance'],
            top_hit['p-value'],
            top_hit['shared-hashes']
        )
