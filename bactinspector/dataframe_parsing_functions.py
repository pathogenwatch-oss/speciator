import os

import numpy as np
import pandas as pd


def create_refseq_species_info_df(local_mash_and_info_file_prefix, ref_and_rep_only=False):
    #  read in species match table
    if local_mash_and_info_file_prefix is not None:
        refseq_species_info_file = local_mash_and_info_file_prefix + '.species.pqt'
    else:
        refseq_species_info_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data',
                                                'all_complete_bacteria_refseq.k21s1000.species.pqt')
    refseq_species_info = pd.read_parquet(refseq_species_info_file)

    # filter only reference and representative genomes if specified as an option
    if ref_and_rep_only:
        refseq_species_info = refseq_species_info.loc[
            (refseq_species_info['refseq_category'] == 'representative genome') |
            (refseq_species_info['refseq_category'] == 'reference genome')
            ]
    return refseq_species_info


def create_refseq_species_metrics_df():
    """
    Create a pandas dataframe based on the per species distance metrics found in the file
    all_complete_bacteria_refseq.k21s1000.species.distance_and_length_metrics.pqt
    """
    # merge with species distances cutoff and length
    refseq_species_metrics_file = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'data',
        'all_complete_bacteria_refseq.k21s1000.species.distance_and_length_metrics.pqt'
    )
    refseq_species_metrics_df = pd.read_parquet(refseq_species_metrics_file)

    # add adjusted distance 0.1 for max_distance > 0.1
    refseq_species_metrics_df['adjusted_max_distance'] = np.where(refseq_species_metrics_df['max_distance'] > 0.1, 0.1,
                                                                  refseq_species_metrics_df['max_distance'])

    # where there are few distances: replace max_length and min_length with 'unknown'
    refseq_species_metrics_df['max_length'] = np.where(refseq_species_metrics_df['max_distance'] == 0, np.nan,
                                                       refseq_species_metrics_df['max_length'])
    refseq_species_metrics_df['min_length'] = np.where(refseq_species_metrics_df['max_distance'] == 0, np.nan,
                                                       refseq_species_metrics_df['min_length'])

    # fill NaNs with 'Unknown'
    refseq_species_metrics_df = refseq_species_metrics_df.fillna('unknown')
    return refseq_species_metrics_df


def add_certainty_to_merged_results_df(results_df,
                                       num_best_matches=10,
                                       allowed_variance=0.1,
                                       allowed_variance_rarer_species=0.5
                                       ):
    # ====== use following tests ======
    # if max_distance is 0 use 0.05 as cutoff
    # if no max_distance  use 0.05 as cutoff
    # compare with max distance + allowed_variance_rarer_species if species has < 10 genomes
    # compare with max distance + allowed_variance_rarer_species if species has >= 10 genomes
    # if no matches
    # if percentage of top 10 hits (normalised for number genomes)
    results_df['result'] = np.where(
        ((results_df['max_distance'] == 0) & (results_df['top_hit_distance'] > 0.05)) |
        (results_df['max_distance'].isna() & results_df['top_hit_distance'] > 0.05) |
        ((results_df['num_genomes'] < 10) & (results_df['top_hit_distance'] > results_df['adjusted_max_distance'] * (
                1 + allowed_variance_rarer_species))) |
        ((results_df['num_genomes'] >= 10) & (
                results_df['top_hit_distance'] > results_df['adjusted_max_distance'] * (1 + allowed_variance))) |
        (results_df['%_of_{0}_best_matches=species'.format(num_best_matches)].isna()) |
        (results_df['%_of_{0}_best_matches=species'.format(num_best_matches)] * num_best_matches / np.minimum(
            num_best_matches, results_df['num_genomes']) < 60),
        'uncertain', 'good'
    )
    # ==================================

    # convert lengths to int
    results_df['max_length'] = pd.to_numeric(results_df['max_length'], errors='ignore', downcast='integer')
    results_df['min_length'] = pd.to_numeric(results_df['min_length'], errors='ignore', downcast='integer')

    # rename lengths
    results_df = results_df.rename(
        columns={
            'max_length': 'maximum_genome_length',
            'min_length': 'minimum_genome_length'
        }
    )
    # drop columns no longer required
    results_df = results_df.drop(columns=[
        'max_distance',
        'mean_distance',
        'std_distance',
        'quartile_75_distance',
        'quartile_25_distance',
        'adjusted_max_distance']
    )
    return results_df
