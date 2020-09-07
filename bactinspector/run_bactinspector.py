#!/usr/bin/env python3
import argparse
import os
import sys
import textwrap

import pkg_resources

from bactinspector.commands import run_check_species, run_closest_match, run_create_species_info, run_info


def is_valid_file(parser, arg):
    if not os.path.isfile(arg):
        parser.error('The file {} does not exist!'.format(arg))
    else:
        # File exists so return the filename
        return arg


def is_valid_dir(parser, arg):
    if not os.path.isdir(arg):
        parser.error('The directory {} does not exist!'.format(arg))
    else:
        # File exists so return the directory
        return arg


class Version(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        print(pkg_resources.require("BactInspectorMax")[0].version)
        sys.exit(0)


def parse_arguments():
    description = textwrap.dedent("""
    A module to determine the most probable species based on sequence in fasta files using refseq and Mash (https://mash.readthedocs.io/en/latest/index.html)
    It will count the species of the top ref seq mash matches and report most frequent.

    In order to use the module:
      • Specify an input directory and output directory (default is current directory)
      • Specify either a 
        • fasta file pattern with -f (e.g "*.fas") or 
        • mash sketch file pattern with -m (e.g "*.msh") if you have already sketched the fasta files
      • By default the top 10 matches will be used. Change this with -n
      • Speed things up by changing the number of parallel processes to match the cores on your computer using -p
      • If mash is not in your PATH specify the directory containing the mash executable with -mp

    If you want to update the genomes used, follow the instructions on https://gitlab.com/antunderwood/bactinspector/wikis/Updating-the-genomes-in-BactInspector
    and use the create_species_info command to make the required file
    """)
    # parse all arguments
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.register('action', 'version', Version)
    parser.add_argument('-v', '--version', action='version', nargs=0, help='print out software version')

    subparsers = parser.add_subparsers(
        help='The following commands are available. Type bactinspector <COMMAND> -h for more help on a specific commands',
        dest='command'
    )

    subparsers.required = True

    # check species sub command
    check_species_command = subparsers.add_parser('check_species',
                                                  help='Check the most frequent matches to a species in refseq'
                                                  )

    check_species_command.add_argument('-i', '--input_dir', help='path to input_directory',
                                       type=lambda x: is_valid_dir(parser, x), default='.')
    check_species_command.add_argument('-o', '--output_dir', help='path to output_directory',
                                       type=lambda x: is_valid_dir(parser, x), default='.')
    check_species_command.add_argument('-p', '--parallel_processes', help='number of processes to run in parallel',
                                       default=1, type=int)
    check_species_command.add_argument('-n', '--num_best_matches', help='number of best matches to return', default=10,
                                       type=int)
    check_species_command.add_argument('-d', '--distance_cutoff', help='mash distance cutoff (default 0.05)',
                                       default=0.05, type=float)
    check_species_command.add_argument('-v', '--allowed_variance',
                                       help='proportion of max_distance allowed over which a result will be marked as uncertain (default 0.1)',
                                       default=0.1, type=float)
    check_species_command.add_argument('-vl', '--allowed_variance_rarer_species',
                                       help='proportion of max_distance allowed over which a result will be marked as uncertain for species which have fewer than 10 representatives in refseq (default 0.5)',
                                       default=0.5, type=float)
    check_species_command.add_argument('-s', '--stdout_summary', help='output a summary of the result to STDOUT',
                                       action='store_true')
    check_species_command.add_argument('-l', '--local_mash_and_info_file_prefix',
                                       help='the path prefix to the mash sketch file and corresponding info file')

    check_species_command.add_argument('-mp', '--mash_path',
                                       help='path to the mash executable. If not provided it is assumed mash is in the PATH')

    filetype_extension = check_species_command.add_mutually_exclusive_group(required=True)
    filetype_extension.add_argument('-f', '--fasta_file_pattern', help='pattern to match fasta files e.g "*.fas"')
    filetype_extension.add_argument('-fq', '--fastq_file_pattern', help='pattern to match fastq files e.g "*.fastq.gz"')
    filetype_extension.add_argument('-m', '--mash_sketch_file_pattern',
                                    help='pattern to match mash sketch files e.g "*.msh"')

    # closest match sub command
    closest_match_command = subparsers.add_parser('closest_match',
                                                  help='Report the closest matches to a set of sequences'
                                                  )

    closest_match_command.add_argument('-i', '--input_dir', help='path to input_directory',
                                       type=lambda x: is_valid_dir(parser, x), default='.')
    closest_match_command.add_argument('-o', '--output_dir', help='path to output_directory',
                                       type=lambda x: is_valid_dir(parser, x), default='.')
    closest_match_command.add_argument('-p', '--parallel_processes', help='number of processes to run in parallel',
                                       default=1, type=int)
    closest_match_command.add_argument('-r', '--ref_and_rep_only',
                                       help='only include reference and representative sequences', action='store_true')
    closest_match_command.add_argument('-l', '--local_mash_and_info_file_prefix',
                                       help='the path prefix to the mash sketch file and corresponding info file')

    closest_match_command.add_argument('-mp', '--mash_path',
                                       help='path to the mash executable. If not provided it is assumed mash is in the PATH')

    filetype_extension = closest_match_command.add_mutually_exclusive_group(required=True)
    filetype_extension.add_argument('-f', '--fasta_file_pattern', help='pattern to match fasta files e.g "*.fas"')
    filetype_extension.add_argument('-fq', '--fastq_file_pattern', help='pattern to match fastq files e.g "*.fastq.gz"')
    filetype_extension.add_argument('-m', '--mash_sketch_file_pattern',
                                    help='pattern to match mash sketch files e.g "*.msh"')

    # make species_info sub command
    create_species_info_command = subparsers.add_parser('create_species_info',
                                                        help='Create species info TSV for locally created mash sketches'
                                                        )
    create_species_info_command.add_argument(
        '-m',
        '--mash_info_file',
        help='path to info file created using mash info -t',
        type=lambda x: is_valid_file(parser, x),
        required=True
    )
    create_species_info_command.add_argument(
        '-r',
        '--refseq_summary_file',
        help='path to refseq assembly summary file downloaded via wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt',
        type=lambda x: is_valid_file(parser, x),
        required=True
    )
    create_species_info_command.add_argument(
        '-b',
        '--bacsort_species_file',
        help='path to bacsort_species_definitions.txt',
        type=lambda x: is_valid_file(parser, x),
        required=True
    )
    create_species_info_command.add_argument(
        '-x',
        '--bacsort_excluded_assemblies_file',
        help='path to bacsort_excluded_assemblies.txt',
        type=lambda x: is_valid_file(parser, x),
        required=True
    )

    # info sub command
    info_command = subparsers.add_parser('info',
                                         help='Provide information about the data in bactinspector'
                                         )
    info_command.add_argument('-t', '--search_term',
                              help='search term to use when searching species within bactinspector', required=True)

    data_source = info_command.add_mutually_exclusive_group(required=True)
    data_source.add_argument('-s', '--summary', help='search the aggregate data', action='store_true')
    data_source.add_argument('-i', '--individual_records', help='search the individual refseq records',
                             action='store_true')

    args = parser.parse_args()
    return args


def choose_command(args):
    if args.command == 'check_species':
        run_check_species(args)
    elif args.command == 'closest_match':
        run_closest_match(args)
    elif args.command == 'create_species_info':
        run_create_species_info(args)
    elif args.command == 'info':
        if args.summary:
            run_info(args.search_term, 'summary')
        elif args.individual_records:
            run_info(args.search_term, 'individual_records')


def main():
    args = parse_arguments()
    choose_command(args)


if __name__ == "__main__":
    main()
