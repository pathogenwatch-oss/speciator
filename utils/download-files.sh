#!/usr/bin/env bash

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
mv assembly_summary.txt bacteria_assembly_summary.txt

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt
mv assembly_summary.txt viral_assembly_summary.txt

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt
mv assembly_summary.txt fungi_assembly_summary.txt

grep "^[^#;]" bacteria_assembly_summary.txt \
 | awk -F '\t' '{print $1, $20}' | sed 's/ftp:/rsync:/' | sed 's&\(.*\)/\(.*\)$&\2 \1/\2/\2_genomic.fna.gz&' \
  > bacteria_links.txt

grep "^[^#;]" viral_assembly_summary.txt \
 | awk -F '\t' '{print $1, $20}' | sed 's/ftp:/rsync:/' | sed 's&\(.*\)/\(.*\)$&\2 \1/\2/\2_genomic.fna.gz&' \
  > viral_links.txt

grep "^[^#;]" fungi_assembly_summary.txt \
 | awk -F '\t' '{print $1, $20}' | sed 's/ftp:/rsync:/' | sed 's&\(.*\)/\(.*\)$&\2 \1/\2/\2_genomic.fna.gz&' \
  > fungi_links.txt

mkdir -p assemblies

# Needs a step to remove any assemblies already present that aren't in the links.

cut -d ' ' -f 3 fungi_links.txt | sed 's/rsync:\/\/ftp\.ncbi\.nlm\.nih.gov\/genomes\///' > fungi_paths.txt

cut -d ' ' -f 3 viral_links.txt | sed 's/rsync:\/\/ftp\.ncbi\.nlm\.nih.gov\/genomes\///' > viral_paths.txt

cut -d ' ' -f 3 bacteria_links.txt | sed 's/rsync:\/\/ftp\.ncbi\.nlm\.nih.gov\/genomes\///' > bacteria_paths.txt

cat fungi_paths viral_paths bacteria_paths > all_paths

rsync -av --no-relative --files-from=all_paths.txt rsync://ftp.ncbi.nlm.nih.gov/genomes assemblies/
#rsync -av --no-relative --files-from=fungi_paths.txt rsync://ftp.ncbi.nlm.nih.gov/genomes assemblies/
#rsync -av --no-relative --files-from=viral_paths.txt rsync://ftp.ncbi.nlm.nih.gov/genomes assemblies/
#rsync -av --no-relative --files-from=bacteria_paths.txt rsync://ftp.ncbi.nlm.nih.gov/genomes assemblies/

find . -maxdepth 1 -name "*.gz" | xargs -I xxx basename xxx | xargs -P 4 -I xxx mash sketch xxx

grep "^[^#;]" assembly_summary.txt | awk -F '\t' '{print $6}' | sort -u > organism_taxids.lst

taxonkit lineage organism_taxids.lst | awk '$2!=""' > organism_lineages.txt
cat organism_lineages.txt \
  | taxonkit reformat \
  | csvtk -H cut -f 1,3 -t \
  | csvtk -H sep -f 2 -s ';' -R -t \
  | csvtk add-header -t -n taxonId,superkingdom,phylum,class,order,family,genus,species \
  | csvtk tab2csv \
  | csvtk cut -f 1,2,7,8 > taxon_names.csv

cat organism_lineages.txt \
  | taxonkit reformat -t \
  | csvtk -H cut -f 1,4 -t \
  | csvtk -H sep -f 2 -s ';' -R -t \
  | csvtk add-header -t -n taxonId,superkingdom,phylum,class,order,family,genus,species \
  | csvtk tab2csv \
  | csvtk cut -f 1,2,7,8 > taxon_codes.csv