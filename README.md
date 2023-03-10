# Speciator
## About
Speciator determines the species of an assembled genome FASTA by using [mash](https://github.com/marbl/Mash) to compare it to a hierarchical library of references.
Speciator is both fast and accurate and will continue to scale well as the number of references increases.

## Installing
Currently, the method for building the library is complex and not scripted.
If you wish to run speciator locally it's best to contact us for a copy of the Docker image.

## How it works
### Reference genomes
The references have been sourced from (1) Kleborate and (2) RefSeq.

### The algorithm
1. The query sequence is compiled into a mash profile.
2. It is then compared to a manually curated subset of genomes (the "Curated" library).
Currently, this is primarily sourced from Kleborate.
3. If a match is not found to one of the curated references, then it is searched against the "Genus" library, which contains a diverse set of representatives for each genus.
4. Provided a match is found to a genus, it is then compared to the comprehensive genus--specific library containing representatives of all species in that genus.
5. If no species is matched then it is searched against the library of representatives without full species classification.

### The match process
Rather than just take the best matching genome, Speciator uses the [Bactinspector](https://gitlab.com/antunderwood/bactinspector) software to select the consensus best match amongst the top close hits.
For some genera, where species classifications are unclear or they have been recently redefined for instance, there can be mislabelled RefSeq genomes (e.g. amongst the Klebsiella).
The consensus approach has proven robust against these kind of errors.

### Acknowledgements
The current version of speciator was written by Corin Yeats, and draws on the Bactinspector code written by Anthony Underwood.
Khalil Abudahab wrote the original version of speciator.

