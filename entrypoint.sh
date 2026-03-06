#!/bin/bash

set -eo pipefail

cat - > /tmp/input.fasta

speciator fasta /tmp/input.fasta
