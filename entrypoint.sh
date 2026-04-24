#!/bin/bash

set -eo pipefail

cat - > /tmp/input.fasta

speciator -L info fasta /tmp/input.fasta
