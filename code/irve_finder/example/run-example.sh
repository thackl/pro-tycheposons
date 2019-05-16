#!/usr/bin/env bash
set -e;

../irve_finder \
 --out example-run \
 --irve irve-v6.hmm \
 --irve-info irve-v6-info.tsv \
 --trna gorg-and-simons-tRNAs-nr99.fixHeader.ffn \
 --pre1 PRE1.hmm \
 --virus Phage_integrase.hmm \
 --faa MIT0604.faa \
 MIT0604.fna MIT0604.gff
