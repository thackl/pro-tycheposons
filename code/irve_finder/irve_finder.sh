#!/usr/env/bin bash

PROTEINS_FA=$1
PROTEIN_POS=$2
HMM=$3
PROFILE_INFO=$4
PREFIX=$5

RSCRIPT_BIN=/nix/store/ymh756v7862r917ia1zyl73krvzd76rh-R-3.5.1-wrapper/bin/Rscript
SCRIPT_DIR=/tank/home/s216121/PostDoc/prochlorococcusIntegrase/analysis/pici
SHARED_DIR=/tank/home/s216121/PostDoc/prochlorococcusIntegrase/analysis/shared

if [ ! -f $HMM.h3i ]
then
    hmmpress $HMM
fi
hmmscan --cpu 4 --tblout $PREFIX.hmmscan.tbl --domtblout $PREFIX.hmmscan.domtbl -o $PREFIX.hmmscan.out $HMM $PROTEINS_FA
grep -v "^#" $PREFIX.hmmscan.tbl | perl -pe 's/ +/\t/g' | cut -f1,3,5,6 | tsv-uniq -f 2 >$PREFIX.hits.tsv

$RSCRIPT_BIN $SCRIPT_DIR/cluster_finder.R $PREFIX.hits.tsv $PROTEIN_POS $PROFILE_INFO $PREFIX.scored_clusters.tsv
$SHARED_DIR/pici-plot.MJA.R $PROTEIN_POS $PREFIX.scored_clusters.tsv $PREFIX.scored_clusters.pdf