#!/usr/env/bin bash

CONTIGS=$1
HMM=$2
PROFILE_INFO=$3
PREFIX=$4

SCRIPT_DIR=/tank/home/s216121/PostDoc/prochlorococcusIntegrase/analysis/pici

prodigal -a $PREFIX.proteins.faa -i $CONTIGS -o $PREFIX.proteins.gb -p meta -m
grep "^>" $PREFIX.proteins.faa | cut -f1,3,5,7 -d" " | perl -pe 's/^>(\S+)_(\d+) /\1\t\1_\2\t/;s/ /\t/g;' >$PREFIX.proteins.tsv
perl -ane '$s="+";$s="-" if($F[4]=~/-/);print "$F[0]\t$F[1]\t$F[2]\t".($F[3]-$F[2]+1)."\t$s\n"' $PREFIX.proteins.tsv >$PREFIX.gene-coord.tsv

bash $SCRIPT_DIR/cluster_finder.sh $PREFIX.proteins.faa $PREFIX.gene-coord.tsv $HMM $PROFILE_INFO $PREFIX