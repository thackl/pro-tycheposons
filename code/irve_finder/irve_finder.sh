#!/usr/env/bin bash
set -e 

PROTEINS_FA=$1
REF_GFF=$2
GENOME_FA=$3
IRVE_HMM=$4
PROFILE_INFO=$5
VIR_HMM=$6
PRE_HMM=$7
TRNA_FA=$8
PREFIX=$9

RSCRIPT_BIN=/nix/store/nlmcmfzv9k8lvbwnlazjncwzafq0hj7f-R-3.5.3-wrapper/bin/Rscript
SCRIPT_DIR=/tank/home/s216121/PostDoc/prochlorococcusIntegrase/docs/IRVE-manuscript/code/irve_finder

echo "Extracting features from gff"
#cat $REF_GFF | ~/software/seq-scripts/bin/gff2bed --feature CDS >$PREFIX.cds.bed
#cat $REF_GFF | ~/software/seq-scripts/bin/gff2bed --feature gap >$PREFIX.gap.bed
#cat $REF_GFF | ~/software/seq-scripts/bin/gff2bed --feature tRNA >$PREFIX.tRNA.bed
#cat $REF_GFF | ~/software/seq-scripts/bin/gff2bed --feature tmRNA >>$PREFIX.tRNA.bed

echo "Searching for IRVE genes"
if [ ! -f $IRVE_HMM.h3i ]
then
    hmmpress $IRVE_HMM
fi
#hmmscan --cpu 4 --tblout $PREFIX.irve.hmmscan.tbl -o /dev/null $IRVE_HMM $PROTEINS_FA
#cat $PREFIX.irve.hmmscan.tbl | ~/software/seq-scripts/bin/hmmer-tbl2tsv | tail -n+3 | cut -f1,3,5,6 | tsv-uniq -f2 >$PREFIX.irve.hits.tsv


echo "Searching for viral genes"
if [ ! -f $VIR_HMM.h3i ]
then
    hmmpress $VIR_HMM
fi
#hmmscan --cpu 2 --tblout $PREFIX.virsorter.hmmscan.tbl -o /dev/null $VIR_HMM $PROTEINS_FA
cat $PREFIX.virsorter.hmmscan.tbl | ~/software/seq-scripts/bin/hmmer-tbl2tsv | tail -n+3 | cut -f1,3,5,6 | tsv-uniq -f2 >$PREFIX.virsorter.hits.tsv


echo "Searching for PRE1"
if [ ! -f $PRE_HMM.h3i ]
then
    hmmpress $PRE_HMM
fi
nhmmscan -o /dev/null --tblout $PREFIX.pre1.nhmmscan.tbl -E 1 $PRE_HMM $GENOME_FA
cat $PREFIX.pre1.nhmmscan.tbl | ~/software/seq-scripts/bin/hmmer-tbl2tsv | tail -n+3 | tsv-select -f 3,7,8,13,12 >$PREFIX.pre1.hits.tsv


echo "Searching for (partial) tRNAs"
if [ ! -f $TRNA_FA.nhr ]
then
    makeblastdb -dbtype nucl -in $TRNA_FA
fi
blastn -num_threads 6 -task blastn -db $TRNA_FA -query $GENOME_FA -reward 1 -penalty -1 -gapopen 2 -gapextend 1 -perc_identity 80 -evalue 10e-2 -outfmt 7 >$PREFIX.tRNA.blastn.tsv
cat $PREFIX.tRNA.blastn.tsv | ~/software/seq-scripts/bin/blast2bed -qa | bedtools merge -delim ";" -c 4,5,6,5 -o collapse,max,distinct,collapse |
  perl -ane '
    @scores = split(";", $F[6]);
    @trnas = split(";", $F[3]);
    @idx = sort{$scores[$b] <=> $scores[$a]}0..$#scores;
    @idx = grep{$scores[$_] == $scores[$idx[0]]}@idx; # max index with ties
    print join("\t", @F[0..2], join(",", do { my %seen; grep { !$seen{$_}++ } map {(split("\\|"))[2]} @trnas[@idx]}), @F[4,5], join(";", @scores[@idx])), "\n"' >$PREFIX.tRNA.blastn.bed
bedtools intersect -wb -v -f .9 -a $PREFIX.tRNA.blastn.bed -b $PREFIX.tRNA.bed >$PREFIX.tRNA.partial.bed
bedtools intersect -wb -f .9 -a $PREFIX.tRNA.blastn.bed -b $PREFIX.tRNA.bed >$PREFIX.tRNA.full.bed
cat <( cut -f 1-6 $PREFIX.tRNA.full.bed | perl -pe 's/$/\tfull/') <( cut -f1-6 $PREFIX.tRNA.partial.bed | perl -pe 's/$/\tpartial/') >$PREFIX.tRNA.hits.bed


echo "Detecting IRVE elements"
$RSCRIPT_BIN $SCRIPT_DIR/irve_finder.R $PREFIX.irve.hits.tsv $PREFIX.cds.bed $PREFIX.tRNA.hits.bed $PROFILE_INFO $PREFIX.scored_clusters.tsv
echo "Plotting results"
$RSCRIPT_BIN $SCRIPT_DIR/irve_plotter.R $PREFIX.irve.hits.tsv $PREFIX.cds.bed $PREFIX.tRNA.hits.bed $PREFIX.scored_clusters.tsv $PREFIX.gap.bed $PREFIX.pre1.hits.tsv $PREFIX.virsorter.hits.tsv $PREFIX.irve.pdf
