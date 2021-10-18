#!/usr/bin/env Rscript
library(tidyverse)
library(magrittr)
library(argparser)

p <- arg_parser("Find IRVEs")
p %<>% add_argument("hit", "hits.bed")
p %<>% add_argument("pos", "CDS.bed")
p %<>% add_argument("meta", "irve-info")
p %<>% add_argument("out", "output .pdf")
args <- parse_args(p)

## argv <- c("irve-finder-aloha-np-v1/25m.pici-like.seqs-irve-raw-hits.tsv",
## "irve-finder-aloha-np-v1/25m.pici-like.seqs-cds.bed",
## "/home/thackl/Code/projects/irve-finder//data/irve-v6-info.tsv",
## "irve-finder-aloha-np-v1/25m.pici-like.seqs-irve-hits.tsv")
## args <- parse_args(p, argv)

hits_0 <- read_tsv(args$hit, col_names=c('protein_id','profile_id','hmmer_evalue','hmmer_score'))
hits_1 <- hits_0 %>% group_by(protein_id) %>%
  # top_n(1) returns multiple rows on tie, summarize to ensure single best
  top_n(1, hmmer_score) %>% summarize_all(first)
pos <- read_tsv(args$pos,c('contig_id','start','end','protein_id','score','strand')) %>% select(-score)
meta <- read_tsv(args$meta)

# keep low scoring hits if fragmented gene (multiple adjacent hits)
hitpos_0 <- left_join(hits_1,pos) %>% left_join(meta) %>%
  arrange(contig_id, start) %>%
  group_by(contig_id) %>%
  # look around and find same profile directly adjacent -> fragmented gene
  mutate(is_fragment = hmmer_score < profile_ga_score & (
           (profile_id == lag(profile_id, 1L, "") & abs(start - lag(end, 1L, Inf)) < 500) |
           (profile_id == lead(profile_id, 1L, "") & abs(end - lead(start, 1L, Inf)) < 500) |
           (profile_id == lag(profile_id, 2L, "") & abs(start - lag(end, 2L, Inf)) < 1000) |
           (profile_id == lead(profile_id, 2L, "") & abs(end - lead(start, 2L, Inf)) < 1000)))
hitpos_1 <- hitpos_0 %>%
  filter(hmmer_score >= .5 * profile_ga_score | is_fragment)

select(hitpos_1, 1:14, is_fragment) %>%
  write_tsv(args$out)

warnings()
