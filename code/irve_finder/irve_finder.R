library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3){
  stop("USAGE: Rscript cluster_finder.R <hit_file> <gene_coord_file> <profile_info_file> <out_file>", call.=FALSE)
}
hit_file <- args[1]
pos_file <- args[2]
meta_file <- args[3]
out_file <- args[4]
#hit_file <- "analysis/pici/finder_test/test.hits.tsv"
#pos_file <- "analysis/pici/finder_test/test.proteinPos.tsv"
#out_file <- "analysis/pici/finder_test/test.clusters.tsv"

max_dist <- 20000
min_size <- 2
max_initial_dist <- 4000

hits <- read_tsv(hit_file, col_names=c('profile','protein_id','hmmer_evalue','hmmer_score'))
pos <- read_tsv(pos_file,c('contig_id','protein_id','pos','len','strand'))
meta <- read_tsv(meta_file)

# try pre-filtering (keep all integrases, other genes only if there is another hit within max_initial_dist (e.g. 5000))
hitpos <- left_join(hits,pos) %>% left_join(meta)
hitpos %<>% group_by(contig_id) %>%
  arrange(contig_id,pos) %>%
  mutate(
    end=pos+len,
    prev_end=lag(end,default=-Inf),
    next_start=lead(pos,default=+Inf),
    closest=pmin(pos-prev_end,next_start-end)
  ) %>%
  filter(class=="int" | closest<=max_initial_dist) %>%
  ungroup %>%
  select(profile,protein_id,hmmer_evalue,hmmer_score,contig_id,pos,len,strand,class,score,length)

clusters <- hitpos %>%
    group_by(contig_id) %>%
    arrange(pos) %>%
    # mutate(diff=c(0,diff(pos)),large_gap=if_else(diff>max_dist,1,0),cluster_id=paste(contig_id,cumsum(large_gap),sep="_cls")) %>%
    mutate(end=pos+len,prev_end=lag(end,default=0),diff=pos-prev_end,large_gap=if_else(diff>max_dist,1,0),cluster_id=paste(contig_id,cumsum(large_gap),sep="_cls")) %>%
    add_count(cluster_id) %>%
    filter(n>=min_size) %>%
    select(-prev_end,-large_gap) %>%
    arrange(contig_id,cluster_id)

# subclustering by int (incl. direction) - not generic enough
#clusters %<>% group_by(cluster_id) %>%
#  mutate(
#    int = if_else(class=="int",1,0),
#    subcluster=cumsum(int),
#    intflip=if_else(int==1 & strand=="-",1,0),
#    subcluster=subcluster-intflip
#  ) %>%
#  ungroup %>%
#  mutate(cluster_id=paste(cluster_id,subcluster,sep="_"))

score <- function(cluster){
  score_sum <- cluster %>% group_by(class) %>% slice(which.max(lencor_score)) %>% pull(lencor_score) %>% sum
  cluster$cluster_score <- score_sum
  cluster
}

scored_clusters <- clusters %>%
  group_by(cluster_id) %>%
  mutate(
    length_deviation=1+abs(len/3-length)/length,
    lencor_score=score/length_deviation, # correct scores for length deviation
    diff=lead(diff), #remove first diff (not meaningful, distance to previous cluster)
    median_dist=median(diff,na.rm=TRUE),
    median_dist=if_else(median_dist<0,0.0,median_dist)
  ) %>%
  split(.$cluster_id) %>%
  map_df(score) %>%
  mutate(
    distcor_cluster_score=cluster_score*(max_dist-median_dist)/max_dist
  ) %>%
  arrange(-distcor_cluster_score, -n, cluster_id, pos)

write_tsv(scored_clusters,out_file)
