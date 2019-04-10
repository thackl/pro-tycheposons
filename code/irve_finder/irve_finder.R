library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3){
  stop("USAGE: Rscript cluster_finder.R <hit_file> <gene_coord.bed> <trna_hits.bed> <profile_info_file> <out_file>", call.=FALSE)
}
hit_file <- args[1]
pos_file <- args[2]
trna_file <- args[3]
meta_file <- args[4]
out_file <- args[5]

hits <- read_tsv(hit_file, col_names=c('protein_id','profile','hmmer_evalue','hmmer_score'))
pos <- read_tsv(pos_file,c('contig_id','start','end','protein_id','score','strand')) %>% select(-score)
meta <- read_tsv(meta_file)
tRNAs <- read_tsv(trna_file, col_names=c('contig_id', 'start', 'end', 'type', 'score', 'strand', 'full'))

hitpos <- left_join(hits,pos) %>% left_join(meta)

integrase_to_element <- function(proteinId, hitpos, tRNAs, upstreamRange=1000, downstreamRange=25000, minLen=5000){
  integrase <- filter(hitpos, protein_id==proteinId)
  score <- 0
  upstreamType <- "None"
  downstreamType <- "None"
  hits_in_range <- tibble()
  if(integrase$strand == "+"){
    upstream_trna <- tRNAs %>% filter(contig_id == integrase$contig_id, start>=integrase$start-upstreamRange, end<=integrase$start) %>% top_n(1, start)
    hits_in_range <- hitpos %>% filter(contig_id == integrase$contig_id, end>=integrase$start-upstreamRange, end<=integrase$end+downstreamRange)
    if(nrow(upstream_trna)!=0){
      startPos <- upstream_trna$start
      upstreamType <- str_c(upstream_trna$type,"-",upstream_trna$full)
      score <- score + 5
    } else {
      startPos <- min(c(integrase$start, hits_in_range %>% filter(class != "int") %>% pull(start)))
    }
    endPos <- integrase$end + minLen
    downstream_trna <- tRNAs %>% filter(contig_id == integrase$contig_id, start>=integrase$start, end<=integrase$end+downstreamRange) %>% top_n(-1, start)
    # only consider downstream_trna if upsteam_trna was found as well
    if(nrow(upsteam_trna)!=0 && nrow(downstream_trna)!=0){
      endPos <- downstream_trna$end
      downstreamType <- str_c(downstream_trna$type,"-",downstream_trna$full)
      score <- score + 2
    } else {
      endPos <- max(endPos, max(hits_in_range$end))
    }
  } else {
    upstream_trna <- tRNAs %>% filter(contig_id == integrase$contig_id, start>=integrase$end, end<=integrase$end+upstreamRange) %>% top_n(-1, start)
    hits_in_range <- hitpos %>% filter(contig_id == integrase$contig_id, start>=integrase$start-downstreamRange, start<=integrase$end+upstreamRange)
    if(nrow(upstream_trna)!=0){
      endPos <- upstream_trna$end
      upstreamType <- str_c(upstream_trna$type,"-",upstream_trna$full)
      score <- score + 5
    } else {
      endPos <- max(c(integrase$end, hits_in_range %>% filter(class != "int") %>% pull(end)))
    }
    startPos <- integrase$start - minLen
    downstream_trna <- tRNAs %>% filter(contig_id==integrase$contig_id, start>=integrase$start-downstreamRange, end<=integrase$end) %>% top_n(1, start)
    if(nrow(upsteam_trna)!=0 && nrow(downstream_trna)!=0){
      startPos <- downstream_trna$start
      downstreamType <- str_c(downstream_trna$type,"-",downstream_trna$full)
      score <- score + 2
    } else {
      startPos <- min(startPos, min(hits_in_range$start))
    }
  }
  hits_in_range <- hitpos %>% filter(contig_id == integrase$contig_id, end>=startPos, start<=endPos)
  score <- score + length(unique(hits_in_range$protein_id)) + length(unique(hits_in_range$class))
  return(tibble(id=integrase$protein_id, start=startPos, upstreamType, end=endPos, downstreamType, strand=integrase$strand, score, contig_id = integrase$contig_id))
}

flag_lower_scoring_overlap <- function(clusters){
  require(plyranges)
  loosers <- clusters %>%
    as_granges(seqnames=contig_id) %>%
    join_overlap_self(maxgap=-1L, minoverlap=100L) %>%
    filter(id != id.overlap, score != score.overlap) %>%
    mutate(looser = if_else(score<score.overlap, id, id.overlap)) %>%
    as_tibble %>%
    pull(looser) %>%
    unique
  return(clusters$id %in% loosers)
}

scored_clusters <- hitpos %>%
  filter(class=="int") %>%
  pull(protein_id) %>%
  map_df(integrase_to_element,hitpos,tRNAs) %>%
  mutate(secondary=flag_lower_scoring_overlap(.)) %>%
  arrange(secondary,-score)
write_tsv(scored_clusters,out_file)

