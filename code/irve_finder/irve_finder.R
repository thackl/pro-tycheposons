library(tidyverse)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 5){
  stop("USAGE: Rscript cluster_finder.R <hit_file> <gene_coord.bed> <trna_hits.bed> <profile_info_file> <out_file>", call.=FALSE)
}
hit_file <- args[1]
pos_file <- args[2]
trna_file <- args[3]
meta_file <- args[4]
out_file <- args[5]
interpret_contig_as_circular_if_smaller_than <- 0
if(length(args) > 5){
  interpret_contig_as_circular_if_smaller_than <- args[6]
}

hits_0 <- read_tsv(hit_file, col_names=c('protein_id','profile_id','hmmer_evalue','hmmer_score'))
hits_1 <- hits_0 %>% group_by(protein_id) %>%
  # top_n(1) returns multiple rows on tie, summarize to ensure single best
  top_n(1, hmmer_score) %>% summarize_all(first)
pos <- read_tsv(pos_file,c('contig_id','start','end','protein_id','score','strand')) %>% select(-score)
meta <- read_tsv(meta_file)
tRNAs <- read_tsv(trna_file, col_names=c('contig_id', 'start', 'end', 'type', 'score', 'strand', 'full'))
contig_lengths <- bind_rows(select(pos, contig_id, end), select(tRNAs, contig_id, end)) %>%
  group_by(contig_id) %>%
  summarize(len=max(end)+500)

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
hitpos <- hitpos_0 %>%
  filter(hmmer_score >= .5 * profile_ga_score | is_fragment)


integrase_to_element <- function(proteinId, hitpos, contig_lengths, tRNAs, targetFunctions, upstreamRange=1000, downstreamRange=25000, minLen=5000){
  integrase <- filter(hitpos, protein_id==proteinId) %>%
    arrange(hmmer_evalue) %>% head(n=1) # top hit per int
  score <- 0
  upstreamType <- "None"
  downstreamType <- "None"
  contig_len <- (filter(contig_lengths, contig_id == integrase$contig_id))$len
  hits_in_range <- tibble()
  trna_on_contig <- tRNAs %>% filter(contig_id == integrase$contig_id)
  hits_on_contig <- hitpos %>% filter(contig_id == integrase$contig_id)
  is_circular <- contig_len < interpret_contig_as_circular_if_smaller_than
  if(is_circular){
    trna_on_contig <- bind_rows(
      mutate(trna_on_contig, start=start-contig_len, end=end-contig_len),
      trna_on_contig,
      mutate(trna_on_contig, start=start+contig_len, end=end+contig_len)
    )
    hits_on_contig <- bind_rows(
      mutate(hits_on_contig, start=start-contig_len, end=end-contig_len),
      hits_on_contig,
      mutate(hits_on_contig, start=start+contig_len, end=end+contig_len)
    )
  }
  direction <- integrase$strand
  # search for closest tRNA within 1kbp (if downstream - flip direction)
  if(integrase$strand == "+"){
    downstream_trna <- tRNAs %>% filter(contig_id == integrase$contig_id, start>=integrase$start, end<=integrase$end+1000) %>% top_n(-1, start)
    if(nrow(downstream_trna)){
      direction <- "-"
    }
  } else {
    downstream_trna <- tRNAs %>% filter(contig_id==integrase$contig_id, start>=integrase$start-1000, end<=integrase$end) %>% top_n(1, start)
    if(nrow(downstream_trna)){
      direction <- "+"
    }
  }
  if(direction == "+"){
    upstream_trna <- trna_on_contig %>% filter(start>=integrase$start-upstreamRange, end<=integrase$start) %>% top_n(1, start)
    hits_in_range <- hits_on_contig %>% filter(end>=integrase$start-upstreamRange, end<=integrase$end+downstreamRange)
    if(nrow(upstream_trna)!=0){
      startPos <- upstream_trna$start
      upstreamType <- str_c(upstream_trna$type,"-",upstream_trna$full)
      score <- score + 5
    } else {
      startPos <- min(c(integrase$start, hits_in_range %>% filter(!(`function` %in% targetFunctions)) %>% pull(start)))
    }
    endPos <- integrase$end + minLen
    downstream_trna <- tRNAs %>% filter(contig_id == integrase$contig_id, start>=integrase$start, end<=integrase$end+downstreamRange) %>% top_n(-1, start)
    # only consider downstream_trna if upstream_trna was found as well
    if(nrow(upstream_trna)!=0 && nrow(downstream_trna)!=0){
      endPos <- downstream_trna$end
      downstreamType <- str_c(downstream_trna$type,"-",downstream_trna$full)
      score <- score + 2
    } else {
      endPos <- max(endPos, max(hits_in_range$end))
    }
  } else {
    upstream_trna <- trna_on_contig %>% filter(start>=integrase$end, end<=integrase$end+upstreamRange) %>% top_n(-1, start)
    hits_in_range <- hits_on_contig %>% filter(start>=integrase$start-downstreamRange, start<=integrase$end+upstreamRange)
    if(nrow(upstream_trna)!=0){
      endPos <- upstream_trna$end
      upstreamType <- str_c(upstream_trna$type,"-",upstream_trna$full)
      score <- score + 5
    } else {
      endPos <- max(c(integrase$end, hits_in_range %>% filter(!(`function` %in% targetFunctions)) %>% pull(end)))
    }
    startPos <- integrase$start - minLen
    downstream_trna <- tRNAs %>% filter(contig_id==integrase$contig_id, start>=integrase$start-downstreamRange, end<=integrase$end) %>% top_n(1, start)
    if(nrow(upstream_trna)!=0 && nrow(downstream_trna)!=0){
      startPos <- downstream_trna$start
      downstreamType <- str_c(downstream_trna$type,"-",downstream_trna$full)
      score <- score + 2
    } else {
      startPos <- min(startPos, min(hits_in_range$start))
    }
  }
  if(!is_circular){
    startPos <- max(startPos, 0)
    endPos <- min(endPos, contig_len)
  } else {
    # make sure the element is no longer than the contig
    if(direction == "+"){
      endPos <- min(endPos, startPos + contig_len)
    } else {
      startPos <- max(startPos, endPos - contig_len)
    }
  }
  hits_in_range <- hits_on_contig %>% filter(end>=startPos, start<=endPos)
  # TODO: use scores from meta
  score <- score + length(unique(hits_in_range$protein_id)) + length(unique(hits_in_range$`function`))
  return(tibble(id=integrase$protein_id, start=startPos, upstreamType, end=endPos, downstreamType, strand=direction, score, contig_id = integrase$contig_id))
}

flag_lower_scoring_overlap <- function(clusters){
  require(plyranges)
  loosers <- clusters %>%
    as_granges(seqnames=contig_id) %>%
    join_overlap_self(maxgap=-1L, minoverlap=100L) %>%
    filter(id < id.overlap) %>%
    mutate(looser = if_else(score<score.overlap, id, id.overlap)) %>%
    as_tibble %>%
    pull(looser) %>%
    unique
  return(clusters$id %in% loosers)
}

target_functions <- c("Tyrosine Recombinase", "Large Serine Recombinase", "Serine Recombinase")
int_ids <- hitpos %>%
  filter(`function` %in% target_functions & irve_score > 0) %>%
  pull(protein_id) %>% unique

scored_clusters_0 <- map_df(int_ids,integrase_to_element,hitpos,contig_lengths,tRNAs,target_functions)

scored_clusters_1 <- scored_clusters_0 %>%
  mutate(secondary=flag_lower_scoring_overlap(.)) %>%
  arrange(secondary,-score)

get_no_int_clusters <- function(hitpos, existing_clusters, max_dist=5000, min_size=3){
  # filter all hits that are on any existing cluster
  already_on_cluster <- join_overlap_inner(hitpos %>% as_granges(seqnames=contig_id), existing_clusters %>% as_granges(seqnames=contig_id)) %>%
    as_tibble %>%
    pull(protein_id) %>%
    unique

  remaining_hits <- hitpos %>% filter(! protein_id %in% already_on_cluster)

  no_int_clusters <- remaining_hits %>%
    group_by(contig_id) %>%
    arrange(contig_id,start) %>%
    mutate(
      prev_end=lag(end,default=-Inf),
      diff=start-prev_end,
      large_gap=if_else(diff>max_dist,1,0),
      cluster_id=paste(contig_id,cumsum(large_gap),sep="_cls")
    ) %>%
    add_count(cluster_id) %>%
    filter(n>=min_size)

  scored_no_int_clusters <- no_int_clusters %>%
    group_by(cluster_id) %>%
    summarize(
      id=dplyr::first(protein_id, order_by=start),
      start=min(start),
      upstreamType="no_int", end=max(end),
      downstreamType="no_int",
      strand="+",
      score=length(module)+length(unique(module)),
      contig_id=dplyr::first(contig_id),
      secondary=FALSE
    ) %>%
    ungroup %>%
    select(-cluster_id)

  return(scored_no_int_clusters)
}

no_int_clusters <- get_no_int_clusters(hitpos, scored_clusters_1)

scored_clusters_2 <- bind_rows(scored_clusters_1, no_int_clusters)

write_tsv(scored_clusters_2,out_file)

warnings()
