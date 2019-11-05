#!/usr/bin/env Rscript
library(tidyverse)
library(thacklr)
library(argparser)

p <- arg_parser("Find IRVEs")
p %<>% add_argument("hit", "hits.bed")
p %<>% add_argument("trna", "tRNA.bed")
p %<>% add_argument("att", "att-blast.tsv")
p %<>% add_argument("contig-length", "contig-length.tsv")
p %<>% add_argument("out", "output .pdf")
p %<>% add_argument("--circle-max", "interpret contig as circular if smaller than this", default=0)
args <- parse_args(p)

# att-sites
weigh_att_scores <- function(x, y, w, df=100, xf=60000, yf=4000){
  ii <- seq_along(x)
  nw <- sapply(ii, function(i){
    d <- 1-abs(x-x[i])/df # dist to neighbor weights
    sum(w[d>0] * d[d>0])
  })
  # weigh by dist to element start
  nw <- nw * (1-abs(y)/yf) #* (1-abs(x-7500)/xf)
  ifelse(nw >= 0, nw, NA)
}

read_att_scores <-function(file, max_bitscore = 200, min_nwscore = 50){
  read_blast(file) %>%
    mutate(protein_id = str_replace(saccver, "\\.[ud]s$", "")) %>%
    group_by(protein_id) %>%
    filter(bitscore < max_bitscore) %>%
    mutate(nwscore = weigh_att_scores((sstart+send)/2, (qstart+qend)/2,
             bitscore)) %>%
    filter(nwscore > min_nwscore) %>%
    select(protein_id, nwscore, ustart=qstart, uend=qend, dstart=sstart,
           dend=send)
}

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
  is_circular <- contig_len < args$circle_max
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
  score <- score +
    sum(filter(hits_in_range, !is_fragment) %>% pull(profile_irve_score)) +
    sum(filter(hits_in_range, is_fragment) %>% pull(profile_irve_score)) * .1 +
    length(unique(hits_in_range$`function`))
  return(tibble(
    contig_id = integrase$contig_id, start=startPos, end=endPos,
    element_id=integrase$protein_id, score, strand=direction,
    upstreamType, downstreamType, int_profile_id=integrase$profile_id,
    profile_ids=paste(hits_in_range$profile_id, collapse=",")))
}

flag_lower_scoring_overlap <- function(clusters){
  require(plyranges)
  loosers <- clusters %>%
    as_granges(seqnames=contig_id) %>%
    join_overlap_self(maxgap=-1L, minoverlap=100L) %>%
    filter(element_id < element_id.overlap) %>%
    mutate(looser = if_else(score<score.overlap, element_id, element_id.overlap)) %>%
    as_tibble %>%
    pull(looser) %>%
    unique
  return(clusters$element_id %in% loosers)
}

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
      contig_id=dplyr::first(contig_id),
      start=min(start),
      end=max(end),
      element_id=dplyr::first(protein_id, order_by=start),
      score=length(module)+length(unique(module)),
      strand="+",
      upstreamType="no_int",
      downstreamType="no_int",
      int_profile_id=NA,
      profile_ids=paste(profile_id, collapse=","),
      secondary=FALSE
    ) %>%
    ungroup %>%
    select(-cluster_id)

  return(scored_no_int_clusters)
}


#setwd("../irve-finder-20190807-204501")
#args <- parse_args(p, c("MIT0604.irve.hits.tsv", "MIT0604.tRNA.hits.bed", "MIT0604.att-m7.tsv", "MIT0604.contig-length.tsv", "bar.tsv"))

hitpos <- read_tsv(args$hit)
tRNAs <- read_tsv(args$trna, col_names=c('contig_id', 'start', 'end', 'type', 'score', 'strand', 'full'))
contig_lengths <- read_tsv(args$contig_length, col_names=c('contig_id', 'len'))
#contig_lengths <- bind_rows(select(pos, contig_id, end), select(tRNAs, contig_id, end)) %>%
#  group_by(contig_id) %>%
#  summarize(len=max(end)+500)
att <- read_att_scores(args$att)

target_functions <- c("Tyrosine Recombinase", "Large Serine Recombinase", "Serine Recombinase")
int_ids <- hitpos %>%
  filter(`function` %in% target_functions & profile_irve_score > 0) %>%
  pull(protein_id) %>% unique

att_sites <- function(){
upstreamRange=1000
downstreamRange=25000
int_id <- int_ids[4]
integrase <- filter(hitpos, protein_id==int_id) %>%
    arrange(hmmer_evalue) %>% head(n=1) # top hit per int

hits_in_range <- if(integrase$strand[1] ==  "+"){
  hitpos %>% filter(contig_id == integrase$contig_id, end>=integrase$start-upstreamRange, end<=integrase$end+downstreamRange)
}else{
  hitpos %>% filter(contig_id == integrase$contig_id, start>=integrase$start-downstreamRange, start<=integrase$end+upstreamRange)
}

hits_in_range %>% select(1,2,start,end)
filter(att, protein_id==int_id) %>%
  mutate(
    ustart = integrase$start[1] - min(2000,integrase$start[1]) + ustart,
    uend = integrase$start[1] - min(2000,integrase$start[1]) + uend,
    dstart = integrase$end[1] + dstart,
    dend = integrase$end[1] + dend,
  size = dend-ustart)
}

scored_clusters_0 <- map_df(int_ids,integrase_to_element,hitpos,contig_lengths,tRNAs,target_functions)
att %>% glimpse
scored_clusters_1 <- scored_clusters_0 %>%
  mutate(secondary=flag_lower_scoring_overlap(.)) %>%
  #  arrange(secondary,-score)
  filter(!secondary) %>% arrange(-score)

no_int_clusters <- get_no_int_clusters(hitpos, scored_clusters_1)

scored_clusters_2 <- bind_rows(scored_clusters_1, no_int_clusters)

write_tsv(scored_clusters_2,args$out)

warnings()


hitpos
