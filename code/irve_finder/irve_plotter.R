library(tidyverse)
library(devtools)
library(thacklr)
library(patchwork)
library(gggenomes)

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3){
  stop("USAGE: Rscript cluster_finder.R <hit_file> <gene_coord.bed> <trna_hits.bed> <profile_info_file> <out_file>", call.=FALSE)
}

v5_hit_file <- args[1] #"test/MIT0604/pici-v5.hits.tsv"
gene_file <- args[2] #"test/MIT0604/MIT0604.cds.bed"
tRNA_file <- args[3] #"test/MIT0604/tRNA.hits.bed"
element_file <- args[4] #"test/MIT0604/scored_clusters.tsv"
gap_file <- args[5] #"test/MIT0604/MIT0604.gap.bed"
pre1_file <- args[6] #"test/MIT0604/pre1.hits.tsv"
virus_file <- args[7] #"test/MIT0604/virsorter.hits.tsv"
out_file <- args[8] #"test/MIT0604/clusters.pdf"

cluster_flank_length <- 1000
contigs_per_page <- 20

# gene coords
genes_0 <- read_tsv(gene_file, col_names=c("contig_id", "start", "end", "protein_id", "score", "strand")) %>%
  mutate(length = end - start + 1)
  
# elements
element_0 <- read_tsv(element_file) %>%
  filter(!secondary)

element_bounds <- element_0 %>%
  transmute(
    element_id=id,
    contig_id,
    cluster_score=score,
    archetype=upstreamType,
    ele_start = start-cluster_flank_length,
    ele_end = end+cluster_flank_length,
  )

element_1 <- element_0 %>%
  select(genome_id=id, contig_id, start, end, strand, score, secondary) %>%
  mutate(genome_id=str_replace_all(genome_id,"-","_"))

# filter down to genes in/next to elements
# element_id => genome_id for plotting!
genes_1 <- left_join(genes_0, element_bounds) %>%
  filter((start > ele_start & start < ele_end) |
           (end < ele_end & end > ele_start)) %>%
           mutate(genome_id=str_replace_all(element_id,"-","_"))

# hits
hits_0 <- read_tsv(v5_hit_file, col_names=c("protein_id","profile","hmmer_evalue","hmmer_score"))

# virus hits
virus_hits_0 <- read_tsv(virus_file, col_names=c("protein_id","virus_profile","virus_evalue","virus_score"))

# now clusters contain all gene - flanking, w/ and w/o hit
genes_2 <- left_join(genes_1, hits_0) %>%
  mutate(class = str_match(profile, "^([^-]+)")[,2]) %>%
  mutate(set = str_match(profile, "-([^-]+)$")[,2]) %>%
  left_join(virus_hits_0, by=c("protein_id"))

## colors ----------------------------------------------------------------------
# focus on int, or others if int is missing
classes <- unique(genes_2$class)
class_colors <- c(
  int="#ffff33",
  capsid="#4daf4a",
  pri_rep="#984ea3",
  rep="#984ea3",
  pri="#984ea3",
  ter="brown2", #"#a65628",
  terS="brown2", #"#a65628",
  php="brown4",
  xis="gold",
  packaging="brown4",
  pif="brown3",
  lig="#ff7f00",
  stl="blue",
  str="blue",
  stl_str_rpr="blue",
  rpr="blue"
)
# uncolored
uncolored <- classes[!(classes %in% names(class_colors))]
class_colors <- c(class_colors, rep("darkblue", length(uncolored)) %>% set_names(uncolored))

## focus & flip ----------------------------------------------------------------
class_order <- names(class_colors) #qc(int, cap, pri, hpA, reg, ter, hpB)

# add evalue bin
genes_3 <- genes_2 %>% 
  mutate(
    gene_focus = str_replace(element_id, "AG_(...)_", "AG-\\1-"),
    evalue_bin = cut(hmmer_evalue, c(Inf,1,1e-1,1e-3,1e-10,-Inf), label=F))

# dummy contigs, one per cluster
contigs_0 <- genes_3 %>% group_by(genome_id, contig_id) %>% 
  summarize(
    length=max(end)+500,
    cluster_score=max(cluster_score, na.rm=T)
  ) %>% ungroup %>%  arrange(-cluster_score)# %>%
#    mutate(plot_part = ceiling(1:n()/contigs_per_page))

# order by score - add "$" as hotfix for reorder's partial matching
cluster_order <- contigs_0 %>% pull(genome_id) %>% paste0("$")

#------------------------------------------------------------------------------# tRNAs
tRNAs_0 <- read_tsv(tRNA_file, col_names=c("contig_id", "start", "end", "type", "score", "strand", "full"))
tRNAs_1 <- contigs_0 %>% select(genome_id, contig_id) %>% unique %>%
  left_join(tRNAs_0, by=c("contig_id")) %>%
  mutate(type = str_replace(type, "tRNA-", ""), strand=if_else(strand %in% c('+','-'), strand, '.'))

#------------------------------------------------------------------------------
# gaps
gaps_0 <- read_tsv(gap_file, col_names=c("contig_id", "start", "end", "gap_id", "score", "strand"))
gaps_1 <- tibble(contig_id=character(0), genome_id=character(0), start=numeric(0), end=numeric(0), strand=character(0))
if(nrow(gaps_0) > 0){
  gaps_1 <- contigs_0 %>% select(genome_id, contig_id) %>% unique %>%
    left_join(gaps_0, by=c("contig_id"))
}

# flip - this impl only makes sense if we expect the integrase to be oriented inwards e.g. INT><PRI
# SaPIs are INT<>PRI
flip <- genes_3 %>%
    filter(protein_id == gene_focus & strand == "-") %>%
    pull(genome_id) %>% unique %>% paste0("$")

pre1_0 <- read_tsv(pre1_file, col_names=c("contig_id", "start", "end", "pre_evalue", "strand"))
pre1_1 <- contigs_0 %>% select(genome_id, contig_id) %>% unique %>% left_join(pre1_0)

plot_contig_data <- function(contig_data, title){
  gg <- gggenomes(contig_data, genes_3) %>% 
    add_features(tRNAs_1, "tRNAs") %>%
    add_features(element_1, "element") %>%
    add_features(gaps_1, "gaps") %>%
    add_features(pre1_1, "PRE1")
  # adjust layout
  gg <- gg %>%
    reorder(cluster_order) %>%
    flip(flip) %>%
    focus(protein_id==gene_focus, plus=50000, center="left")
  # plot geoms
  gg <- gg +
    ggtitle(title) +
    # islands
    #geom_feature(aes(y=y+.30, yend=y+.30), data=use(islands), size=2, color="burlywood2") +
    geom_feature(data=use(gaps), size=1.5, color="grey70") +
    geom_feature(aes(y=y+.30, yend=y+.30), data=use(element, !secondary), size=2, color="cornflowerblue") +
    geom_feature(aes(y=y+.30, yend=y+.30), data=use(element, secondary), size=2, color="pink") +
    # profile source & evalue
    geom_point(aes((x+xend)/2, y+.20, alpha=virus_score), use(genes, !is.na(virus_score)), color="black") +
    #geom_point(aes((x+xend)/2, y+.3), use(genes, !is.na(set)), color="black", shape=21, size=2) +
    # genes
    geom_gene(aes(fill=class), color="grey50", arrowhead_width=grid::unit(3,"mm"),
              arrowhead_height=grid::unit(3,"mm"),
              arrow_body_height=grid::unit(3,"mm")) +
    geom_text(aes((x+xend)/2, y=y-.3, label=profile), use(genes), angle=30, hjust=0, vjust=1, size=3) +
    # score
    geom_text(aes(-2500,y+.5,label=format(cluster_score,digits=3)), use(contigs)) +
    # tRNAs and PRE1
    geom_feature(data=use(tRNAs, full=="partial"), size=5, color="blue") +
    geom_feature(data=use(tRNAs, full=="full"), size=5, color="red") +
    geom_text(aes((x+xend)/2, y=y-.3, label=type), use(tRNAs), angle=30, hjust=0, vjust=1, size=3) +
    #geom_feature(data=use(tRNAs, type=="tmRNA"), size=5, color="darkorchid") +
    geom_segment(aes(x=x, xend=x, y=y-.40, yend=y-.10), data=use(PRE1), arrow=arrow(length=unit(.15,"cm"), type="closed"), linejoin='mitre', size=.5, color="darkgreen") +
    scale_fill_manual(values=class_colors) +
    scale_color_discrete(na.value=NA) +
    coord_cartesian(xlim=c(-2500,25000))# +
    #theme(legend.position="none")
  gg  
}

plot_contig_data_pages <- function(contig_data, title, genomes_per_page=20){
  contig_data_1 <- contig_data %>%  arrange(-cluster_score) %>%
    mutate(plot_part = ceiling(1:n()/genomes_per_page))
  plot_parts <- sort(unique(contig_data_1$plot_part))
  gg_pages <- purrr::map(plot_parts, function(pp){
    plot_contig_data(filter(contig_data_1,plot_part==pp),paste(title,pp))
  })
  return(gg_pages)
}

pdf(out_file, title="IRVEs", width=20, height=14)
plot_contig_data_pages(contigs_0, "IRVEs")
dev.off()
