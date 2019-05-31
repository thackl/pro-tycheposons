library(tidyverse)
library(devtools)
library(thacklr)
library(patchwork)
library(gggenomes)
library(viridis)

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 9){
  stop("USAGE: Rscript cluster_finder.R <hit_file> <gene_coord.bed> <trna_hits.bed> <element_file> <gap.bed> <pre1.hits.bed> <virus.hits.tsv> <irve-info.tsv> <out.pdf>", call.=FALSE)
}


v5_hit_file <- args[1] #"test/MIT0604/pici-v5.hits.tsv"
gene_file <- args[2] #"test/MIT0604/MIT0604.cds.bed"
tRNA_file <- args[3] #"test/MIT0604/tRNA.hits.bed"
element_file <- args[4] #"test/MIT0604/scored_clusters.tsv"
gap_file <- args[5] #"test/MIT0604/MIT0604.gap.bed"
pre1_file <- args[6] #"test/MIT0604/pre1.hits.tsv"
virus_file <- args[7] #"test/MIT0604/virsorter.hits.tsv"
meta_file <- args[8]
cog_file <- args[9]
island_file <- args[10]
out_file <- args[11] #"test/MIT0604/clusters.pdf"

cluster_flank_length <- 3000
contigs_per_page <- 20

# gene coords
genes_0 <- read_tsv(gene_file, col_names=c("contig_id", "start", "end", "protein_id", "score", "strand")) %>%
  mutate(length = end - start + 1)

# elements
element_0 <- read_tsv(element_file) %>% filter(!secondary)

element_bounds <- element_0 %>%
  transmute(
    element_id,
    contig_id,
    cluster_score=score,
    archetype=upstreamType,
    ele_start = start-cluster_flank_length,
    ele_end = end+cluster_flank_length,
    ele_strand = strand
  )

element_1 <- element_0 %>%
  select(genome_id=element_id, contig_id, start, end, strand, score, secondary) %>%
  mutate(genome_id=str_replace_all(genome_id,"-","_"))

contig_lengths <- bind_rows(select(genes_0, contig_id, end)) %>%
  group_by(contig_id) %>%
  summarize(len=max(end)+500)

wrap_features_for_circular <- function(features, contig_lengths){
  features_1 <- left_join(features, select(contig_lengths,contig_id,len))
  features_2 <- bind_rows(
    mutate(features_1, start=start-len, end=end-len),
    features_1,
    mutate(features_1, start=start+len, end=end+len)
  )
  return(features_2)
}

# filter down to genes in/next to elements
# element_id => genome_id for plotting!
genes_1 <- wrap_features_for_circular(genes_0, contig_lengths)
genes_2 <- left_join(genes_1, element_bounds) %>%
  filter((start > ele_start & start < ele_end) |
           (end < ele_end & end > ele_start)) %>%
           mutate(genome_id=str_replace_all(element_id,"-","_"))

# hits
hits_0 <- read_tsv(v5_hit_file, col_names=c("protein_id","profile","hmmer_evalue","hmmer_score"), col_types="ccnn") %>%
  group_by(protein_id) %>% arrange(hmmer_evalue) %>%
  summarize_all(first) # best hit only

# virus hits
virus_hits_0 <- read_tsv(virus_file, col_names=c("protein_id","virus_profile","virus_evalue","virus_score"), col_types="ccnn") %>%
  group_by(protein_id) %>% arrange(-virus_score) %>%
  summarize_all(first) # best hit only

cogs_0 <- read_tsv(cog_file) %>% select(protein_id = gene_id, cog_num_strains)
islands_0 <- read_tsv(island_file)


# now clusters contain all gene - flanking, w/ and w/o hit
genes_3 <- left_join(genes_2, hits_0) %>%
  mutate(class = str_match(profile, "^([^_]+)")[,2]) %>%
  mutate(set = str_match(profile, "_([^-]+)$")[,2]) %>%
  mutate(profile_expr = str_replace(profile, "_[^-]+-?(.*)", "[\\1]")) %>%
  left_join(virus_hits_0, by=c("protein_id")) %>%
  left_join(cogs_0)


#print(genes_2)
#genes_2 %<>% split(.$element_id) %>%
#  keep(~any(str_detect(.x$profile, "alpA"), na.rm=TRUE)) %>%
#  bind_rows

## colors ----------------------------------------------------------------------
meta_0 <- read_tsv(meta_file)
profile_colors <- meta_0 %>% select(profile_id, plot_color) %>% deframe

## # focus on int, or others if int is missing
## classes <- unique(genes_2$class)
## class_colors <- c(
##   # Integration/Excision
##   int="#ffff33",
##   intY="#ffff33",
##   intS="#ffff33",
##   capsid="#4daf4a",
##   # Replication
##   pri_rep="#984ea3",
##   rep="#984ea3",
##   pri="#984ea3",
##   hel="#984ea3",
##   lig="#984ea3",
##   # Packaging/Interference
##   ter="brown2", #"#a65628",
##   terS="brown2", #"#a65628",
##   php="brown4",
##   xis="gold",
##   packaging="brown4",
##   pif="brown3",
##   # Regulation
##   stl="blue",
##   str="blue",
##   stl_str_rpr="blue",
##   rpr="blue"
## )
## # uncolored
## uncolored <- classes[!(classes %in% names(class_colors))]
## class_colors <- c(class_colors, rep("darkblue", length(uncolored)) %>% set_names(uncolored))

## focus & flip ----------------------------------------------------------------
# add evalue bin
genes_4 <- genes_3 %>%
  mutate(
    gene_focus = str_replace(element_id, "AG_(...)_", "AG-\\1-"),
    evalue_bin = cut(hmmer_evalue, c(Inf,1,1e-1,1e-3,1e-10,-Inf), label=F))

# dummy contigs, one per cluster
contigs_0 <- genes_4 %>% group_by(genome_id, contig_id) %>%
  summarize(
    length=max(end)+500,
    cluster_score=max(cluster_score, na.rm=T),
    seed_profile = (profile[protein_id == gene_focus])[1]
  ) %>%
    ungroup %>%
    arrange(-cluster_score)# %>%
#    mutate(plot_part = ceiling(1:n()/contigs_per_page))

islands_1 <- contigs_0 %>% select(genome_id, contig_id) %>%
  left_join(select(islands_0, -genome_id))

contig_lengths <- left_join(select(contigs_0, genome_id, contig_id),contig_lengths)
contig_ends <- bind_rows(
  select(contig_lengths, genome_id, contig_id) %>% mutate(start=-500, end=0, strand="+"),
  select(contig_lengths, genome_id, contig_id, end=len) %>% mutate(start=end-500, strand="+")
) %>%
  left_join(transmute(element_bounds, genome_id=str_replace_all(element_id,"-","_"), ele_start, ele_end)) %>%
  filter(start>=ele_start, end<=ele_end)

# order by score - add "$" as hotfix for reorder's partial matching
cluster_order <- contigs_0 %>% pull(genome_id) %>% paste0("$")

#------------------------------------------------------------------------------# tRNAs
aa3to1 <- function(aa3){
  aa1 <- read_tsv(
"aa_name	aa3	aa1
Alanine	Ala	A
Arginine	Arg	R
Asparagine	Asn	N
Aspartic acid	Asp	D
Cysteine	Cys	C
Glutamic acid	Glu	E
Glutamine	Gln	Q
Glycine	Gly	G
Histidine	His	H
Isoleucine	Ile	I
Leucine	Leu	L
Lysine	Lys	K
Methionine	Met	M
Phenylalanine	Phe	F
Proline	Pro	P
Serine	Ser	S
Threonine	Thr	T
Tryptophan	Trp	W
Tyrosine	Tyr	Y
Valine	Val	V
Selenocysteine	seC	X
") %>% select(-1) %>% deframe
  return(aa1[aa3])
}

expressify_tRNA <- function(tRNA_labels){
  str_split(tRNA_labels, "[-()]", simplify=T) %>%
    as_tibble() %>% mutate(V2 = aa3to1(V2)) %>%
    mutate(tRNA_exp = ifelse(V1 == "tmRNA", V1, paste0(V2,"[", V3, "]"))) %>%
    pull(tRNA_exp)
}

tRNAs_0 <- read_tsv(tRNA_file, col_names=c("contig_id", "start", "end", "type", "score", "strand", "full"))
tRNAs_1 <- wrap_features_for_circular(tRNAs_0, contig_lengths)
tRNAs_2 <- element_bounds %>% transmute(genome_id=str_replace_all(element_id,"-","_"), contig_id, ele_start, ele_end) %>%
  left_join(tRNAs_1, by=c("contig_id")) %>%
  filter(start>=ele_start, end<=ele_end) %>%
  mutate(
    expr = expressify_tRNA(type),
    strand=if_else(strand %in% c('+','-'), strand, '.')
  ) %>% unique # not sure why but some tRNAs get multipied

#------------------------------------------------------------------------------
# gaps
gaps_0 <- read_tsv(gap_file, col_names=c("contig_id", "start", "end", "gap_id", "score", "strand"))
gaps_1 <- tibble(contig_id=character(0), genome_id=character(0), start=numeric(0), end=numeric(0), strand=character(0))
if(nrow(gaps_0) > 0){
  gaps_1 <- wrap_features_for_circular(gaps_0, contig_lengths)
  gaps_1 <- element_bounds %>% transmute(genome_id=str_replace_all(element_id,"-","_"), contig_id, ele_start, ele_end) %>%
    left_join(gaps_1, by=c("contig_id")) %>%
    filter(start>=ele_start, end<=ele_end)
}

# flip - this impl only makes sense if we expect the integrase to be oriented inwards e.g. INT><PRI
# SaPIs are INT<>PRI
flip <- genes_4 %>%
    filter(protein_id == gene_focus & ele_strand == "-") %>%
    pull(genome_id) %>% unique %>% paste0("$")

pre1_0 <- read_tsv(pre1_file, col_names=c("contig_id", "start", "end", "pre_evalue", "strand"))
pre1_1 <- wrap_features_for_circular(pre1_0, contig_lengths)
pre1_2 <- element_bounds %>% transmute(genome_id=str_replace_all(element_id,"-","_"), contig_id, ele_start, ele_end) %>%
  left_join(pre1_1, by=c("contig_id")) %>%
  filter(start>=ele_start, end<=ele_end)

plot_contig_data <- function(contig_data, title, genomes_per_page=20){
  gg <- gggenomes(contig_data, genes_4, scale = list("lab", limits=c(genomes_per_page+.5,.5))) %>%
    add_features(tRNAs_2, "tRNAs") %>%
    add_features(element_1, "element") %>%
    add_features(gaps_1, "gaps") %>%
    add_features(contig_ends, "contig_ends") %>%
    add_features(pre1_2, "PRE1") %>%
    add_features(islands_1, "islands")

  # adjust layout
  gg <- gg %>%
    reorder(cluster_order) %>%
    flip(flip) %>%
    focus(protein_id==gene_focus, plus=50000, center="left", restrict_to_contig=FALSE)
  # plot geoms
  gg <- gg +
    ggtitle(title) +
    # islands
    geom_feature(aes(y=y+.25, yend=y+.25), data=use(islands), size=2, color="grey70") +
    geom_feature(data=use(gaps), size=1, color="black") +
    geom_feature(aes(y=y+.3, yend=y+.3), data=use(element, !secondary), size=.5, color="black") +
    #geom_feature(aes(y=y+.25, yend=y+.25), data=use(element, secondary), size=1, color="pink") +
    # profile source & evalue
    #geom_feature(aes(y=y+.40, yend=y+.40, alpha=virus_score), use(genes, !is.na(virus_score)), color="black", show.legend=FALSE, size=1) +
    #geom_point(aes((x+xend)/2, y+.3), use(genes, !is.na(set)), color="black", shape=21, size=2) +
    # genes
    geom_gene(aes(fill=profile), color="black", arrowhead_width=grid::unit(3,"mm"),
              arrowhead_height=grid::unit(3,"mm"),
              arrow_body_height=grid::unit(3,"mm"), show.legend=FALSE) +
    geom_feature(aes(y=y+.25, yend=y+.25, color=pmax(sqrt(cog_num_strains), sqrt(10))), use(genes), size=1) +
    #geom_text(aes((x+xend)/2-100, y=y-.3, label=str_extract(profile, "^[^_]+")), use(genes), angle=30, hjust=0, vjust=1, size=3) +
    geom_text(aes(pmin(x,xend)+0.2*abs(x-xend), y=y-.3, label=profile_expr), use(genes, !str_detect(profile, "^rep")), angle=30, hjust=0, vjust=1, size=3, parse=TRUE) +
    # score
    geom_text(aes(0, y+.45,label=format(cluster_score,digits=3)), use(contigs), hjust=1) +
    # contig start and end
    geom_segment(aes(x=(x+xend)/2-2, xend=(x+xend)/2-2, y=y-.4, yend=y+.4), use(contig_ends), linetype=1, color="royalblue3", size=.5) +
    # tRNAs and PRE1
    geom_feature(data=use(tRNAs, full=="partial"), size=6, color="grey30") +
    geom_feature(data=use(tRNAs, full=="full"), size=6, color="black") +
    geom_text(aes(pmin(x,xend)-100, y=y-.35, label=expr), use(tRNAs), angle=30, hjust=0, vjust=1, size=3, parse=TRUE) +
    #geom_feature(data=use(tRNAs, type=="tmRNA"), size=5, color="darkorchid") +
    geom_segment(aes(x=x, xend=x, y=y-.11, yend=y-.10), data=use(PRE1), arrow=arrow(length=grid::unit(1.5,"mm"), type="closed"), linejoin='mitre', size=.2, color="black") +
    scale_fill_manual(values=profile_colors, guide=FALSE) +
    scale_color_viridis("COG freq.", direction=-1, breaks=sqrt(c(0, 25, 100, 225, 400, 625)), labels=c(0, 25, 100, 225, 400, 625)) +
    expand_limits(cog_num_strains = c(10,623)) +
    #theme(legend.position="none") +
    coord_cartesian(xlim=c(-2500,35000))
  gg
}

plot_contig_data_pages <- function(contig_data, title, genomes_per_page=20){
  contig_data_1 <- contig_data %>%
    arrange(seed_profile, -cluster_score) %>%
    mutate(plot_part = ceiling(1:n()/genomes_per_page))
  plot_parts <- sort(unique(contig_data_1$plot_part))
  gg_pages <- purrr::map(plot_parts, function(pp){
    plot_contig_data(filter(contig_data_1,plot_part==pp),paste(title,pp), genomes_per_page=genomes_per_page)
  })
  return(gg_pages)
}

contigs_0 %<>% add_count(seed_profile, name="seed_n")
contigs_l <- contigs_0 %>% split(ifelse(contigs_0$seed_n >=5 & str_detect(contigs_0$seed_profile, "^(YR|LSR)"), contigs_0$seed_profile, "ungrouped"))

if(length(contigs_l) > 1){
  contigs_lg <- contigs_l[names(contigs_l) != "ungrouped"]
  contigs_l <- c(contigs_lg[rev(order(map_int(contigs_lg, nrow)))], contigs_l["ungrouped"])
}

pdf(out_file, title="IRVEs", width=20, height=15)
map2(contigs_l, paste(names(contigs_l),'(', map(contigs_l, nrow), ')'),  ~plot_contig_data_pages(.x, .y, genomes_per_page=20))

#plot_contig_data_pages(contigs_0, "IRVEs")
dev.off()
