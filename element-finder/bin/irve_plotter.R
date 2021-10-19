library(argparser)
library(tidyverse)
library(devtools)
library(thacklr)
library(patchwork)
library(gggenomes)
library(viridis)
library(ggrepel)

p <- arg_parser("Plot IRVEs")
p %<>% add_argument("hits", "irve.bed")
p %<>% add_argument("gene", "CDS.bed")
p %<>% add_argument("att", "att.tsv")
p %<>% add_argument("trna", "tRNA.bed")
p %<>% add_argument("elements", "scored_clusters.tsv")
p %<>% add_argument("gaps", "gap.bed")
p %<>% add_argument("pre1", "pre1.tsv")
p %<>% add_argument("virus", "virus.tsv")
p %<>% add_argument("meta", "irve-info.tsv")
p %<>% add_argument("out", "output .pdf filename",  default="irves.pdf")
p %<>% add_argument("--cogs", "pro-623-cogs.tsv")
p %<>% add_argument("--islands", "pro-623-allmaps-islands-novt.tsv")
p %<>% add_argument("--cluster-flank-length", "size of region flanking the element", default = 3000)
p %<>% add_argument("--elements-per-page", "number of elements to plot per page", default = 20)

args <- parse_args(p)

# gene coords
genes_0 <- read_tsv(args$gene, col_names=c("contig_id", "start", "end", "protein_id", "score", "strand")) %>%
  mutate(length = end - start + 1)

# elements
element_0 <- read_tsv(args$elements) %>% filter(!secondary)

element_bounds <- element_0 %>%
  transmute(
    element_id,
    contig_id,
    cluster_score=score,
    archetype=upstreamType,
    ele_start = start-args$cluster_flank_length,
    ele_end = end+args$cluster_flank_length,
    ele_strand = strand
  )

element_1 <- element_0 %>%
  select(genome_id=element_id, contig_id, start, end, strand, score, secondary)

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
           mutate(genome_id=element_id)
# hits
hits_0 <- read_tsv(args$hits, col_names=c("protein_id","profile","hmmer_evalue","hmmer_score"), col_types="ccnn") %>%
  group_by(protein_id) %>% arrange(hmmer_evalue) %>%
  summarize_all(first) # best hit only

# virus hits
## virus_hits_0 <- read_tsv(args$virus, col_names=c("protein_id","virus_profile","virus_evalue","virus_score"), col_types="ccnn") %>%
##   group_by(protein_id) %>% arrange(-virus_score) %>%
##   summarize_all(first) # best hit only

if(!is.na(args$cogs))
  cogs_0 <- read_tsv(args$cogs) %>% select(protein_id = gene_id, cog_num_strains)

if(!is.na(args$islands))
  islands_0 <- read_tsv(args$islands)


# now clusters contain all gene - flanking, w/ and w/o hit
genes_3 <- left_join(genes_2, hits_0[0:4]) %>%
  mutate(class = str_match(profile, "^([^_]+)")[,2]) %>%
  mutate(set = str_match(profile, "_([^-]+)$")[,2]) %>%
  mutate(profile_expr = str_replace(profile, "_[^-]+-?(.*)", "[\\1]"))
 #%>% left_join(virus_hits_0, by=c("protein_id"))


#print(genes_2)
#genes_2 %<>% split(.$element_id) %>%
#  keep(~any(str_detect(.x$profile, "alpA"), na.rm=TRUE)) %>%
#  bind_rows

## colors ----------------------------------------------------------------------
meta_0 <- read_tsv(args$meta)
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
    gene_focus = element_id,
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

if(!is.na(args$islands)){
  islands_1 <- contigs_0 %>% select(genome_id, contig_id) %>%
    left_join(select(islands_0, -genome_id))
}

contig_lengths <- left_join(select(contigs_0, genome_id, contig_id),contig_lengths)
contig_ends <- bind_rows(
  select(contig_lengths, genome_id, contig_id) %>% mutate(start=-500, end=0, strand="+"),
  select(contig_lengths, genome_id, contig_id, end=len) %>% mutate(start=end-500, strand="+")
) %>%
  left_join(transmute(element_bounds, genome_id=element_id, ele_start, ele_end)) %>%
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

tRNAs_0 <- read_tsv(args$trna, col_names=c("contig_id", "start", "end", "type", "score", "strand", "full"))
tRNAs_1 <- wrap_features_for_circular(tRNAs_0, contig_lengths)
tRNAs_2 <- element_bounds %>% transmute(genome_id=element_id, contig_id, ele_start, ele_end) %>%
  left_join(tRNAs_1, by=c("contig_id")) %>%
  filter(start>=ele_start, end<=ele_end) %>%
  mutate(
    expr = expressify_tRNA(type),
    strand=if_else(strand %in% c('+','-'), strand, '.')
  ) %>% unique # not sure why but some tRNAs get multipied


# att-sites
weigh_att_scores <- function(x, y, w, df=50, xf=60000, yf=4000){
  ii <- seq_along(x)
  nw <- sapply(ii, function(i){
    d <- 1-abs(x-x[i])/df # dist to neighbor weights
    sum(w[d>0] * d[d>0])
  })
  # weigh by dist to element start
  nw <- nw * (1-abs(y)/yf) #* (1-abs(x-7500)/xf)
  ifelse(nw >= 0, nw, NA)
}

read_att_scores <-function(file, max_score = 400, min_score = 60, max_dist = 40){
  d0 <- read_blast(file) %>%
    mutate(protein_id = str_replace(saccver, "\\.[ud]s$", "")) %>%
    group_by(protein_id) %>%
    arrange(protein_id, sstart) %>%
    mutate(
      sdist = pmin(sstart,send) - lag(pmax(sstart,send),default=0),
      schain = cumsum(ifelse(sdist <= max_dist,0,1))) %>%
    arrange(protein_id, schain, qstart) %>%
    mutate(
      qdist = pmin(qstart,qend) - lag(pmax(qstart,qend),default=0),
      qchain = cumsum(ifelse(qdist <= max_dist,0,1))
    )

  d1 <- d0 %>%
    group_by(protein_id, schain, qchain) %>%
    summarize(
      hitnum = ifelse(n() > 4, 4, n()),
      dstart = min(sstart,send),
      dend = max(sstart,send),
      ustart = min(qstart,qend),
      uend = max(qstart,qend),
      att_score = sum(bitscore)) %>%
        filter(att_score < max_score) %>%
    group_by(protein_id) %>%
    #mutate(nwscore = weigh_att_scores((dstart+dend)/2, (ustart+uend)/2, bitscore)) %>%
    filter(att_score > min_score) %>%
    top_n(3, -ustart) %>%
    arrange(protein_id, -att_score) %>%
    mutate(att_color = row_number()) %>%
      select(protein_id, hitnum, att_score, att_color, ustart, uend,
             dstart, dend)
  d1
}

#d0 %>% filter(protein_id == "AG-355-I20_00922") %>% print(n=50)


#args <- list(att = "/home/thackl/Research/Pro-complex/IRVEs/irve-finder/results/irve-finder-v6-r7/pro-623-allmaps-att-m7.tsv")
#file <- args$att

att <- read_att_scores(args$att) %>%
  left_join(genes_4) %>%
  rename(int_start=start, int_end=end)


att_us <- att %>%
  mutate(
  start = ifelse(strand == "+", int_start - pmax(ustart, uend), int_end + pmin(ustart, uend)),
  end = ifelse(strand == "+", int_start - pmin(ustart, uend), int_end + pmax(ustart, uend)))
att_ds <- att %>% mutate(
  start = ifelse(strand == "+", int_end + pmin(dstart, dend), int_start - pmax(dstart, dend)),
  end = ifelse(strand == "+", int_end + pmax(dstart, dend), int_start - pmin(dstart, dend)))

filter(att_ds, protein_id == "MIT0604_1_1279") %>% glimpse

#------------------------------------------------------------------------------
# gaps
gaps_0 <- read_tsv(args$gaps, col_names=c("contig_id", "start", "end", "gap_id", "score", "strand"))
gaps_1 <- tibble(contig_id=character(0), genome_id=character(0), start=numeric(0), end=numeric(0), strand=character(0))
if(nrow(gaps_0) > 0){
  gaps_1 <- wrap_features_for_circular(gaps_0, contig_lengths)
  gaps_1 <- element_bounds %>% transmute(genome_id=element_id, contig_id, ele_start, ele_end) %>%
    left_join(gaps_1, by=c("contig_id")) %>%
    filter(start>=ele_start, end<=ele_end)
}

# flip - this impl only makes sense if we expect the integrase to be oriented inwards e.g. INT><PRI
# SaPIs are INT<>PRI
flip <- genes_4 %>%
    filter(protein_id == gene_focus & ele_strand == "-") %>%
    pull(genome_id) %>% unique %>% paste0("$")

## pre1_0 <- read_tsv(args$pre1, col_names=c("contig_id", "start", "end", "pre_evalue", "strand"))
## pre1_1 <- wrap_features_for_circular(pre1_0, contig_lengths)
## pre1_2 <- element_bounds %>% transmute(genome_id=element_id, contig_id, ele_start, ele_end) %>%
##   left_join(pre1_1, by=c("contig_id")) %>%
##   filter(start>=ele_start, end<=ele_end)

plot_contig_data <- function(contig_data, title, genomes_per_page=20){
  gg <- gggenomes(contig_data, genes_4, scale = list("lab", limits=c(genomes_per_page+.5,.5))) %>%
    add_features(tRNAs_2, "tRNAs") %>%
    add_features(element_1, "element") %>%
    add_features(gaps_1, "gaps") %>%
    add_features(contig_ends, "contig_ends") %>%
    #add_features(pre1_2, "PRE1") %>%
    add_features(att_us, "att_us") %>%
    add_features(att_ds, "att_ds")

  if(!is.na(args$islands))
    gg %<>% add_features(islands_1, "islands")

  # adjust layout
  gg <- gg %>%
    reorder(cluster_order) %>%
    flip(flip) %>%
    focus(protein_id==gene_focus, plus=50000, center="left", restrict_to_contig=FALSE)
  # plot geoms
  gg <- gg + ggtitle(title)
  # islands
  if(!is.na(args$islands))
    gg <- gg + geom_feature(aes(y=y+.25, yend=y+.25), data=use(islands), size=2, color="grey70")

  gg <- gg + geom_feature(data=use(gaps), size=1, color="black") +
    geom_feature(aes(y=y+.3, yend=y+.3), data=use(element, !secondary), size=.5, color="black") +
    #geom_feature(aes(y=y+.25, yend=y+.25), data=use(element, secondary), size=1, color="pink") +
    # profile source & evalue
    #geom_feature(aes(y=y+.40, yend=y+.40, alpha=virus_score), use(genes, !is.na(virus_score)), color="black", show.legend=FALSE, size=1) +
    #geom_point(aes((x+xend)/2, y+.2, shape=virus_profile), use(genes, !is.na(virus_score)), color="blue", size=2, show.legend=FALSE) +
    # genes
    geom_gene(aes(fill=profile), color="black", arrowhead_width=grid::unit(3,"mm"),
              arrowhead_height=grid::unit(3,"mm"),
              arrow_body_height=grid::unit(3,"mm"), show.legend=FALSE)
  if(!is.na(args$cogs))
    gg <- gg + geom_feature(aes(y=y+.25, yend=y+.25, color=pmax(sqrt(cog_num_strains), sqrt(10))), use(genes), size=1)

   #geom_text(aes((x+xend)/2-100, y=y-.3, label=str_extract(profile, "^[^_]+")), use(genes), angle=30, hjust=0, vjust=1, size=3) +
  gg <- gg +
    geom_text(aes(pmin(x,xend)+0.2*abs(x-xend), y=y-.3, label=profile_expr), use(genes, !str_detect(profile, "^rep")),
              angle=30, hjust=0, vjust=1, size=3, parse=TRUE) +
    # score
    geom_text(aes(0, y+.45,label=format(cluster_score,digits=3)), use(contigs), hjust=1) +
    # contig start and end
    geom_segment(aes(x=(x+xend)/2-2, xend=(x+xend)/2-2, y=y-.4, yend=y+.4), use(contig_ends), linetype=1, color="royalblue3", size=.5) +
    # tRNAs and PRE1
    geom_feature(data=use(tRNAs, full=="partial"), size=6, color="grey30") +
    geom_feature(data=use(tRNAs, full=="full"), size=6, color="black") +
    geom_text(aes(pmin(x,xend)-100, y=y-.35, label=expr), use(tRNAs), angle=30, hjust=0, vjust=1, size=3, parse=TRUE) +
    #geom_feature(data=use(tRNAs, type=="tmRNA"), size=5, color="darkorchid") +
#    geom_segment(aes(x=x, xend=x, y=y-.11, yend=y-.10), data=use(PRE1), arrow=arrow(length=grid::unit(1.5,"mm"), type="closed"), linejoin='mitre', size=.2, color="black") +
    # att
    geom_segment(aes(x=x, xend=xend, y=y-(att_color/10), yend=y-(att_color/10), color=att_color), data=use(att_us), size=4) +
    geom_segment(aes(x=x, xend=xend, y=y-(att_color/10), yend=y-(att_color/10), color=att_color), data=use(att_ds), size=4) +
    geom_text_repel(aes(x=x, y=y-.2, label=round(att_score)), data=use(att_ds), size=1.5, min.segment.length = 0.1) +
    # scales
    scale_fill_manual(values=profile_colors, guide=FALSE) +
    #scale_color_viridis("COG freq.", direction=-1, breaks=sqrt(c(0, 25, 100, 225, 400, 625)), labels=c(0, 25, 100, 225, 400, 625)) +
    scale_color_distiller(palette="Spectral", direction=1) +
    expand_limits(cog_num_strains = c(10,623)) +
    coord_cartesian(xlim=c(-2500,30000))
  gg
}

plot_contig_data_pages <- function(contig_data, title, genomes_per_page=20){
  contig_data_1 <- contig_data %>%
    arrange(seed_profile, -cluster_score) %>%
    mutate(plot_part = ceiling(1:n()/genomes_per_page))
  plot_parts <- sort(unique(contig_data_1$plot_part))
  print(contig_data_1)

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

pdf(args$out, title="IRVEs", width=20, height=15)
map2(contigs_l, paste(names(contigs_l),'(', map(contigs_l, nrow), ')'),  ~plot_contig_data_pages(.x, .y, genomes_per_page=args$elements_per_page))

#plot_contig_data_pages(contigs_0, "IRVEs")
dev.off()
