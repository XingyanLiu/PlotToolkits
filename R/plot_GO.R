library("ggplot2")
library("dplyr")
library("stringr")
library("forcats")
library("clusterProfiler")
library("enrichplot")

setwd("E:/lxy_pro/003")
source("RFunxGO.R")

allStages = c(
  "2cell", "4cell", "8cell", "32cell", "256cell", "B",
  "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0", "L2"
)

###############################################################
fdir = "E:/lxy_pro/003/amph/GOresults/1119-amphRef/GO"
figdir = "E:/lxy_pro/003/amph/GOresults/1119-amphRef/figs"
dir.create(figdir, recursive = T)

"======= loding GO tables into a list ========" %>% message()
i = 2
GOTables = list()
for (i in seq(allStages)){
  stg1 = allStages[i]
  stg2 = allStages[i + 1]
  tag0 = sprintf("%s_%s", stg1, stg2)
  fni = sprintf("GO_%s.csv", tag0)
  fpi = sprintf("%s/%s", fdir, fni)
  if (file.exists(fpi)){
    gotable = read.csv(fpi, as.is = T)
    # gotable$group = stg
    GOTables[[stg2]] = gotable
    sprintf("file loaded (%s)", fpi) %>% message()
  }else{
    sprintf("file not exists, skipped (%s)", fpi) %>% message()
  }
}

###############################################################################
"making some modifications on the GO tables, including filtering and selection"
for (nm in names(GOTables)){
  GOTables[[nm]]$group = nm
}
" Concatenate GO tables of different groups (for GO comparison later) " %>% message()
tb_cated = bind_rows(GOTables)
tb_cated$group = factor(
  tb_cated$group, levels = allStages[allStages %in% tb_cated$group])
tb_cated$group %>% table()
str(tb_cated)


##############################[ single table ]#################################
"handling a single table, preparing for barplot" %>% message()
ont_choices = c("BP", "CC", "MF")
stg = names(GOTables)[2]
for (stg in names(GOTables)){
  n_show_bar = 25
# ont = ont_choices[1]
  for (ont in ont_choices){
    tb_plt = GOTables[[stg]] %>% 
      filter(ONTOLOGY == ont) %>%
      group_by(ONTOLOGY) %>% 
      top_n(n_show_bar, - pvalue)
    str(tb_plt)
    
    "================ plot GO terms (one group) ===============" %>% message()
    tag_bar = sprintf("%s_%s", ont, stg)
    tt = sprintf("Top enriched GO terms (%s)", stg)
    fn_bar = sprintf("%s/bar_%s.pdf", figdir, tag_bar)
    
    source("RFunxGO.R")
    plot_GO_bar(tb_plt, title = tt, fn = fn_bar, wscale = 1.2)
  }
}

##############################[ multiple groups ]#################################
" FIRST taking the names of the top GO terms, then subsetting the table! "
n_show_compare = 10
groups_compare = tb_cated$group %>% unique()
groups_compare = c("4cell", "8cell", "32cell", "256cell", "B", "G3")[-c(1, 4)]
cut_pv = 0.02
ont_choices = c("BP", "CC", "MF")
ont = ont_choices[1]
source("RFunxGO.R")

for (ont in ont_choices){
  
goid_use = tb_cated %>% 
  filter(ONTOLOGY == ont, group %in% groups_compare) %>%
  group_by(group) %>% 
  top_n(n_show_compare, wt = - pvalue)
goid_use = goid_use$ID

" subsetting table, filterring by p-values "
tb_plt = filter(tb_cated, 
                ID %in% goid_use, 
                group %in% groups_compare, 
                pvalue < cut_pv)
tb_plt$group %>% table
str(tb_plt)

# "=== removing terms with tooo looong description ===" %>% message()
# terms_rmv = c("RNA splicing, via transesterification reactions with bulged adenosine as nucleophile")
# tb_plt = tb_plt %>% filter(! Description %in% terms_rmv)

"============= plot GO terms (multiple groups) ============" %>% message()
sizeby = c("GeneRatio", "Count")[1]
colorby = c("pvalue", "p.adjust")[1]

tt_c = sprintf("Top enriched GO terms")
tag_compare = sprintf("%s-pv%.3g_sz%s-%s", ont, cut_pv, sizeby, paste0(groups_compare, collapse = "&"))
fn_compare = sprintf("%s/compare_%s.pdf", figdir, tag_compare)

plot_GO_compare(
  tb_plt, x = c("group", "Cluster", "cluster")[1],
  y = "Description",
  title = tt_c, sizeby = sizeby, colorby = colorby,
  fn = fn_compare,
  xrotation = 45,
  wscale = 0.8, # 0.85 for most
  hscale = 1
)
}









