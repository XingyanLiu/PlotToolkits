library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(RColorBrewer)
library(cowplot)
library(Matrix)
library(matrixStats)
library(dplyr)


MAINDIR = "E:/lxy_pro/003"   #"E:/others/003_scr"
setwd(MAINDIR)


# ===================================================================
#        setting direcories 
# ===================================================================

refdir = "resources/ref_markers"
DATADIR = sprintf("test_data/devTree/%s", "20201014") # "20200707"
# DATADIR = "E:/Users/xyliu/data003/amph/tree"

figdir = sprintf("%s/figs%s", DATADIR, '-1120')
if(!dir.exists(figdir)){dir.create(figdir, recursive = T)}


# ===================================================================
#         loading 1) tree-structure and 2) group-expressions
# ===================================================================

message("loading edges of developmental tree...")
tag = ""         #"sepPCA_"
tailtag = c("", "-labeled", "-rmvE15")[2]     #"-rmvE15"
df_struct0 = read.csv(sprintf("%s/%stree_struct%s.csv", DATADIR, tag, tailtag), 
                     na.strings = '', 
                     as.is = T) # c("parent", "node", "label"))
message("loading average expressions and expresson proportions...")
df_expr = read.csv(sprintf("%s/%savg_expr_all.csv", DATADIR, tag), row.names = 1)
df_exprprop = read.csv(sprintf("%s/%sexpr_prop_all.csv", DATADIR, tag), row.names = 1)
Genes = rownames(df_expr)

str(df_struct0)
# ============== prepare for tree plot ===============
source("RFunxTreePlot.R")

#####[1] re-order the nodes for tree-branch order
rownames(df_struct0) = df_struct0$node
nodes_early = c(paste0("E6_", c(2, 0, 1, 3)), 
              paste0("E7_", c(1, 5, 2, 3, 4, 0, 6)),
              paste0("E8_", c(2, 1, 3, 5, 4, 0, 6)),
              paste0("E9_", c(4, 1, 8, 5, 7, 3, 6, 2, 0, 9)))
nodes_late = subset(df_struct0, ! stage %in% c("E6", "E7", "E8", "E9"))$node
nodes_ordered = c(nodes_early, nodes_late)
df_struct = df_struct0[nodes_ordered, ]

#####[2] rename stages
levels(df_struct$stage) <-
  levels(df_struct$stage) %>% sapply(function(x){StageNameMap[[x]]})

str(df_struct)

##### basic structure
tree0 = makePhyloFromDf(df_struct)
str(tree0)
ggtree(tree0) + geom_tiplab() + geom_nodelab()


# =====================[ stages ]=========================
source("RFunxTreePlot.R")

treedt = makeTreeDataStage(df_struct, stg_levels = StageNames)
treedt@data
ggtree(treedt) + geom_nodelab()

##### formal code ======================
p01 = plotLabeledStageTree(treedt, anno_list = LineageAnno0, 
                           xmax = 13, #colors_stg = 
                           filename = sprintf("%s/devTree_labeled.pdf", figdir))
p020 = plotStageTreeLR(treedt, anno_list = LineageAnno0, 
                       size_tiplab=2, width = 7, height = 4.5,
                       xmax = 14, offset.text = 0.16,
                       filename = sprintf("%s/devTree_LR.pdf", figdir))
p022 = plotStageTreeLR1(treedt, anno_list = LineageAnno0,
                        size_tiplab=2, width = 6.5, height = 4.5,
                        xmax = 14,
                        offset.text = 0.16,
                        filename = sprintf("%s/devTree_LR1.pdf", figdir))
p023 = plotStageTreeLR1(treedt, anno_list = LineageAnno0,
                        size_tree = 0.6, color_tree = 'grey30',
                        size_tiplab=2, width = 6.5, height = 4.5,
                        xmax = 14,
                        offset.text = 0.16,
                        filename = sprintf("%s/devTree_LR1-light.pdf", figdir))

# p021 = plotStageTree(treedt, size_tiplab=2,# height = 5,
#                      filename = sprintf("%s/devTree.pdf", figdir))





# =============================================
# temp: selected genes of interest
# =============================================

figdir = file.path(DATADIR, 'figs')
dir.create(figdir)

Annos = readRDS("resources/Annos.rds")
str(Annos)
rownames(Annos)[1:5]

# figdir3 = file.path(DATADIR, 'figs-stage')
# dir.create(figdir3)
# # gids = c("bf_00014181", "bf_00021102", "bf_00001993", "bf_00000872")
# gids = c('bf_00000872', 'bf_00001993', 'bf_00023369', 'bf_00003252',
#           'bf_00000722', 'bf_00014181', 'bf_00006882', 'bf_00010279',
#           'bf_00000536', 'bf_00013295', 'bf_00021102', 'bf_00005799',
#           'bf_00000272', 'bf_00009749', 'bf_00009656', 'bf_00006993',
#           'bf_00026470', 'bf_00004621', 'bf_00015478')

##################[ Case-1: for only one list of genes ]##################
# load selected gene ids (or names)
.dir = "E:/lxy_pro/003/amph/GOresults/1015-sublineage(MAST)"
.fp = sprintf("%s/res-cutFC0.2/__selected_gene_idnames-E7_5.csv", .dir)
gids = read.csv(.fp, as.is=T)
gids = gids$ID %>% unique()
tagname = "selected_genes-E7_5"

fn_save = sprintf("%s/exprTrees_%s.pdf", figdir, tagname)
plotTreesBothLR(gids, tree0, df_expr, df_exprprop,
                anno_list_lineage = LineageAnno,
                annos=Annos, cln = "tt",
                fname=NULL) -> plts
{
  pdf(fn_save, width = 9, height = 6)
  invisible(lapply(plts, print))
  dev.off()
}


# ###
# lin_gene_ids = filter(LineageDEGs_use, lineage == "E7_5.E12_1")$ID
# gids = intersect(gids, lin_gene_ids)

##################[ Case-2: for multiple lists of genes ]##################
listGeneIDs = list(
  # tmp = "bf_00001993",
  CHRD = "bf_00020985",
  LRIG3 = "bf_00008183",
  tail_bud_xen = c("bf_00025359", "bf_00004503", "bf_00012984", 
                   "bf_00013232", "bf_00013231", "bf_00011511"),
  CDX = c("bf_00012984", "bf_00011511"),
  PGC = c("bf_00004829", "bf_00016526", "bf_00018121", "bf_00005532"),
  stem_cells = c("bf_00006715", "bf_00013097", "bf_00019574"),
  neural_ectoderm = c("bf_00009455", "bf_00006429"),# E7_4
  notochord = c("bf_00001923", "bf_00006621", "bf_00020985"),
  mesoderm = c("bf_00014595", "bf_00023044", "bf_00004969", "bf_00024930"),
  endoderm = c("bf_00000272", "bf_00024930", "bf_00005799", "bf_00000536"),
  epithelial_ectoderm = c("bf_00003323", "bf_00000489")
  # sox2 = c("bf_00017975", "bf_00001714", "bf_00001715"),
  # Foxp1_4 = "bf_00008448",
  
) 
for (tlname in names(listGeneIDs)[1: 2]){ #"Foxp1_4"
  
  gids = listGeneIDs[[tlname]]
  fname = sprintf("%s/exprTrees_%s.pdf", figdir, tlname)
  plotTreesBothLR(gids, tree0, df_expr, df_exprprop,
                  anno_list_lineage = LineageAnno,
                  annos=Annos, cln = "tt",
                  fname=fname) -> plts
  
}

























