library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(RColorBrewer)
library(cowplot)
library(Matrix)
library(matrixStats)
library(dplyr)
library(stringr)



MAINDIR = "E:/lxy_pro/003"   #"E:/others/003_scr"
setwd(MAINDIR)
source("RFunxTreePlot.R")


# ===================================================================
#        setting direcories 
# ===================================================================

DATADIR = sprintf("test_data/devTree/%s", "20201014") # "20200707"
# DATADIR = "E:/Users/xyliu/data003/amph/tree"
date_tag = c("20210115", "20201014", "20200707")[1]
DATADIR = sprintf("test_data/devTree/%s", date_tag) # "20200707"

figdir = sprintf("%s/figs%s", DATADIR, '')
if(!dir.exists(figdir)){dir.create(figdir, recursive = T)}


# ===================================================================
#         loading 1) tree-structure and 2) group-expressions
# ===================================================================

message("loading edges of developmental tree...")
tailtag = c("", "-formal", "-labeled", "-rmvE15")[2]     #"-rmvE15"
df_struct0 = read.csv(sprintf("%s/tree_struct%s.csv", DATADIR, tailtag), 
                     na.strings = '', 
                     as.is = T) # c("parent", "node", "label"))
rownames(df_struct0) = df_struct0$node
str(df_struct0)

message("loading average expressions and expresson proportions...")
df_expr = read.csv(sprintf("%s/avg_expr_all.csv", DATADIR), row.names = 1)
df_exprprop = read.csv(sprintf("%s/expr_prop_all.csv", DATADIR), row.names = 1)
Genes = rownames(df_expr)

# ============== prepare for tree plot ===============
source("RFunxTreePlot.R")

#####[1] re-order the nodes for tree-branch order
nodes_early = c(paste0("B_", c(1, 0, 2)), 
              paste0("G3_", c(1, 5, 2, 3, 4, 0, 6)),
              paste0("G4_", c(2, 1, 3, 5, 4, 0, 6)),
              paste0("G5_", c(4, 1, 8, 5, 7, 3, 6, 2, 0, 9)))
nodes_late = subset(df_struct0, ! stage %in% c("B", "G3", "G4", "G5"))$node
nodes_ordered = c(nodes_early, nodes_late)
df_struct = df_struct0[nodes_ordered, ]
head(df_struct0)

#####[2] rename stages (unnecessary)
head(df_struct)

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
root_lb = "32cell"
p01 = plotLabeledStageTree(treedt, anno_list = LineageAnno, 
                           xmax = 13, #colors_stg = 
                           filename = sprintf("%s/devTree_labeled.pdf", figdir))
p020 = plotStageTreeLR(treedt, anno_list = LineageAnno, 
                       size_tiplab=2, width = 7, height = 4.5,
                       xmax = 14, offset.text = 0.16,
                       filename = sprintf("%s/devTree_LR.pdf", figdir))
p022 = plotStageTreeLR1(treedt, anno_list = LineageAnno,
                        root_lb = root_lb,
                        size_tiplab=2, width = 6.5, height = 4.5,
                        xmax = 14,
                        offset.text = 0.16,
                        filename = sprintf("%s/devTree_LR1.pdf", figdir))
p023 = plotStageTreeLR1(treedt, anno_list = LineageAnno,
                        root_lb = root_lb,
                        size_tree = 0.6, color_tree = 'grey30',
                        size_tiplab=2, width = 6.5, height = 4.5,
                        xmax = 14,
                        # filename = sprintf("%s/devTree_LR1-light.pdf", figdir),
                        offset.text = 0.16
                        )

# p021 = plotStageTree(treedt, size_tiplab=2,# height = 5,
#                      filename = sprintf("%s/devTree.pdf", figdir))
################################################################################
source("RFunxTreePlot.R")
xmax = 9
stages = c("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0")
# xticks = c(stages, rep("", xmax - length(stages)))

p_fig1 = plotStageTreeLR1(treedt, anno_list = LineageAnnoDetailed,
                        root_lb = root_lb,
                        size_node = 5.5,
                        size_tree = 0.6, color_tree = 'black',#'grey30',
                        alpha_point = 1,
                        size_tiplab=2, width = 6.5, height = 4.5,
                        xmax = 14,
                        # filename = sprintf("%s/devTree_LR1-light.pdf", figdir),
                        offset.text = 0.16,
                        default_colr = 'grey20',
                        fontsize=3.2,
)
p_fig1
p_fig11 = p_fig1 + geom_nodelab(
  aes(x = x, 
      label=str_split_fixed(label, '_', 2)[, 2]), 
  color='black', #'grey5',
  size=3) +
  geom_tiplab(
    aes(x = x - 0.12, 
        label=str_split_fixed(label, '_', 2)[, 2]), 
    color='black', #'grey5',
    size=3) +
  scale_x_discrete(limits=c(1:xmax), labels = stages) +
  coord_cartesian(clip = 'off') +
  theme_tree2(plot.margin=margin(6, 150, 6, 6)) +
  theme(legend.position="left")
p_fig11
hlightClade(p_fig11, anno_list=LineageAnno, 
            anno_color_list=LineageColor, #list(), 
            # barsize=1, 
            extend=0.4,
            alpha = 0.1,
            # align = T, offset = 0.4,
            default_colr = 'grey90'
) -> p_fig111
ggsave(filename = sprintf("%s/devTree_LR1-light-labeled0.pdf", figdir),
       width = 6.5, height = 4.5,
       plot = p_fig11)
ggsave(filename = sprintf("%s/devTree_LR1-light-labeled.pdf", figdir),
       width = 6.5, height = 4.5,
       plot = p_fig111)

# =============================================
# temp: selected genes of interest
# =============================================

# figdir = file.path(DATADIR, 'figs')
# dir.create(figdir)

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
"load selected gene ids (or names)" %>% message()
.dir = "E:/lxy_pro/003/amph/GOresults/1015-sublineage(MAST)"
.fn = "res-cutFC0.2/__selected_gene_idnames-E7_5.csv"
tagname = "selected_genes--E7_5"
.figdir = figdir

.dir = "E:/lxy_pro/003/amph/genes"
.fn = "space_gene_ids.txt"
tagname = "space"
.figdir = figdir

fn_save = sprintf("%s/exprTrees_%s-.pdf", .figdir, tagname)

.fp = sprintf("%s/%s", .dir, .fn)
print(.fp)
gids = read.csv(.fp, as.is=T)
gids = gids$ID %>% unique()

plotTreesBothLR(gids, 
                tree0, df_expr, df_exprprop,
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



##########################################

figdir_ALL = "F:/amph_exprTree_ALL"
all_gids = rownames(df_expr)
length(all_gids)
.gid = all_gids[1]

for (.gid in all_gids){
  .gnm = gsub('/', '_', x=Annos[.gid, 'gene_short_name'])
  .fn_save = sprintf("%s/exprTree-%s(%s).pdf", figdir_ALL, .gid, .gnm)
    if (file.exists(.fn_save)){
      sprintf("skipped for %s(%s) has already been plotted.", .gid, .gnm) %>% message()
      next
    }
  # plotTreesBothLR(.gid, 
  #                 tree0, df_expr, df_exprprop,
  #                 anno_list_lineage = LineageAnno,
  #                 annos=Annos, cln = "tt",
  #                 fname=.fn_save
  #                 )
  plotExprTreeLR_full(tree0, .gid, df_expr, df_exprprop,
                      anno_list = LineageAnnoDetailed,
                      annos=Annos, cln = "tt",
                      size_tree = 1.,
                      size_nodelb = 3,
                      size_node_range=c(1, 7.5),
                      size_xticks = 12,
                      nodelb_offset = - 0.,
                      tiplb_offset = - 0.1,
                      xticks=c("B", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0"),
                      filename = .fn_save,
                      fontsize=4 # lineage annotations
  ) -> p
  # p

}

.gmax = apply(df_expr, 1, max)
sum(.gmax == 0)

hist(log1p(.gmax))

##########################################







