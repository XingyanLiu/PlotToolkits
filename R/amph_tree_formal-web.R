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
source("RFunxTreePlot.R")


# ===================================================================
#        setting direcories 
# ===================================================================

date_tag = c("20210115", "20201014", "20200707")[1]
DATADIR = sprintf("test_data/devTree/%s", "20201014") # "20200707"



# ===================================================================
#         loading 1) tree-structure and 2) group-expressions
# ===================================================================
Annos = readRDS("resources/Annos.rds")
str(Annos)

message("loading edges of developmental tree...")
tag = ""         #"sepPCA_"
tailtag = c("", "-formal", "-labeled", "-rmvE15")[2]     #"-rmvE15"
df_struct0 = read.csv(sprintf("%s/%stree_struct%s.csv", DATADIR, tag, tailtag), 
                      na.strings = '', 
                      as.is = T) # c("parent", "node", "label"))
message("loading average expressions and expresson proportions...")
df_expr = read.csv(sprintf("%s/%savg_expr_all.csv", DATADIR, tag), row.names = 1)
df_exprprop = read.csv(sprintf("%s/%sexpr_prop_all.csv", DATADIR, tag), row.names = 1)
Genes = rownames(df_expr)

str(df_struct0)

##########################################################
"changing the column names of the expression dataframes"
library(stringr)
colnms = colnames(df_expr) %>% sapply(function(nm){
  .tmp = str_split_fixed(nm, "_", 2)
  nm = paste(StageNameMap[[.tmp[, 1]]], .tmp[, 2], sep = "_")
})

colnames(df_expr) = colnms
colnames(df_exprprop) = colnms


# ============== prepare for tree plot ===============
source("RFunxTreePlot.R")

#####[1] re-order the nodes for tree-branch order
rownames(df_struct0) = df_struct0$node
nodes_early = c(paste0("B_", c(2, 0, 1, 3)), 
                paste0("G3_", c(1, 5, 2, 3, 4, 0, 6)),
                paste0("G4_", c(2, 1, 3, 5, 4, 0, 6)),
                paste0("G5_", c(4, 1, 8, 5, 7, 3, 6, 2, 0, 9)))
nodes_late = subset(df_struct0, ! stage %in% c("B", "G3", "G4", "G5"))$node
nodes_ordered = c(nodes_early, nodes_late)
df_struct = df_struct0[nodes_ordered, ]

#####[2] rename stages
levels(df_struct$stage) <-
  levels(df_struct$stage) %>% sapply(function(x){StageNameMap[[x]]})

str(df_struct)

##### basic structure
tree0 = makePhyloFromDf(df_struct)

########################[ jpeg ]#######################
"jpeg"
source("RFunxTreePlot.R")

figdir_ALL = "F:/amph_exprTree_jpeg"
all_gids = rownames(df_expr)
length(all_gids)
.gid = all_gids[1]

for (.gid in all_gids){
  .gnm = gsub('/', '_', x=Annos[.gid, 'gene_short_name'])
  .fn_save = sprintf("%s/exprTree-%s(%s).jpeg", figdir_ALL, .gid, .gnm)
  
  if (file.exists(.fn_save)){
    sprintf("skipped for %s(%s) has already been plotted.", .gid, .gnm) %>% message()
    next
  }
  plotExprTreeLR_full(tree0, .gid, df_expr, df_exprprop,
                      anno_list = LineageAnnoDetailed,
                      annos=Annos, cln = "tt",
                      size_tree = 1.,
                      size_nodelb = 5,
                      size_node_range=c(3, 12.5),
                      size_xticks = 15,
                      xticks=c("", "G3", "G4", "G5", "G6", "N0", "N1", "N3", "L0"),
                      height=600,
                      filename = .fn_save,
                      fontsize=5 # lineage annotations
                      ) -> p
  
}








