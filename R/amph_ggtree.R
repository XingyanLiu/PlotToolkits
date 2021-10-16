# pkgs = c("nlme")
# install.packages(pkgs[1])
# BiocManager::install("ggtree")
library(ape)
library(tidytree)
library(ggplot2)
library(ggtree)
library(RColorBrewer)


# ============== setting direcories ============
DATADIR = "E:/others/003_scr/test_data/devTree"
# DATADIR = "E:/Users/xyliu/data003/amph/tree"
figdir = file.path(DATADIR, 'figs')
if(!dir.exists(figdir)){dir.create(figdir, recursive = T)}

# ============== load tree ============
df_struct = read.csv(sprintf("%s/tree_struct.csv", DATADIR), 
                     na.strings = '', 
                     as.is = c("parent", "node", "label")) # 
df_expr = read.csv(sprintf("%s/avg_expr_all.csv", DATADIR), row.names = 1)
df_exprprop = read.csv(sprintf("%s/expr_prop_all.csv", DATADIR), row.names = 1)

tb_tree = as_tibble(df_struct[, c("parent", "node", "label")]) %>% filter(!is.na(parent))
tb_extdata = as_tibble(df_struct[, c("label", "stage", "stage_int")])
tb_extdata$stage = factor(tb_extdata$stage, levels = paste0("E", 1:15), ordered = T)
# tb_extdata$length = 1


tree0 = as.phylo(tb_tree)
tree0$edge.length = 1
ggtree(tree0) + geom_tiplab() + geom_nodelab()

# prepare for stage-tree
treedt = full_join(x = as_tibble(tree0),
                   y = tb_extdata, # y = tb_data_expr,
                   by='label') %>% as.treedata()
# # prepare for expression display
# treedt_expr = full_join(x = as_tibble(tree0),
#                    y = tb_expr,
#                    by='label') %>% as.treedata()
# str(treedt)
# treedt@data -> tmp
# =======================================================
# ===================== plot tree =======================
# =======================================================
# treedt = my_tree2
# colored by staage
# levels(p$data$stage)
colors_stg = brewer.pal(11, 'Spectral')
if(FALSE){colors_stg = rev(colors_stg)}
colors_stg1 = colorRampPalette(colors_stg)(15)
# === (0) basic layout ===
p = ggtree(treedt, size=1.5, color='lightgrey', ladderize = F)#, layout = 'slanted')

# ==== (1) node points ===
p0 = p + 
  geom_point(aes(color=stage), size = 4, alpha=1) +
  scale_color_manual(values = colors_stg1) +
  # scale_color_gradientn(colors = colors_stg) +
  geom_tiplab(aes(label=label), color='darkgrey', 
              align=TRUE, hjust = -0.4, size=3) + 
  # geom_nodelab(aes(x=x-0.3, label=label), color='darkgrey', size=3) + 
  theme(legend.position="right") +
  xlim(0, 15) # leave more space for tip labels
ggsave(p0, file=sprintf("%s/devTree.pdf", figdir), width=10, height=8)


# ==== (2) filled node labels ====
p1 = p + 
  geom_label(aes(label=label, fill = stage_int), alpha = 0.9, color='black', size=2.5) +
  # scale_color_manual(values = colors_stg1) + 
  scale_fill_gradientn("stage", colors = colors_stg) +
  # geom_tiplab(aes(label=label), color='darkgrey', align=TRUE, hjust = -0.5, size=3) + 
  theme_tree(fgcolor='grey', legend.position="right") +
  xlim(0, 14.5) # leave more space for tip labels
ggsave(p1, file=sprintf("%s/devTree_labeled.pdf", figdir), width=10, height=8)


# ==== (3) circular tree ====
p2 = ggtree(treedt, size=1.5, color='lightgrey', layout = 'circular', ladderize = T) + 
  geom_point(aes(color=stage), size = 4, alpha=1) +
  scale_color_manual(values = colors_stg1) +
  # scale_color_gradientn(colors = colors_stg) +
  geom_tiplab2(aes(label=label), color='darkgrey',
               align=TRUE, hjust = -0.2, size=4) + 
  # geom_nodelab(aes(x=x-0.3, label=label), color='darkgrey', size=3) + 
  theme_tree(legend.position="right") +
  xlim(0, 15.5) # leave more space for tip labels
ggsave(p2, file=sprintf("%s/devTree_cricular.pdf", figdir), width=10, height=8)


# ========== (4) up-down ===========
p00 = p + 
  geom_point(aes(color=stage), size = 4.5, alpha=1) +
  scale_color_manual(values = colors_stg1) +
  geom_tiplab(aes(label=label), 
              color='darkgrey', 
              angle=90,
              hjust=1, offset=-0.2,
              align=TRUE, size=3) + 
  # theme_tree(fgcolor='grey', legend.position="right") +
  scale_x_reverse() + coord_flip() +
  theme(legend.position="right") 
ggsave(p00, file=sprintf("%s/devTree_updown.pdf", figdir), width=10, height=8)



p + 
  geom_point(aes(fill=stage), color='grey', size = 4.5, alpha=1, shape=21) +
  scale_fill_manual(values = colors_stg1) +
  # scale_color_gradientn(colors = colors_stg) +
  geom_tiplab(aes(label=label), color='darkgrey', 
              angle=90,
              hjust=1, offset=-0.2,
              align=TRUE, size=3) + 
  scale_x_reverse() + coord_flip() +
  theme(legend.position="right") 
ggsave(file=sprintf("%s/devTree_updown1.pdf", figdir), width=10, height=8)

# ======================================================
# ============== gene expression on tree ===============
# ======================================================
chord_gene = "bf_00020985"
candi_genes = c(
  'bf_00003294',
  'bf_00019235',
  'bf_00004470',
  # 'bf_00001725',
  'bf_00007978',
  'bf_00008183'
)

# ===== prepare for expression display ======
# settings
colors_expr = brewer.pal(7, "Reds")
colors_expr[1] = "#F0F0F0"

# gene_name = "gene1"
# tb_expr = tibble(label = df_struct$label, gene1 = rnorm(dim(df_struct)[1])) # for testing
gene_name = chord_gene
gene_name = candi_genes[6]

if(gene_name %in% rownames(df_expr)){
tb_expr = as_tibble(t(df_expr[gene_name, ]))
tb_expr$proportion = t(df_exprprop[gene_name, ])[, 1]
tb_expr$label = colnames(df_expr)
treedt_expr = full_join(x = as_tibble(tree0),
                        y = tb_expr,
                        by='label') %>% as.treedata()


# === (0) basic layout ===
mapping_expr = parse(text = sprintf("aes(color=%s, size=proportion)", gene_name))
p_expr = ggtree(treedt_expr, size=1, color='lightgrey',#mapping = eval(mapping_expr), 
                ladderize = FALSE)#, continuous = TRUE
}else{
  message(sprintf("No expressions of %s!, STOP!", gene_name))
  }
# === node points ===
# mapping_expr_nd = parse(text = sprintf("aes(label=label, x = branch + 0.2, color=%s)", gene_name))
{
p_expr1 <-
p_expr + 
  geom_point(eval(mapping_expr), alpha=1) +
  scale_color_gradientn(colors = colors_expr) +
  geom_tiplab(mapping = aes(label=label), color='grey', 
              align = TRUE,
              size=3, hjust = -0.3) +
  # geom_nodelab(aes(label=label, x = branch + 0.1, y=y+0.2), color='grey', size=3) +
  # geom_nodelab(eval(mapping_expr_nd), size=3) +
  xlim(0, 15) +
  # theme_tree2(fgcolor='lightgrey', legend.position="right")
  theme(legend.position="right") 
ggsave(p_expr1, file=sprintf("%s/exprTree_%s.pdf", figdir, gene_name), width=10, height=8)
}
##=======  label the nodes =======
{
mapping_expr_lb = parse(text = sprintf("aes(label=label,fill=%s)", gene_name))
p_expr2 <-
p_expr + 
  scale_color_gradientn(colors = colors_expr) +
  geom_label(mapping = eval(mapping_expr_lb), color='black', size=3) +
  scale_fill_gradientn(colors = colors_expr) +
  # scale_x_reverse() + coord_flip() +
  # theme_tree2(fgcolor='lightgrey', legend.position="right")
  theme(legend.position="right")
# scale_fill_brewer(palette = 'Spectral') for discrete values
ggsave(p_expr2, file=sprintf("%s/exprTree_labeled_%s.pdf", figdir, gene_name), width=10, height=8)
}

# === node points (Up-down) ===
# mapping_expr_nd = parse(text = sprintf("aes(label=label, x = branch + 0.2, color=%s)", gene_name))
{
  p_expr11 <-
    p_expr + 
    geom_point(eval(mapping_expr), alpha=1) +
    scale_color_gradientn(colors = colors_expr) +
    geom_tiplab(aes(label=label), 
                color='grey', 
                angle=90,
                hjust=1, offset=-0.2,
                align=TRUE, size=3) + 
    scale_x_reverse() + coord_flip() +
    # geom_nodelab(aes(label=label, x = branch + 0.1, y=y+0.2), color='grey', size=3) +
    # geom_nodelab(eval(mapping_expr_nd), size=3) +
    theme(legend.position="right") 
  ggsave(p_expr11, file=sprintf("%s/exprTree_updown_%s.pdf", figdir, gene_name), width=10, height=8)
}











##======= label on the branches (too ugly...) =======
p_expr + 
  geom_point(eval(mapping_expr), alpha=0.99) +
  scale_color_gradientn(colors = colors_expr) +
  geom_tiplab(mapping = aes(label=label, x = branch, y=y+0.2), color='grey', size=3) +
  geom_nodelab(aes(label=label, x = branch + 0.1, y=y+0.2), color='grey', size=3) +
  # geom_nodelab(eval(mapping_expr_nd), size=3) +
  xlim(0, 15) +
  # theme_tree2(fgcolor='lightgrey', legend.position="right")
  theme(legend.position="right")







#==================draft==================

display.brewer.pal(7, 'Reds')
display.brewer.pal(9, 'Greys')

brewer.pal(9, 'Greys')[2]




tr00 = rtree(5)
str(tr00)
tr0 = as_tibble(tr00)
my_tbl0 = tibble(parent = paste0('E', tr0$parent),
         node = paste0('E', tr0$node),
         edge.length = tr0$branch.length, #ifelse(tr0$parent == tr0$node, NA, 1),
         label = paste0('E', tr0$node), #paste0('t', tr0$node),
         weight = tr0$branch.length
         )


my_tbl = my_tbl0 %>% filter(!is.na(weight)) # important!

my_tree0 = as.phylo(my_tbl)
str(my_tree0)
my_tree0$edge.length = 1#my_tbl$branch.length
my_tree0

ggtree(my_tree0) + geom_tiplab() + geom_nodelab()

as_tibble(my_tree0) %>% as.phylo()





















