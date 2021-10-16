
library(ggplot2)
library(dplyr)


df0 = subset(obj@meta.data, stage_name != "B")

# labels_stage = df0$stage_name
# labels_lineage = df0$lineage
# cross_tb = table(labels_lineage, labels_stage)
# cross_tb %>% as.matrix() #%>% as_tibble()

df = tibble(
  stage = factor(df0$stage_name, levels = StageOrd),
  lineage = factor(df0$lineage)#, levels=LineageOrdAll[-c(1:4)])
)
df

ylb = c("cell percentage", "cell counts")[1]
p = ggplot(df) + geom_bar(aes(x = stage, fill=lineage), position = "fill")
p + scale_fill_manual(values = LineageColors) +
  ylab(ylb) +
  theme_classic()
ggsave(sprintf("%s/%s_stage_lineage.pdf", figdir, ylb),
       width = 5, height = 3)

figdir

