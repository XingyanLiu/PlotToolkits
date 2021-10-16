library("dplyr")
library("stringr")
library("clusterProfiler")
library("forcats")



map2GOsymbols = function(genes, db = "org.Hs.eg.db",
                         fromType="SYMBOL", toType="ENTREZID", 
                         icol = 1){
  mapping_df = bitr(genes, 
                    fromType=fromType, toType=toType, 
                    OrgDb=db)
  
  if(icol == 'all'){
    return(mapping_df)
  }else{
    mapping_df[, icol]
  }
}


toGOsymbolsZebr = function(gene_symbols,
                           upper = c("nd5","cox3", "cox1"),
                           db = "org.Dr.eg.db", #ignored
                           icol = 1){
  gene_symbols = gene_symbols %>% tolower()
  gene_symbols %>% sapply(function(x){
    ifelse(x %in% upper, toupper(x), x)
  })
  
  map2GOsymbols(gene_symbols, db=db, icol=icol)
}



toGOsymbolsXen = function(gene_symbols, 
                            db = "org.Xl.eg.db",
                            icol = 1){
  
  gene_symbols = tolower(gene_symbols)
  gene_symbols1 = paste0(gene_symbols, ".L")
  gene_symbols2 = paste0(gene_symbols, ".S")
  
  syms_candi = c(gene_symbols, gene_symbols1, gene_symbols2)
  map2GOsymbols(syms_candi, db = db, icol = icol)
  
  
}



toGOsymbolsXen2 = function(gene_symbols, 
                            db = "org.Xl.eg.db",
                            icol = 1,
                            exts = c(".L", ".S")){
  
  gene_symbols = tolower(gene_symbols)
  gene_symbols1 = paste0(gene_symbols, exts[1])
  
  syms_candi0 = c(gene_symbols, gene_symbols1)
  sym2id0 = bitr(syms_candi0, 
                 fromType="SYMBOL", toType="ENTREZID", 
                 OrgDb=db)
  
  unmapped1 = setdiff(gene_symbols1, sym2id0[, 1])
  if (length(unmapped1) >= 1){
    syms_candi1 = str_replace(unmapped1,
                              pattern = exts[1], replacement = exts[2])
    # print(syms_candi1)
    sym2id1 = bitr(syms_candi1, 
                   fromType="SYMBOL", toType="ENTREZID", 
                   OrgDb=db)
    sym2id = bind_rows(sym2id0, sym2id1)
  }else{
    sym2id = sym2id0
  }
  
  if(icol == 'all'){
    return(sym2id)
  }else{
    sym2id[, icol]
  }
  
}

################################################################################


quickGO = function(
  gene_syms,
  OrgDb,# = org.Hs.eg.db,
  tag = "",
  tag_title = "",
  ont_choices = c("ALL", "BP", "MF", "CC"),
  keyType = c("SYMBOL", "ENSEMBL"),
  do_simplify = TRUE,
  do_save_plots = TRUE,
  save_rds = TRUE,
  save_csv = TRUE,
  n_show_bar = 30,
  n_show_dot = 30,
  n_show_net = 6,
  resdir = "tempGOresults",
  resdir_rds = NULL,
  figdir = NULL,
  ...){
  
  if (save_rds && is.null(resdir_rds)){
    resdir_rds = file.path(resdir, "res.rds")
    dir.create(resdir_rds, recursive = T)
  }
  if (do_save_plots && is.null(figdir)){
    figdir = file.path(resdir, "figs")
    dir.create(figdir, recursive = T)
  }
  
  results = list()
  
  for (ont in ont_choices){
    print(ont)
    ftag = sprintf("%s-%s", ont, tag)
    ego = enrichGO(gene_syms, 
                   OrgDb = OrgDb, 
                   keyType = keyType, #c("SYMBOL", "ENSEMBL"),
                   ont = ont,
                   qvalueCutoff = 0.2,
                   pvalueCutoff = 0.2,
                   # minGSSize = 10,
                   maxGSSize = 500,
                   pool = F,
                   readable = F,
                   ...)
    
    print("analysis completed~")
    ego_plt = ego
    
    if (ont != "ALL" && do_simplify){
      message("simplifying...(removing redundant terms)")
      sego = simplify(ego)
      ego_plt = sego
      ### saving results
      if (save_rds){
        saveRDS(sego, sprintf("%s/sego_%s.rds", resdir_rds, ftag))
      }
      if(save_csv){
        write.csv(sego@result, 
                  sprintf("%s/enrich-%s.csv", resdir, ftag), 
                  row.names = F)
      }    
    }
    
    # ====[ saving results (un-simplified) ]====
    ### RDS
    if (save_rds){
      saveRDS(ego, sprintf("%s/ego_%s.rds", resdir_rds, ftag))
    }
    ### CSV
    if (save_csv){
      fn_restb = sprintf("%s/enrich-%s.csv", resdir, ftag)
      write.csv(ego@result, 
                fn_restb, 
                row.names = F)
      message(sprintf("enrich results table saved into\n\t%s", fn_restb))
    }    
    # =================[ plotting ]=================
    restb = ego@result
    print(dim(restb))
    print(summary(restb$p.adjust))
    print(summary(restb$pvalue))
    
    if (do_save_plots && min(restb$p.adjust) < ego@pvalueCutoff){
      
      tt = sprintf("Top %d enriched GO terms %s", n_show_bar, tag_title)
      p_bar = barplot(ego_plt, showCategory = n_show_bar) + ggtitle(tt)
      tt = sprintf("Top %d enriched GO terms %s", n_show_dot, tag_title)
      p_dot = dotplot(ego_plt, showCategory = n_show_dot) + ggtitle(tt)
      
      fn_go_plts = sprintf("%s/enrichVis%s.pdf", figdir, ftag)
      pdf(file = fn_go_plts, width = 11, height = 6)
      invisible(lapply(list(p_bar, p_dot), print))
      dev.off()
      # n_show_net = 6
      tt_net = sprintf("Top %d enriched GO terms %s", n_show_net, tag_title)
      p_net <- cnetplot(ego_plt, showCategory = n_show_net) + ggtitle(tt_net)
      
      ggsave(filename = sprintf("%s/enrichNet%s.pdf", figdir, ftag),
           plot = p_net,
           width = 6.5, height = 6)
      
      print(figdir)
    }

  }
  return(results)
  }



quick_compareGO = function(
  genelist,
  OrgDb, #  = org.Hs.eg.db,
  fun="enrichGO",
  tag = "",
  tag_title = "",
  ont_choices = c("BP", "MF", "CC"),
  keyType = c("SYMBOL", "ENSEMBL"),
  do_simplify = TRUE,
  do_save_plots = TRUE,
  save_rds = TRUE,
  save_csv = TRUE,
  n_show_dot = 30,
  resdir = "tempGOresults",
  resdir_rds = NULL,
  figdir = NULL,
  ...){
  
  n_groups = length(genelist)
  if (save_rds && is.null(resdir_rds)){
    resdir_rds = file.path(resdir, "res.rds")
    dir.create(resdir_rds, recursive = T)
  }
  if (do_save_plots && is.null(figdir)){
    figdir = file.path(resdir, "figs")
    dir.create(figdir, recursive = T)
  }
  
  for (ont in ont_choices){
    cgo = compareCluster(genelist, 
                         fun=fun,
                         keyType = keyType,
                         OrgDb = OrgDb,
                         ont=ont, 
                         pvalueCutoff=0.2,
                         qvalueCutoff = 0.2)
    
    print(table(cgo@compareClusterResult$Cluster))
    
    ### simplify
    if (do_simplify){
      message("simplifying...(removing redundant terms)")
      scgo = simplify(cgo)
      print(table(scgo@compareClusterResult$Cluster))
    }
    
    # =============[ plot ]==============
    
    tt = sprintf("Top enriched GO terms compared between groups (%s)", ont)
    
    p_dot = dotplot(cgo, showCategory=n_show_dot, title=tt) 
    n_terms = nlevels(p_dot$data$Description)
    ggsave(filename = sprintf("%s/dotCompare%s-%s.pdf", figdir, ont, tag),
           plot = p_dot,
           width = min(6.5 + n_groups * 0.8, 49),
           height = min(1.5 + n_terms * 0.2, 49))
    
    p_dot2 = dotplot(scgo, showCategory=n_show_dot, title=tt) 
    n_terms2 = nlevels(p_dot2$data$Description)
    ggsave(filename = sprintf("%s/_dotCompareSp%s-%s.pdf", figdir, ont, tag),
           plot = p_dot2,
           width = min(6.5 + n_groups * 0.8, 49),
           height = min(1.5 + n_terms2 * 0.2, 49))
    
    # emapplot(cgo) # error...
    # ==============[ save results ]===============
    saveRDS(cgo, sprintf("%s/cgo_%s_%s.rds", resdir_rds, ont, tag))
    saveRDS(scgo, sprintf("%s/scgo_%s_%s.rds", resdir_rds, ont, tag))
    
    ## saving .csv table
    fn_restb = sprintf("%s/compare-%s-%s.csv", resdir, ont, tag)
    restb = cgo@compareClusterResult
    write.csv(restb, 
              fn_restb, 
              row.names = F)
    write.csv(scgo@compareClusterResult, 
              file = sprintf("%s/compareSp-%s-%s.csv", resdir, ont, tag), 
              row.names = F)
    
    message(sprintf("enrich results table saved into\n\t%s", fn_restb, tag))
    
  }
  return()
}


########################################
string2num = function(s, alow){
  # x = as.numeric(s)
  # if (any(is.na(x))){
  tmp = str_split_fixed(s, pattern = "/", n = 2)
  x1 = as.numeric(tmp[, 1]) 
  x2 = as.numeric(tmp[, 2])
  x2 = ifelse(is.na(x2), yes = 1, no = x2)
  x = x1 / x2
  # }
  x
}

################################################################################
" functions for plot "

candi.colors = list(
  "reds" = c("#e8e8e8", "#FA9FB5", "#d7385e", "#821752"),
  "blackreds" = c("#132743", "#d7385e"),
  "default" = c("red", "blue")
)

plot_GO_bar = function(
  df_plt, 
  x = "Count",
  y = "Description",
  colorby = c("pvalue", "p.adjust")[1],
  xtt = "gene count", xrotation = 0,
  title = "Top enriched GO terms",
  fontsize_tt = 12,
  fn = NULL,
  colors = rev(candi.colors[["blackreds"]]),
  wscale = 1, hscale = 1,
  ...){
  # `df_plt` should be formed like the results from "clusterProfiler"
  p_bar = ggplot(
    df_plt, aes(
      x = .data[[x]], 
      y = fct_reorder(.data[[y]], - pvalue), 
      fill = .data[[colorby]])
    ) +
    geom_bar(stat = "identity") + 
    scale_fill_gradient(
      low=colors[1], high=rev(colors)[1], guide=guide_colorbar(reverse = T)) +
    ylab(NULL) + xlab(xtt) +
    ggtitle(title) +
    # theme_classic()
    theme_bw() 
  if (xrotation != 0){hjust = 1}else{hjust = 0}
    p_bar = p_bar + 
      theme(axis.text.x=element_text(angle=xrotation, hjust = hjust), 
            plot.title = element_text(size = fontsize_tt)
            )#size=13, 
  if (! is.null(fn)){
    w = min(6.5 * wscale, 49)
    h = min((1 + dim(df_plt)[1] * 0.125) * hscale, 49)
    ggsave(fn, p_bar, width = w, height = h)
    sprintf("figure saved: %s", fn) %>% message()
  }
  p_bar
}


plot_GO_compare = function(
  df_plt,
  x = NULL, # c("group", "Cluster", "cluster")[1]
  y = "Description",
  title = "Top enriched GO terms",
  sizeby = c("GeneRatio", "Count")[1],
  colorby = c("pvalue", "p.adjust")[1],
  colors = rev(candi.colors[["blackreds"]]),
  xtt = NULL, xrotation = 0,
  fontsize_tt = 12,
  legend.position="right",
  legend.box = c("horizontal", "vertical")[2],
  fn = NULL,
  transpose=F,
  wscale = 1, hscale = 1,
  ...
){
  # `df_plt` should be formed like the results from "clusterProfiler"
  if (is.null(x)){
  candi.x = c("group", "Cluster", "cluster", "Group")
  for (nm in candi.x){
    if (nm %in% colnames(df_plt)){
      x = nm
    }
  }}
  p_dot = ggplot(df_plt) +
    geom_point(aes(
      x = .data[[x]], 
      y = fct_reorder(.data[[y]], as.integer(.data[[x]]), .fun=min, .desc = T),
      color = .data[[colorby]], 
      size = string2num(.data[[sizeby]])
    )) + 
    scale_color_gradient(
      colorby, low=colors[1], high=rev(colors)[1], 
      guide=guide_colorbar(reverse = T)) +
    scale_size_continuous(sizeby) +
    ylab(NULL) + xlab(xtt) +
    ggtitle(title) +
    theme_bw() 
  if (xrotation != 0){hjust = 1}else{hjust = NULL}
  p_dot = p_dot + 
    theme(axis.text.x=element_text(angle=xrotation, hjust = hjust), 
          plot.title = element_text(size = fontsize_tt),
          legend.position=legend.position,
          legend.box = legend.box)#size=13, 
  n_groups = df_plt[[x]] %>% unique() %>% length()
  n_terms = df_plt[[y]] %>% unique() %>% length()
  if (transpose){
    p_dot = p_dot + coord_flip()
    h0 = 1 + n_groups * 0.125
    w0 = 6 + n_terms * 0.7
  }else{
    h0 = 1 + n_terms * 0.125
    w0 = 6 + n_groups * 0.7
  }
  if (! is.null(fn)){
    w = min(w0 * wscale, 49)
    h = min(h0 * hscale, 49)
    ggsave(fn, p_dot, width = w, height = h)
    sprintf("figure saved: %s", fn) %>% message()
  }
  p_dot
}





















