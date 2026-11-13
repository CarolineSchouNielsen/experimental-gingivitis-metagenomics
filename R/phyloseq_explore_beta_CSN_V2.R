# Modified from Florentin Constancias' original function found here: https://raw.githubusercontent.com/fconstancias/KU_Caroline/refs/heads/main/code/functions/phyloseq_functions.R
# The function have been modified at four graphical parameters: 
# Line 127. Point size have been decreased from 3 to 2.25.
# Line 129. Arrow thickness have been increased from 0.05 to 0.1. 
# Line 129. Arrow linetype has been changed from "longdash" to "dashed".
# Line 275. Line color has been changed from "grey50" to "grey30". 
phyloseq_explore_beta <- function(ps_up = ps_up,
                                  beta = beta,
                                  m = "PCoA",
                                  distance_for_more = "rAitchison",
                                  color_group = "Time",
                                  shape_group = "Sample",
                                  tax_rank_fit  = "Genus",
                                  metadata_sel_envfit = c("delta_mean_plaque", "delta_mean_bleeding"),
                                  pval_cutoff = 0.05,
                                  top_r_envfit = 15,
                                  adj_method = "fdr",
                                  alpha = NULL,
                                  col_pal = time_pal,
                                  fill_pal = time_pal,
                                  path_group = "interaction(Sample,Subject)",
                                  facet_formula = "Sample ~ .",
                                  axis1 = 1,
                                  axis2 = 2,
                                  seed = 123,
                                  permanova_terms = c("Time", "cluster_Dtp2"),
                                  metadata_dist_boxplot = NULL,
                                  plot_dist_boxplot = TRUE,
                                  strata = "none",
                                  run_phyloseq_TW = TRUE,
                                  run_pairwise_permanova = TRUE,
                                  perform_pcoa_per_variable = TRUE,
                                  perm = 999)
{
  # 1. Robust runtime override of pairwise.adonis2 to force flat vector extraction
  pairwise.adonis2_robust <- function(x, data, strata = NULL, nperm=999, ... ) {
    ststri <- ifelse(is.null(strata),'Null',strata)
    fostri <- as.character(x)
    x1 <- x
    lhs <- eval(x1[[2]], environment(x1), globalenv())
    environment(x1) <- environment()
    rhs <- x1[[3]]
    x1[[2]] <- NULL
    rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE)
    co <- combn(unique(as.character(rhs.frame[,1])),2)
    nameres <- c('parent_call')
    for (elem in 1:ncol(co)){
      nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
    }
    res <- vector(mode="list", length=length(nameres))
    names(res) <- nameres
    res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3], ', strata =',ststri, ', permutations',nperm ))
    for(elem in 1:ncol(co)){
      if(inherits(eval(lhs),'dist')){
        xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]), rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
      }else{
        xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
      }
      mdat1 <- data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
      if(length(rhs) == 1){
        xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))
      }else{
        xnew <- as.formula(paste('xred' , paste(rhs[-1],collapse= as.character(rhs[1])), sep='~'))}
      if(is.null(strata)){
        ad <- adonis2(xnew,data=mdat1, ... )
      }else{
        perm <- how(nperm = nperm)
        setBlocks(perm) <- mdat1[[ststri]] # Force vector extraction
        ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )
      }
      res[nameres[elem+1]] <- list(ad[1:5])
    }
    class(res) <- c("pwadstrata", "list")
    return(res)
  }
  
  # Register robust function in global environment and package namespace if loaded
  assign("pairwise.adonis2", pairwise.adonis2_robust, envir = .GlobalEnv)
  if ("pairwiseAdonis" %in% loadedNamespaces()) {
    try(assignInNamespace("pairwise.adonis2", pairwise.adonis2_robust, ns = "pairwiseAdonis"), quietly = TRUE)
  }
  
  suppressMessages({
    suppressWarnings({
      
      out <- list() # Initialize output list
      
      #### 0. Sanitize Metadata to Base R data.frame
      if (!is.null(sample_data(ps_up, errorIfNULL = FALSE))) {
        meta_clean <- as.data.frame(as(sample_data(ps_up), "data.frame"))
        class(meta_clean) <- "data.frame"
        sample_data(ps_up) <- sample_data(meta_clean)
      }
      
      #### 1. Compute beta diversity and organize data
      # Subset each distance matrix in beta to match the subsetted ps_up samples
      if (is.list(beta)) {
        samples_to_keep <- sample_names(ps_up)
        beta <- lapply(beta, function(d) {
          d_mat <- as.matrix(d)
          common_samples <- intersect(samples_to_keep, colnames(d_mat))
          as.dist(d_mat[common_samples, common_samples, drop = FALSE])
        })
      }
      
      #### 2. Ordination
      ps_up %>% 
        phyloseq_plot_bdiv(dlist = beta, ps_rare = ., m = m, seed = seed, axis1 = axis1, axis2 = axis2) -> out$plot_list
      
      # Generate PCoA ordination plots with specified color and shape groups
      out$plot_list %>%
        phyloseq_plot_ordinations_facet(color_group = color_group, shape_group = shape_group, alpha = alpha) +
        scale_color_manual(name = "", values = col_pal, na.value = "black") +
        scale_fill_manual(name = "", values = fill_pal, na.value = "black") + 
        theme_linedraw() + theme(legend.position = "right") +
        facet_wrap(distance ~ ., scales = "free", nrow = 3) -> out$pcoas
      
      # Extract explained variance from ordination plot
      out$plot_list %>%
        phyloseq_ordinations_expl_var() %>% 
        mutate_if(is.numeric, ~round(., 3)) -> out$expl_var
      
      #### 3. Detailed Ordination for distance_for_more selections
      out$plot_list[[distance_for_more]] %>% 
        gginnards::delete_layers(. , "GeomPoint") -> out$plot_list[[distance_for_more]] 
      
      out$plot_list[[distance_for_more]] + geom_point(size = 2.25, aes_string(colour = color_group, shape = shape_group, alpha = alpha)) +
        geom_path(data = out$plot_list[[distance_for_more]]$data %>% arrange(Subject), aes_string(group = path_group),
                  arrow = arrow(angle = 30, length = unit(0.15, "inches"), ends = "last", type = "open"), linetype = "dashed", size = 0.1) +
        theme_light() + scale_color_manual(name = "", values = col_pal, na.value = "black") +
        scale_fill_manual(name = "", values = fill_pal, na.value = "black") +
        facet_grid(as.formula(facet_formula), scales = "fixed", space = "fixed", switch = "y") + theme_linedraw() + theme(strip.placement = "outside") + coord_fixed() -> out$PCOA
      
      # Extract and hide legend for the ordination plot
      out$PCOA %>% 
        get_plotandlegend() -> out$PCOA
      
      if(length(top_r_envfit) > 0){
        #### 4. Env Fit (optional)
        phyloseq_add_taxa_vector_fix(dist = beta[[distance_for_more]],
                                     phyloseq =  ps_up ,
                                     figure_ord = out$PCOA$p + facet_null(),
                                     m = "PCoA",
                                     pval_cutoff = pval_cutoff,
                                     top_r = top_r_envfit,
                                     tax_rank_plot = tax_rank_fit,
                                     vector_color = "green4",
                                     taxnames_rm = c("unknown", "Unclassified"),
                                     seed = seed,
                                     perm = perm,
                                     adj_method = adj_method,
                                     transform = "compositional") -> env_fit_tax
        
        out$env_fit_tax <- env_fit_tax
      }
      if(!is.null(metadata_sel_envfit)){
        phyloseq_add_metadata_vector_fix(dist = beta[[distance_for_more]],
                                         phyloseq = ps_up ,
                                         figure_ord = out$PCOA$p + facet_null(),
                                         m = "PCoA",
                                         pval_cutoff = pval_cutoff,
                                         top_r = top_r_envfit,
                                         metadata_sel = metadata_sel_envfit,
                                         vector_color = "salmon3",
                                         fact = 1,
                                         seed = seed,
                                         perm = perm,
                                         adj_method = adj_method,
                                         norm_method = "center_scale",
                                         na_rm = TRUE) -> env_fit_meta
        
        out$env_fit_meta <- env_fit_meta
      }
      
      #### 5. PERMANOVA Analysis
      if(!is.null(permanova_terms)){
        
        tmp_out1 = NULL; tmp_out2 = NULL
        for (terms_margins in c("terms", "margin")) {
          form = paste0(permanova_terms, collapse = " * ")
          ps_up %>% lapply(beta, FUN = phyloseq_adonis2, physeq = ., formula = form, nrep = perm, terms_margins = terms_margins) %>%
            bind_rows(.id = "Distance") %>% mutate("Group" = as.vector(form), "terms_margin" = terms_margins) %>%
            bind_rows(., tmp_out1) -> tmp_out1
          
          form = paste0(permanova_terms, collapse = " + ")
          
          # Dynamically resolve strata into a permute::how() block design passed directly to nrep
          nrep_arg <- perm
          strata_arg <- strata
          if (is.character(strata) && length(strata) == 1 && strata != "none") {
            if (strata %in% colnames(sample_data(ps_up))) {
              block_vector <- as.factor(as.character(sample_data(ps_up)[[strata]]))
              nrep_arg <- permute::how(blocks = block_vector, nperm = perm)
              strata_arg <- "none" # Use "none" to bypass the wrapper's subsetting logic safely
            }
          }
          
          ps_up %>% lapply(beta, FUN = phyloseq_adonis2, physeq = ., formula = form, nrep = nrep_arg, terms_margins = terms_margins, strata = strata_arg) %>%
            bind_rows(.id = "Distance") %>% mutate("Group" = as.vector(form), "terms_margin" = terms_margins) %>%
            bind_rows(., tmp_out2) -> tmp_out2
        }
        bind_rows(tmp_out1, tmp_out2) -> out$perm
      }
      
      if(isTRUE(run_pairwise_permanova)){
        #### 6. Pairwise PERMANOVA
        for (compare_header in permanova_terms) {
          sample_data_df <- as.data.frame(sample_data(ps_up))
          
          if (is.factor(sample_data_df[[compare_header]]) || is.character(sample_data_df[[compare_header]])) {
            ps_up %>% 
              lapply(beta, FUN = physeq_pairwise_permanovas_adonis2, physeq = ., compare_header = compare_header, n_perm = perm, strata = strata, terms_margins = "terms") %>%
              bind_rows(.id = "Distance") %>% 
              mutate("Group" = as.vector(compare_header)) %>% 
              bind_rows(., out$pw_perm) -> out$pw_perm
            
            if (run_phyloseq_TW) {
              lapply(beta, FUN = phyloseq_TW, physeq = ps_up, variable = compare_header, nrep = perm) %>%
                bind_rows(.id = "Distance") %>% 
                mutate("Group" = as.vector(compare_header)) %>% 
                bind_rows(., out$tw_perm) -> out$tw_perm
            }
          }
        }
      }
      
      if(isTRUE(perform_pcoa_per_variable))
      {
        #### 7. Generate PCoA plots for each variable in permanova_terms
        for (compare_header in permanova_terms) {
          print(compare_header)
          phyloseq_generate_pcoa_per_variables(tmp = ps_up, group = compare_header, m = "PCoA", dist = beta,
                                               color_group = color_group, shape_group = shape_group, alpha = alpha,
                                               col_pal = col_pal, fill_pal = fill_pal) -> out[[compare_header]]
        }
      }
      
      ####-------- distance boxplot
      if(!is.null(metadata_dist_boxplot)){
        lapply(beta, FUN = phyloseq_distance_boxplot, p = ps_up, d = "Sample") -> dist_bx
        
        ps_up %>%
          sample_data()%>%
          data.frame() %>%
          dplyr::select(any_of(metadata_dist_boxplot)) %>%
          rownames_to_column("Var") %>% 
          drop_na() -> meta_sel
        
        meta_sel_2 <- meta_sel
        colnames(meta_sel) <- paste0(colnames(meta_sel) , "_1")
        colnames(meta_sel_2) <- paste0(colnames(meta_sel_2) , "_2")
        
        dist_df = NULL
        for (d in names(dist_bx)) {
          dist_bx[[d]]$matrix %>% 
            mutate("Distance" = d) %>% 
            bind_rows(.,dist_df) -> dist_df
        }
        
        dist_df %>% 
          left_join(meta_sel, by = c("Var1" = "Var_1")) %>% 
          left_join(meta_sel_2, by = c("Var2" = "Var_2")) -> dist_df
        
        out$dist_df <- dist_df
        
        if(isTRUE(plot_dist_boxplot)){
          dist_df %>% 
            dplyr::filter(Sample_1 == Sample_2,
                          Time_1 == "TP1",
                          Subject_1 == Subject_2) %>% 
            arrange(Distance, Subject_2, Time_2) %>% 
            ggplot(data = ., aes_string(x="Time_2", y="value")) +
            geom_boxplot(outlier.colour = NA, alpha=0.7, aes_string(fill = "Time_2")) +
            geom_jitter(size=1, position = pd, aes_string(color = "Subject_1")) +
            geom_line(aes_string(group="Subject_1"), position = pd, linetype = "dashed", color = "grey30", linewidth = 0.08) +
            facet_grid(as.formula(paste0("Distance ~ Sample_1")), scales = "free_y", space = "fixed", switch = "y") +
            scale_color_manual(name = "", values = sub_pal, na.value = "black") +
            scale_fill_manual(name = "", values = time_pal, na.value = "black") +
            theme_light() + ylab("Distance to Baseline") + xlab(NULL) + 
            theme(axis.text.x = element_blank()) + theme_linedraw() + theme(strip.placement = "outside") -> dist_box
          
          dist_box %>% 
            ggpubr::get_legend(.) %>% 
            ggpubr::as_ggplot(.) -> dist_box_leg
          
          dist_box + theme(legend.position = "none") -> dist_box
          
          out$dist_box <- dist_box
          out$dist_box_leg <- dist_box_leg
        }
      }
      
      return(out)
    })
  })
}
