
fgsea_multiple_samples <- function(ranks, pathways){
    # generate score from fgsea function for each sample in ranks list
    # fgsea requires single gene rank differences vector
    ES_df_list <- list()
    print('------Running fGSEA--------')
    for (sample in names(ranks)) {
        sample_ranks <- ranks %>% select(sample)
        sample_ranks <- sample_ranks[!is.na(sample_ranks), , drop = F]
        print(dim(sample_ranks))
        enrichment_score <- fgsea(pathways,
            setNames(sample_ranks[, sample], rownames(sample_ranks)),,
            scoreType = "std", minSize = 10, nPermSimple = 5000)
        ES_df_list[[sample]] <- enrichment_score
    }
    return(ES_df_list)
}


#-----------------FGSEA wrapper function-----------------
fgsea_wrapper <- function(output_dir, ranks_input, metadata, genesets, analysis_name, input_type, tool_colour){

    annot_colouring <- data.frame(row.names = colnames(ranks_input)) %>%
        mutate(annotation = metadata$annotation[match(colnames(ranks_input), metadata[, c("filename", "annotation")]$filename)]) %>% 
        dplyr::arrange(annotation)
    
    ranks_input_ordered <- ranks_input %>% 
        dplyr::select(rownames(annot_colouring))

    annotation_list <- unique(metadata$annotation)
    NES_scores_pvalfilt_full <- data.frame(GOI_set = sort(names(genesets)))
    NES_scores_full <- data.frame(GOI_set = sort(names(genesets)))
    PADJ_scores_full <- data.frame(GOI_set = sort(names(genesets)))

    invisible(lapply(c('ES_matrix', 'ES_plots'), function(x) dir.create(file.path(output_dir, x), recursive = TRUE)))

    # run for single annotation groups for direct comparison of scores and plotting
    for (annotation_group in annotation_list){
        
        annotation_group_og <- annotation_group
        # simpler/shorter name saving for later (also avoid file saving issue with spaces)
        if (grepl("\\(", annotation_group) & grepl("\\)", annotation_group)) {
            annotation_group <- sub(".*\\((.*)\\).*", "\\1", annotation_group)
        }
        annotation_group <- gsub(' ', '_', gsub('/', '-', gsub(' / ', '-', annotation_group)))
        print(annotation_group)

        metadata_of_annotation <- metadata %>% 
            filter(annotation == annotation_group_og)
        ranks_input_filt <- ranks_input_ordered[, intersect(colnames(ranks_input_ordered), metadata_of_annotation$filename), drop = FALSE]

        # fgsea scoring 
        fgsea_scores <- fgsea_multiple_samples(ranks_input_filt, genesets)
        if (length(fgsea_scores) == 0) {
            print(paste('-----SKIPPING:', annotation_group, 'no scores computed-----'))
            next
        }

        #plot heatmap of NES scores
        GOI_sets_mod <- subset(genesets, names(genesets) %in% intersect(names(genesets),  fgsea_scores[[1]][,pathway]))
        # create score df per gene set across samples
        NES_scores_pvalfilt <- data.frame(GOI_set = names(GOI_sets_mod))
        NES_scores <- data.frame(GOI_set = names(GOI_sets_mod))
        PADJ_scores <- data.frame(GOI_set = names(GOI_sets_mod))

        for (i in 1:length(fgsea_scores)) {
            NES_scores_pvalfilt <- cbind(NES_scores_pvalfilt, fgsea_scores[[i]][, "NES"])
            NES_scores <- cbind(NES_scores, fgsea_scores[[i]][, "NES"])
            PADJ_scores <- cbind(PADJ_scores, fgsea_scores[[i]][, "padj"])
            colnames(NES_scores_pvalfilt) <- c(colnames(NES_scores_pvalfilt)[1:i], names(fgsea_scores)[i])
            colnames(NES_scores) <- c(colnames(NES_scores)[1:i], names(fgsea_scores)[i])
            colnames(PADJ_scores) <- c(colnames(PADJ_scores)[1:i], names(fgsea_scores)[i])
        }
        NES_scores_pvalfilt <- column_to_rownames(NES_scores_pvalfilt, "GOI_set")
        NES_scores <- column_to_rownames(NES_scores, "GOI_set")
        PADJ_scores <- column_to_rownames(PADJ_scores, "GOI_set")

        NES_scores_pvalfilt <- NES_scores_pvalfilt[order(rownames(NES_scores_pvalfilt)), , drop = FALSE] # genesets in alphabetical order
        NES_scores <- NES_scores[order(rownames(NES_scores)), , drop = FALSE]
        PADJ_scores <- PADJ_scores[order(rownames(PADJ_scores)), , drop = FALSE]

        write.table(rownames_to_column(PADJ_scores, var = 'GOI_Set'),
                    file.path(output_dir, 'ES_matrix', paste0(annotation_group, "-PADJ.tsv")),
                    sep = "\t", quote = FALSE, row.names = FALSE)

        # Filter the values of NES (converted to NA) that do not have a padj value < 0.05
        padj_mask <- PADJ_scores <= 0.05          # TRUE where significant
        padj_mask[is.na(padj_mask)] <- FALSE      # treat NA padj as non-significant
        NES_scores_pvalfilt[!padj_mask] <- NA

        if (any(padj_mask)) {
            plot_GSEApheatmap_wNAs(NES_scores_pvalfilt,
                file.path(output_dir, 'ES_plots', paste0(annotation_group, '.png')),
                paste0(annotation_group_og, " - fgsea_", input_type),
                tool_colour)
        }
        # add scores for annotation group to evergrowing full matrix
        NES_scores_pvalfilt <- rownames_to_column(NES_scores_pvalfilt, var = 'GOI_set')
        NES_scores <- rownames_to_column(NES_scores, var = 'GOI_set')
        PADJ_scores <- rownames_to_column(PADJ_scores, var = 'GOI_set')

        NES_scores_pvalfilt_full <- merge(NES_scores_pvalfilt_full, NES_scores_pvalfilt, by = "GOI_set", all = TRUE)
        NES_scores_full <- merge(NES_scores_full, NES_scores, by = "GOI_set", all = TRUE)
        PADJ_scores_full <- merge(PADJ_scores_full, PADJ_scores, by = "GOI_set", all = TRUE)

    }
    NES_scores_pvalfilt_full <- column_to_rownames(NES_scores_pvalfilt_full, var = 'GOI_set')
    NES_scores_full <- column_to_rownames(NES_scores_full, var = 'GOI_set')
    PADJ_scores_full <- column_to_rownames(PADJ_scores_full, var = 'GOI_set')

    NES_scores_pvalfilt_full <- NES_scores_pvalfilt_full[order(rownames(NES_scores_pvalfilt_full)), , drop = FALSE]
    NES_scores_full <- NES_scores_full[order(rownames(NES_scores_full)), , drop = FALSE]
    PADJ_scores_full <- PADJ_scores_full[order(rownames(PADJ_scores_full)), , drop = FALSE]


    write.table(rownames_to_column(NES_scores_pvalfilt_full, var = 'GOI_Set'),
                file.path(output_dir, paste0(analysis_name, "-fullNES_padjfilt.tsv")),
                sep = '\t', quote = FALSE, row.names = FALSE)
    write.table(rownames_to_column(NES_scores_full, var = 'GOI_Set'),
                file.path(output_dir, paste0(analysis_name, "-fullNES.tsv")),
                sep = '\t', quote = FALSE, row.names = FALSE)
    write.table(rownames_to_column(PADJ_scores_full, var = 'GOI_Set'),
                file.path(output_dir, 'ES_matrix', paste0(analysis_name, "-PADJ.tsv")),
                sep = '\t', quote = FALSE, row.names = FALSE)
    
    # plot without legend
    plot_GSEApheatmap_wNAs(NES_scores_full,
        file.path(output_dir, 'ES_plots', paste0(analysis_name, ".png")),
        paste0(analysis_name, " - fgsea_", input_type),
        tool_colour, wannotation = annot_colouring)
    # with legend
    plot_GSEApheatmap_wNAs(NES_scores_full,
        file.path(output_dir, 'ES_plots', paste0(analysis_name, "_annot.png")),
        paste0(analysis_name, " - fgsea_", input_type),
        tool_colour, wannotation = annot_colouring, wlegend = TRUE)
}

plot_GSEApheatmap_wNAs <- function(ES_matrix, png_name, plot_title, tool_colour, wannotation = NA, wlegend = FALSE){
    # function to plot all GSEA pheatmaps the same way, handling NAs in ES matrix outputs
    ES_matrix[ES_matrix == "---"] <- NA
    for (i in seq_along(ES_matrix)) {
        ES_matrix[[i]] <- suppressWarnings(as.numeric(ES_matrix[[i]]))
    }
    
    ES_matrix <- as.matrix(ES_matrix)
    
    scale_columns <- function(mat) {
        out <- mat  # keep original values unless scaling is appropriate
        
        for (j in seq_len(ncol(mat))) {
            col <- mat[, j]
            
            # if all values are NA → leave column untouched
            if (all(is.na(col))) next
            
            m <- mean(col, na.rm = TRUE)
            s <- sd(col, na.rm = TRUE)
            
            # If no variation: leave column as-is
            if (is.na(s) || s == 0) {
                next
            } else {
                out[, j] <- (col - m) / s
            }
        }
        out
    }
    
    scaled_ES <- scale_columns(ES_matrix)

    # --- Create guaranteed-unique breaks ---
    data_min <- min(scaled_ES, na.rm = TRUE)
    data_max <- max(scaled_ES, na.rm = TRUE)

    if (data_min == data_max) {
        data_min <- data_min - 0.1
        data_max <- data_max + 0.1
    }

    myColor <- colorRampPalette(c("white", "white", tool_colour))(101)
    breaks <- seq(data_min, data_max, length.out = length(myColor) + 1)

    png(png_name)
    pheatmap::pheatmap(scaled_ES,
                    show_rownames = T, show_colnames = F,
                    treeheight_row = 0, treeheight_col = 0,
                    cluster_cols = F, cluster_rows = F, scale = 'none',
                    color = myColor, breaks = breaks,
                    main = plot_title, fontsize_row = 5,
                    annotation_col = wannotation,
                    annotation_legend = wlegend)
    dev.off()
}
