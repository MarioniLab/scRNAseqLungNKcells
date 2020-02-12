## functions for single-cell RNAseq analyses of lung NK cells

plot_2D_generic <- function(d1, d2, cols='black', cols_fac=NULL, cols_names=levels(cols_fac), shapes_fac=NULL, shapes_names=levels(shapes_fac), xlab, ylab, cex.legend=1, ...){
  ## plots a 2D representation of data with factors determining point colours and shapes
  ## d1 = dimension 1 of the data
  ## d2 = dimension 2 of the data
  ## cols = colors to use - must be as long as levels of cols_fac or colors will be re-used
  ## cols_fac = the factor vector used to determine color categories
  ## cols_names = the names of the levels of cols_fac for the legend (defaults to levels(cols_fac))
  ## shapes_fac = the factor used to determine shape categories
  ## shapes_names = the names of the levels of shapes_fac for the legend (defaults to levels(shapes_fac))
  ## xlab = x-axis label
  ## ylab = y-axis label
  ## ... other arguments to pass to plot or par
  par(oma=c(0,0,0,5), xpd=FALSE, bty='l', cex=1, ...)
  set.seed(100)
  ordered <- sample(length(d1), length(d1), replace=FALSE)
  if((!is.null(cols_fac)) & (!is.null(shapes_fac))){
    plot(d1[ordered], d2[ordered], 
         col=cols[as.numeric(cols_fac)][ordered], 
         pch=as.numeric(shapes_fac)[ordered], 
         xlab=xlab, ylab=ylab, ...)
  }else if(!is.null(cols_fac)){
    plot(d1[ordered], d2[ordered], 
         col=cols[as.numeric(cols_fac)][ordered], 
         pch=20, 
         xlab=xlab, ylab=ylab, ...)
  }else if(!is.null(shapes_fac)){
    plot(d1[ordered], d2[ordered], 
         col=cols, 
         pch=as.numeric(shapes_fac)[ordered], 
         xlab=xlab, ylab=ylab, ...)
  }else if(cols!='black'){
    plot(d1[ordered], d2[ordered], 
         col=cols[ordered], 
         pch=20, 
         xlab=xlab, ylab=ylab, ...)
  }
  if((!is.null(cols_fac)) & (!is.null(shapes_fac))){
    par(xpd=NA)
    legend(par("usr")[2], par("usr")[4], legend=c(levels(cols_fac), levels(shapes_fac)), 
           col=c(cols, rep('black', length(shapes_names))),
           pch=c(rep(20, times=length(cols_names)), c(1:length(shapes_names))),
           bty='n', cex=cex.legend)
  }else if(!is.null(cols_fac)){
    par(xpd=NA)
    legend(par("usr")[2], par("usr")[4], legend=cols_names, col=cols, pch=20, bty='n', cex=cex.legend, pt.cex=3)
  }else if(!is.null(shapes_fac)){
    par(xpd=NA)
    legend(par("usr")[2], par("usr")[4], legend=shapes_names, col='black', pch=c(1:length(shapes_names)), 
           bty='n', cex=cex.legend, pt.cex=3)
  }
}

gene_plot <- function(gene, mat, red.dim, ...){
  ## generates a reduced dimensionality plot from a SingleCellExperiment object colored by specified gene
  ## gene = gene symbol - corresponding to rowData(mat)$Symbol entry
  ## mat = SingleCellExperiment object
  ## red.dim = which reduced dimensionality representation to be plotted
  if(gene %in% rowData(mat)$Symbol){
    cols_g <- viridis(300)[cut(exprs(mat)[which(rowData(mat)$Symbol %in% gene),], breaks=300)]
    plot(reducedDim(mat, type=red.dim), pch=20, col=cols_g, 
         main=paste0(red.dim, ' colored by log normalized expression of ', gene), 
         xlab=paste0(red.dim, ' 1'), ylab=paste0(red.dim, ' 2'), ...)
  } else {
    print(paste0(gene, ' is not expressed'))
  }
}

gene_box <- function(gene, mat, col_cat, ...){
  ## generates a boxplot of expression of a selected gene split by a phenotypic data category from 
  ### a SingleCellExperiment object 
  ## gene = gene symbol - corresponding to rowData(mat)$Symbol entry
  ## mat = SingleCellExperiment object
  ## col_cat = colData column on which to split the data for boxplot
  if(gene %in% rowData(mat)$Symbol){
    boxplot(split(exprs(mat)[which(rowData(mat)$Symbol == gene),], colData(mat)[,col_cat]), 
            main=gene, ylab='log expression', ...)
  } else {
    print(paste0(gene, ' is not expressed'))
  }
}

gene_heatmap <- function(gene_vector, mat, col_cat_df=NA, col_cat_cols=NA, clust_dist='euclidean', 
                         ...){
  ## generates a heatmap from a SingleCellExperiment object, using the list of genes supplied and with
  ### sidebars colored by desired cell characteristics
  ## gene_vector = character vector of gene symbols desired in heatmap
  ## mat = SingleCellExperiment object
  ## col_cat_df = dataframe with the same number of rows as mat has columns, each column is a 
  ### category for the colored sidebars
  ## col_cat_cols = list of named vectors of color schemes for annotation categories, pheatmap 
  ### default used if not specified
  ## clust_dist = distance metric used for clustering rows and columns
  ## ... = other arguments to pass to pheatmap
  
  require(pheatmap)
  
  missing_genes <- gene_vector[!gene_vector %in% rowData(mat)$Symbol]
  if(length(missing_genes) > 0){
    print(paste0('These genes are not in the dataset: ', missing_genes, collapse=' '))
  }
  if(length(missing)==length(gene_vector)){
    print('No genes present in dataset.')
  }else{
    
    
    mat <- mat[which(rowData(mat)$Symbol %in% gene_vector),]
    mat <- mat[!duplicated(rowData(mat)$Symbol),]
    print(rowData(mat)$Symbol)
    rownames(mat) <- rowData(mat)$Symbol
    mat <- as.matrix(exprs(mat))
    colnames(mat) <- c(1:ncol(mat))
    print(dim(mat))
    
    rownames(col_cat_df) <- colnames(mat)
    print(head(col_cat_df))
    
    pheatmap(mat, annotation_col=col_cat_df, annotation_colors=col_cat_cols, 
             clustering_distance_cols=clust_dist_c,
             clustering_distance_rows=clust_dist, 
             show_rownames=TRUE,
             show_colnames=FALSE,
             ...)
  }
}

score_violin <- function(scores, break_factor, break_factor_cols, 
                         score_name="", break_factor_name=""){
  ## generates a violin plot from pre-calculated module scores
  ## score_colname = name of colData column containing the score of interest
  ## break_factor = factor vector for splitting cells into groups 
  ## break_factor_cols = colors for the different factors
  ## score_name = y-axis label
  ## break_factor_name = x-axis label

  require(ggplot2)
  df <- data.frame(score=scores, pheno=break_factor)
  ggplot(df, aes(x=pheno, y=score, fill=break_factor)) + 
    geom_violin(trim=FALSE)+
    scale_fill_manual(values=break_factor_cols) +
    geom_boxplot(width=0.1, fill='white') +
    stat_summary(fun.y="median", geom="point", color='grey35', shape = 95, size=15) +
    labs(title="",x=break_factor_name, y = score_name) +
    theme_classic() +  theme(legend.position = "none") 
}

