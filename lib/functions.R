
# R functions


ggstat <- function(data, mapping, comparisons = NULL, type = sd, barwidth = 0.1, shape = 4, size = 4, stroke = 1.5, add_labels = NULL, ...){
  data <- data.frame(data)
  x <- mapping[["x"]]
  y <- mapping[["y"]]
  group_by <- mapping[["x"]]
  type <- rlang::enquo(type)
  
  statsdf <- subset(data, !is.na(rlang::as_name(x)) & !is.na(rlang::as_name(y))) %>% 
    group_by(!!x) %>% summarize(mean = mean(!!y, na.rm = TRUE), var = var(!!y, na.rm = TRUE), sd = sd(!!y, na.rm = TRUE), n = n())
  
  statsdf %<>% mutate(se = sqrt(var/n))
  statsdf %<>% mutate(ymin = mean - !!type, ymax = mean + !!type)
  ggs <- ggplot2::ggplot(data = statsdf,
                         mapping = aes(x = !!x, y = mean)) +
    geom_errorbar(aes(ymin = ymin, max = ymax), width = barwidth) + 
    ggplot2::geom_point(data = data, mapping = mapping, shape = shape, size = shape, stroke = stroke) +
    
    ylab(rlang::as_name(y))
  
  getbw <- function(ggs, barwidth) {
    diff(datamisc:::gg_getLimits(ggs)$x) * barwidth
  }
  
  ggs$layers[[which(datamisc:::gg_getGeoms(ggs) == "GeomErrorbar")]]$geom_params$width <- getbw(ggs, barwidth)
  
  if (!is.null(comparisons)){
    ggs <- get_stat_comparisons(ggs, data, x, y, comparisons, ymax = max(datamisc:::gg_getLimits(ggs)$y), add_labels = add_labels, ...)
  }
  
  
  ggs
}


get_stat_comparisons <- function(gg, data, x, y, comparisons, ymax, add_labels = NULL, ...){
  
  data <- subset(data, !is.na(rlang::as_name(x)) & !is.na(rlang::as_name(y))) 
  
  groups <- split(data, data[[rlang::as_name(x)]]) |> lapply(FUN = dplyr::pull, rlang::as_name(y))
  names(comparisons) <- sapply(comparisons, function(comp) sort(comp) |> paste(collapse = "vs"))
  testres <- sapply(comparisons, function(comp){
    wilcox.test(groups[[comp[1]]], groups[[comp[2]]], ...)$p.value
  })
  
  compdf <- expand.grid(unique(unlist(comparisons)), unique(unlist(comparisons))) |> apply(1, sort) |> t()
  l <- compdf |> apply(1, unique) |> sapply(length)
  compdf <- unique(compdf[l == 2,]) |> data.frame()
  colnames(compdf) <- c("x","xend")
  compdf$comparison <- as.data.frame(t(compdf[,c("xend", "x")])) |> sapply(FUN = function(x) sort(x) |> paste(collapse = "vs")) |> unname()
  compdf$p <- testres[compdf$comparison] |> p.adjust(method = "fdr")
  
  x <- union(compdf$x,compdf$xend)
  xlev <- setNames(as.numeric(factor(x)), x)
  compdf$dist <- abs(xlev[compdf$x] - xlev[compdf$xend])
  
  if (nrow(compdf) > 1){
    compdf$yoffset <- rank(compdf$dist, ties.method = "random") |> scale(center = TRUE) |> as.numeric()
  } else {
    compdf$yoffset <- 1
  }
  
  compdf$y <- ymax
  compdf$y <- compdf$y + compdf$y*0.1*compdf$yoffset
  ydist <- as.numeric(dist(c(min(0, ymax), max(compdf$y))))
  
  compdf$text <- datamisc::pval_format(compdf$p)
  if (!is.null(add_labels)){
    compdf$text <- paste0(compdf$text, " ", add_labels)
  }
  
  compdf$xmean <- rowMeans(cbind(xlev[compdf$x], xlev[compdf$xend]))
  
  gg + geom_segment(data = compdf, mapping = aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
    geom_text(data = compdf, mapping = aes(x = xmean, y = y, label = text), nudge_y = ydist*0.1, inherit.aes = FALSE, size = 14/ggplot2::.pt)  
  
}







scale_colors_custom <- function(colors = NULL, breaks = NULL, breaks.lab = NULL, n = 2, ...){
  
  if (is.null(colors)){
    if (n == 2){
      colors <- c(rgb(0.94,0.94,0.94), "#5106b2")
    } else if (n == 3){
      colors <- c("#5106b2", rgb(0.94,0.94,0.94), "#f2c80e")
    }
  }
  
  if (is.null(breaks)){
    ggplot2::scale_color_gradientn(colors = colors, breaks = breaks.lab, ...)
  } else {
    ggplot2::scale_color_gradientn(colors = colors, values = scales::rescale(c(min(breaks), max(breaks))), limits = breaks, breaks = breaks.lab, ...)
  }
}

hm2g <- function(hm){
  grid::grid.grabExpr(ComplexHeatmap::draw(hm))
}


filter_rownames <- function(data, pattern = NULL, ...){
  data[grepl(pattern, rownames(data), ...),, drop = FALSE]
}


filter_colnames <- function(data, pattern = NULL, ...){
  data[, grepl(pattern, colnames(data), ...), drop = FALSE]
}


getGEOdata <- function(acc, dir = "~/myScratch/GEO", ...){
  stopifnot(requireNamespace("GEOquery"))
  
  geo <- list()
  geo$coldata <- GEOquery::getGEO(acc, destdir = dir)[[1]] |> pData()
  geo$suppfiles <- GEOquery::getGEOSuppFiles(acc, baseDir = dir, ...) |> rownames()
  names(geo$suppfiles) <- basename(geo$suppfiles)
  
  geo
}


logtrans <- function(data, p = 1, base = 2){
  log(data + p, base = base)
}


runTEST <- function(data, design, contrasts, pull = "p.value", FUN = wilcox.test, p.adjust.method = "fdr", ...){
  
  pval <- sapply(contrasts, function(tmpcontrast){
    
    tmpdesign <- design[design[[tmpcontrast[1]]] %in% tmpcontrast[-1],, drop = FALSE]
    tmpdata <- data[,rownames(tmpdesign)]
    
    apply(tmpdata, 1, function(x){
      FUN(x ~ tmpdesign[[tmpcontrast[1]]], ...)[[pull]]
    })
    
  }) |> data.frame()
  
  padj <- lapply(pval, p.adjust, method = p.adjust.method) |> as.data.frame()
  rownames(padj) <- rownames(pval)
  
  list(pval = pval, padj = padj)
}


limmamarkers <- function(data, design, group = "Celltype", formula = ~ group, ...){
  data <- data[,rownames(design)]
  ct <- unique(design[[group]])
  markers <- lapply(setNames(ct, ct), function(ctx){
    design$group <- ifelse(design[[group]] == ctx, ctx, "other")
    res <- runLIMMA(data, design, formula = formula, contrasts = list(de = c("group", ctx, "other")))
    subset(res$de, ...)$id
  })
  markers
}


getmarkers <- function(data, design, group = "Celltype", formula = ~ group, ...){
  ct <- unique(design[[group]])
  markers <- lapply(setNames(ct, ct), function(ctx){
    design$group <- ifelse(design[[group]] == ctx, ctx, "other")
    res <- runDESeq2(data, design, formula = formula, contrasts = list(de = c("group", ctx, "other")))
    subset(res$results[[1]], ...)$gene
  })
  markers
}


cheatmap <- function(data, rowdf = NULL, coldf = NULL, scale = FALSE, cluster_rows = NULL, cluster_cols = NULL,
                      rowdf_side = "left", coldf_side = "top", rowdf_legend = TRUE, coldf_legend = TRUE,
                      row_split = NULL, col_split = NULL,
                      legend_border = "black", anno_border = "black",
                      fontsize = 12, rowcex = NULL, colcex = 1,
                      col.up = "red", col.mid = "white", col.dn = "blue",
                      heatpal = NULL, border = NULL, title = NULL, colors = NULL,
                      inf = F, na = 0, mat = NULL, markoob = FALSE, markshape = 4, marksize = NULL, na_col = "grey", maxchar = 35, legend_params = list(), ...){
  
  
  ### Data ----
  datacall <- substitute(data)
  if (is.null(title)){
    if (grepl("row|col", scale, ignore.case = TRUE)){
      title <- "z-score"
    } else {
      title <- deparse1(datacall)
    }
  }
  if (title == FALSE) title <- " "
  
  heatdata <- eval(datacall, envir = parent.frame())
  heatdata <- data.matrix(heatdata)
  heatdata <- matScale(heatdata, rows = grepl("row", scale, ignore.case = TRUE), cols = grepl("col", scale, ignore.case = TRUE))
  
  ### Clustering ----
  if (is.null(row_split) & is.null(col_split)){
    clust <- clusterData(heatdata, rows = cluster_rows, cols = cluster_cols, inf = inf, na = na)
    if (is.null(clust$rows)) clust$rows <- FALSE
    if (is.null(clust$cols)) clust$cols <- FALSE
  } else {
    if (is.null(cluster_rows)) cluster_rows <- TRUE
    if (is.null(cluster_cols)) cluster_cols <- TRUE
    clust <- list()
    clust$rows <- cluster_rows
    clust$cols <- cluster_cols
  }
  
  
  ### Colors ----
  
  # scale colors
  if (is.null(heatpal)){
    heatpal_colors <- datamisc:::getColorScale(heatdata, col.up = col.up, col.mid = col.mid, col.dn = col.dn)
    heatpal <- circlize::colorRamp2(breaks = heatpal_colors, colors = names(heatpal_colors))
  }
  
  # annotation colors
  docol <- setdiff(unlist(lapply(list(coldf, rowdf), colnames)), names(colors))
  addcol <- NULL
  if (length(docol) > 0){
    if (!is.null(coldf)) all <- coldf
    if (!is.null(rowdf)) all <- rowdf
    if (!is.null(coldf) & !is.null(rowdf)) all <- dplyr::full_join(coldf, rowdf, by = character())
    addcol <- getColors(all[,docol,drop = FALSE])
  }
  
  colors <- c(colors, addcol)
  
  colors <- lapply(colors, function(tmp){
    if (!is.function(tmp)){
      tmp[!is.na(tmp) & !is.na(names(tmp))]
    } else {
      tmp
    }
  })
  
  # cell border
  if (is.null(border)){
    if (nrow(heatdata) < 100 & ncol(heatdata) < 100){
      border <- grid::gpar(col = rgb(1,1,1), lwd = grid::unit(1, "pt"))
    } else {
      border <- grid::gpar(col = NA)
    }
  } else {
    if (any(border == TRUE)){
      border <- grid::gpar(col = rgb(1,1,1), lwd = grid::unit(1, "pt"))
    } else if (length(border) > 1){
      border <- grid::gpar(col = border[is.na(as.numeric(border))], lwd = grid::unit(as.numeric(border)[!is.na(as.numeric(border))], "pt"))
    } else {
      border <- grid::gpar(col = NA)
    }
  }
  
  
  ### Annotations ----
  
  # Legends
  # see ?`color_mapping_legend,ColorMapping-method`
  legend_params <- c(list(title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
                        legend_height = grid::unit(0.2, "npc"),
                        border = legend_border,
                        labels_gp = grid::gpar(fontsize = fontsize)), legend_params)
  
  
  # Row annotation
  rowAnn <- NULL
  if (!is.null(rowdf)) rowAnn <- datamisc:::getCXanno(df = rowdf[rownames(heatdata),, drop = FALSE],
                                           colors = colors,
                                           anno_border = anno_border,
                                           side = rowdf_side,
                                           legend = rowdf_legend,
                                           legend_params = legend_params)
  
  # Column annotation
  colAnn <- NULL
  if (!is.null(coldf)) colAnn <- datamisc:::getCXanno(coldf[colnames(heatdata),, drop = FALSE],
                                           colors = colors,
                                           anno_border = anno_border,
                                           side = coldf_side,
                                           legend = coldf_legend,
                                           legend_params = legend_params)
  
  
  # Cell annotation
  
  if (markoob == TRUE & is.null(mat)){
    mat <- matrix(data = FALSE, nrow = nrow(heatdata), ncol = ncol(heatdata), dimnames = dimnames(heatdata))
    mat[heatdata < min(heatpal_colors)] <- TRUE
    mat[heatdata > max(heatpal_colors)] <- TRUE
  }
  
  if (is.null(marksize)) marksize <- fontsize * 0.6
  cellFUN <- NULL
  if (!is.null(mat)){
    if (!is.logical(mat)) stop("'Mat' must be a logical indicator of whether cells should be marked!")
    cellmat <- mat[rownames(heatdata), colnames(heatdata)]
    cellFUN <- function(j, i, x, y, width, height, fill){
      if (naf(cellmat[i,j] == TRUE)){ grid::grid.points(x, y, pch = markshape, size = unit(marksize, "pt")) }
    }
  }
  
  
  ### Heatmap ----
  
  dimnames(heatdata) <- lapply(dimnames(heatdata), function(x) cutstr(x, maxchar = maxchar))
  
  if (is.null(rowcex) & nrow(heatdata) > 10*ncol(heatdata)) rowcex <- 1/log10(nrow(heatdata))
  if (is.null(rowcex)) rowcex <- 1
  if (is.null(colcex)) colcex <- 1
  
  hm <- ComplexHeatmap::Heatmap(name = title,
                                matrix = heatdata,
                                row_names_max_width = grid::unit(0.3, "npc"),
                                column_names_max_height = grid::unit(0.3, "npc"),
                                column_title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
                                rect_gp = border,
                                na_col = na_col,
                                left_annotation = rowAnn,
                                top_annotation = colAnn,
                                row_names_gp = grid::gpar(fontsize = fontsize * rowcex),
                                column_names_gp = grid::gpar(fontsize = fontsize * colcex),
                                col = heatpal,
                                row_split = row_split,
                                column_split = col_split,
                                heatmap_legend_param = legend_params,
                                cluster_rows = clust$rows,
                                cluster_columns = clust$cols,
                                cell_fun = cellFUN,
                                ...)
  
  hm
}


get_abs_conc <- function(TE, ref){
  mM <- metData(TE)$concentration.of.target.mM |> setNames(rownames(metData(C13medium)))
  nocells <- subset(TE, is.na(Celltype))
  au_medium <- assay(nocells, "norm", type = "met") |> rowMeans()
  
  au_to_mM <- mM / au_medium
  
  assay(TE, "conc", type = "met") <- assay(TE, "norm", type = "met") * au_to_mM
  assay(TE, "conc") <- assay(TE, "clean") * au_to_mM[isoData(TE)$metabolite]
  TE
}


get_exchanges <- function(TE, assay = "norm", type = "met", assay.new = "exchange"){
  
  nocells <- subset(TE, is.na(Celltype))
  
  medium_nocells <- assay(nocells, assay, type = type)
  medium_cells <- assay(TE, assay, type = type)
  
  batch <- colData(TE)[,c("Celltype", "Batch")]
  colnames(medium_nocells) <- batch[colnames(medium_nocells),]$Batch
  exch <- medium_cells - medium_nocells[,batch[colnames(medium_cells),]$Batch]
  
  assay(TE, assay.new, type = type) <- exch
  TE
}


get_test_exchanges <- function(TE, assay = "norm", type = "met", assay.new = "exchange"){
  
  nocells <- subset(TE, is.na(Celltype))
  
  medium_nocells <- assay(nocells, assay, type = type)
  medium_cells <- assay(TE, assay, type = type)
  
  batch <- colData(TE)[,c("Celltype", "Batch")]
  colnames(medium_nocells) <- batch[colnames(medium_nocells),]$Batch
  
  list(medium_cells = medium_cells, medium_nocells = medium_nocells[,batch[colnames(medium_cells),]$Batch])
}


get_top_mets <- function(data, n = Inf, IS = "myr_d27"){
  
  allmeans <- rowMeans(abs(data), na.rm = FALSE)
  allmeans <- allmeans[!is.na(allmeans)]
  
  allmeans <- sort(allmeans, decreasing = TRUE)
  
  data[setdiff(names(allmeans), IS), ] |> head(n = n)
}


