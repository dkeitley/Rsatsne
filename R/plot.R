
#' Extracts dataframe of SATSNE coordinates with colData observation.
#' @export
prepareSATSNEPlot <- function(out, A_data, B_data, colour_by = "dataset") {

  names(out) <- c("Y1", "Y2", "costs1", "costs2", "itercosts")

  A_df <- data.frame(out$Y1, row.names = colnames(A_data))
  colnames(A_df) <- c("SATSNE_1", "SATSNE_2")
  A_df$dataset <- "A"

  B_df <- data.frame(out$Y2, row.names = colnames(B_data))
  colnames(B_df) <- c("SATSNE_1", "SATSNE_2")
  B_df$dataset <- "B"

  if(colour_by != "dataset") {
    if(!(colour_by %in% colnames(colData(A_data))) || !(colour_by %in% colnames(colData(B_data)))) {
      stop("Input colour_by is not in colData.")
    } else {
      A_df <- cbind(A_df, colData(A_data)[,colour_by])
      colnames(A_df)[4] <- colour_by

      B_df <- cbind(B_df, colData(B_data)[,colour_by])
      colnames(B_df)[4] <- colour_by
    }
  }

  df <- rbind(A_df, B_df)
  df <- df[sample(nrow(df)),] # Plot points in random order

  return(df)

}


#' Plots SATSNE embedding.
#' @importFrom ggrastr geom_point_rast
#' @export
plotSATSNE <- function(out, A_data, B_data, colour_by = "dataset",
                       point_size = NULL, point_shape = 16,
                       palette = NULL, legend_pos = "right",
                       raster_points = TRUE, raster.dpi = 300) {

  df <- prepareSATSNEPlot(out, A_data, B_data, colour_by)
  p <- ggplot(df, aes_string(x = "SATSNE_1", y = "SATSNE_2", color = colour_by))

  if(is.null(point_size)) { point_size <- nrow(df)*0.000005 }
  if(raster_points) { p <- p + ggrastr::geom_point_rast(raster.dpi = raster.dpi, size = point_size, shape = point_shape)}
  else { p <- p + geom_point(size = point_size, shape = point_shape) }

  if(!is.null(palette)) {
    p <- p + scale_color_manual(values = palette)
  }

  p <- p + theme_linedraw() + theme(legend.position = legend_pos,
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.text=element_blank(),
                                    axis.ticks=element_blank())

  return(p)

}


#' Highlights colData observation values in SATSNE embedding.
#' @importFrom ggrastr geom_point_rast
#' @export
highlightSATSNE <- function(out, A_data, B_data, colour_by = "dataset",
                            highlight = NULL, alpha = 0.6,
                            point_size = NULL, point_shape = 16,
                            palette = NULL, legend_pos = "right",
                            raster_points = TRUE, raster.dpi = 300) {

  df <- prepareSATSNEPlot(out, A_data, B_data, colour_by)

  if(is.null(point_size)) { point_size <- nrow(df)*0.000005 }
  p <- ggplot(df[!(df[[colour_by]] %in% highlight),], aes_string(x = "SATSNE_1", y = "SATSNE_2"))

  if(raster_points) { p <- p + ggrastr::geom_point_rast(raster.dpi = raster.dpi,
                                                        size = point_size, shape = point_shape,
                                                        color = "grey", alpha = alpha)}

  else { p <- p + geom_point(size = point_size, shape = point_shape,
                             color = "grey", alpha = alpha) }

  p <- p + geom_point(data = df[df[[colour_by]] %in% highlight,], aes_string(x = "SATSNE_1", y = "SATSNE_2", color = colour_by),
                      size = point_size, shape = point_shape)

  if(!is.null(palette)) {
    p <- p + scale_color_manual(values = palette)
  }

  p <- p + theme_linedraw() + theme(legend.position = legend_pos,
                                    panel.grid.major = element_blank(),
                                    panel.grid.minor = element_blank(),
                                    axis.text=element_blank(),
                                    axis.ticks=element_blank())

  return(p)

}


.calcEuclidean <- function(x1,y1,x2,y2) (x1-y1)**2 + (x2-y2)**2



#' @export
plotAnnotatedSATSNE <- function(out, A_data, B_data, colour_by="celltype", ncell_filt=5,
                                size=1, label_force=15, text_size = 2, line_size = 0.5,
                                palette=NULL) {

  col_data <- prepareSATSNEPlot(out, A_data, B_data, colour_by)
  col_data$group <- col_data[[colour_by]]

  # Get mean position for each group
  mean_data <- col_data %>% group_by(group) %>% summarize_at(.vars = vars(SATSNE_1,SATSNE_2),.funs = c("mean"))
  mean_data <- as.data.frame(mean_data[complete.cases(mean_data),])
  rownames(mean_data) <- mean_data$group

  # Get position of closest cell to group mean
  label_pos <- col_data %>% group_by(group) %>%  filter(
    .calcEuclidean(SATSNE_1, mean_data[group,"SATSNE_1"], SATSNE_2, mean_data[group,"SATSNE_2"]) ==
      min(.calcEuclidean(SATSNE_1, mean_data[group,"SATSNE_1"], SATSNE_2, mean_data[group,"SATSNE_2"])))

  # Filter annotations with less than ncell_filt cells
  freqs <- table(col_data[[colour_by]])
  label_pos <- label_pos[label_pos$group %in% names(freqs[freqs > ncell_filt]),]

  # Repels labels from rest of points
  col_data$group <- ""
  label_pos <- rbind(label_pos, col_data)

  # Wrap long labels
  label_pos$group_wrapped <- stringr::str_wrap(label_pos$group , width = 10)

  p <- ggplot(col_data, aes_string(x="SATSNE_1",y="SATSNE_2",colour=colour_by)) +
    ggrastr::geom_point_rast(size=size) +
    geom_text_repel(data=label_pos, aes(x=SATSNE_1, y=SATSNE_2,label=group_wrapped,segment.colour=group),color="black",
                    min.segment.length = 0,box.padding = 0.5,max.overlaps=Inf,size=text_size,force=label_force,
                    segment.size = line_size) +
    coord_cartesian(clip = "off") +
    scale_colour_discrete(aesthetics=c("color","segment.colour"),drop=TRUE,
                          breaks=names(freqs[freqs > ncell_filt])) +
    theme_linedraw()+
    theme(legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank())

  return(p)

}


