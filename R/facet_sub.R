# Creating a layout specification
layout <- function(data, params) {
  data.frame(PANEL = c(1L, 2L), SCALE_X = 1L, SCALE_Y = 1L)
}

# Mapping data into panels
mapping <- function(data, layout, params) {
  if (is.null(data) || nrow(data) == 0) {
    return(cbind(data, PANEL = integer(0)))
  }
  tmp <- cbind(data, PANEL = 1L)
  # tmp2 <- tmp[1,]
  # tmp2[] <- 0
  # tmp2$PANEL <- 2L
  # print(rbind(tmp,tmp2))
  # rbind(tmp,tmp2)
}


# Laying out the panels
render <- function(panels, layout, x_scales, y_scales, ranges, coord, data,
                   theme, params) {
  # Place panels according to settings
  if (params$position %in% c("left","right")) {
    # Put panels in matrix and convert to a gtable
    if(params$position == "right"){
      panels <- matrix(panels, ncol = 2)
    }else if(params$position == "left"){
      panels <- matrix(c(panels[2],panels[1]), ncol = 2)
    }


    panel_table <- gtable::gtable_matrix("layout", panels,
                                         widths = unit(params$panelWidth, "null"),
                                         heights = unit(1, "null"), clip = "on")
    # Add spacing according to theme
    panel_spacing <- if (is.null(theme$panel.spacing.x)) {
      theme$panel.spacing
    } else {
      theme$panel.spacing.x
    }
    panel_table <- gtable::gtable_add_col_space(panel_table, panel_spacing)
  } else {
    # check position
    if(params$position == "bottom"){
      panels <- matrix(c(panels[1],panels[2]), ncol = 1)
    }else if(params$position == "top"){
      panels <- matrix(c(panels[2],panels[1]), ncol = 1)
    }


    panel_table <- gtable::gtable_matrix("layout", panels,
                                         widths = unit(1, "null"),
                                         heights = unit(params$panelWidth, "null"), clip = "on")
    panel_spacing <- if (is.null(theme$panel.spacing.y)) {
      theme$panel.spacing
    } else {
      theme$panel.spacing.y
    }
    panel_table <- gtable::gtable_add_row_space(panel_table, panel_spacing)
  }

  # Name panel grobs so they can be found later
  panel_table$layout$name <- paste0("panel-", c(1, 2))

  # ============================================================================
  # Construct the axes
  # ============================================================================
  axes <- render_axes(ranges[1], ranges[1], coord, theme,
                      transpose = TRUE)

  # Add axes around each panel
  panel_pos_h <- panel_cols(panel_table)$l
  panel_pos_v <- panel_rows(panel_table)$t
  axis_width_l <- unit(grid::convertWidth(
    grid::grobWidth(axes$y$left[[1]]), "cm", TRUE), "cm")
  axis_width_r <- unit(grid::convertWidth(
    grid::grobWidth(axes$y$right[[1]]), "cm", TRUE), "cm")

  ## We do it reverse so we don't change the position of panels when we add axes
  for (i in rev(panel_pos_h)) {
    panel_table <- gtable::gtable_add_cols(panel_table, axis_width_r, i)
    panel_table <- gtable::gtable_add_grob(panel_table,
                                           rep(axes$y$right, length(panel_pos_v)), t = panel_pos_v, l = i + 1,
                                           clip = "off")

    if(params$position %in% c("top","bottom")){
      panel_table <- gtable::gtable_add_cols(panel_table, axis_width_l, i - 1)

      if(params$position == "top"){
        len <- 1; t <- panel_pos_v[2]
      }else if(params$position == "bottom"){
        len <- 1; t <- panel_pos_v[1]
      }else{
        len <- length(panel_pos_v); t <- panel_pos_v
      }

      panel_table <- gtable::gtable_add_grob(panel_table,
                                             rep(axes$y$left, len), t = t, l = i,
                                             clip = "off")
    }else{
      if(params$position == "left"){
        itmp <- panel_pos_h[2]
      }else{
        itmp <- panel_pos_h[1]
      }

      panel_table <- gtable::gtable_add_cols(panel_table, axis_width_l, itmp - 1)
      panel_table <- gtable::gtable_add_grob(panel_table,
                                             rep(axes$y$left, length(panel_pos_v)), t = panel_pos_v, l = itmp,
                                             clip = "off")
      break
    }
  }


  ## Recalculate as gtable has changed
  panel_pos_h <- panel_cols(panel_table)$l
  panel_pos_v <- panel_rows(panel_table)$t
  axis_height_t <- unit(grid::convertHeight(
    grid::grobHeight(axes$x$top[[1]]), "cm", TRUE), "cm")
  axis_height_b <- unit(grid::convertHeight(
    grid::grobHeight(axes$x$bottom[[1]]), "cm", TRUE), "cm")


  for (i in rev(panel_pos_v)) {
    if(params$position %in% c("right","left")){
      panel_table <- gtable::gtable_add_rows(panel_table, axis_height_b, i)

      if(params$position == "right"){
        len <- 1; l <- panel_pos_h[1]
      }else if(params$position == "left"){
        len <- 1; l <- panel_pos_h[2]
      }else{
        len <- length(panel_pos_h); l <- panel_pos_h
      }

      panel_table <- gtable::gtable_add_grob(panel_table,
                                             rep(axes$x$bottom, len), t = i + 1, l = l,
                                             clip = "off")
    }else{
      if(params$position == "top"){
        itmp <- panel_pos_v[2]
      }else{
        itmp <- panel_pos_v[1]
      }

      panel_table <- gtable::gtable_add_rows(panel_table, axis_height_b, itmp)
      panel_table <- gtable::gtable_add_grob(panel_table,
                                             rep(axes$x$bottom, length(panel_pos_h)), t = itmp + 1, l = panel_pos_h[1],
                                             clip = "off")

      break
    }

    panel_table <- gtable::gtable_add_rows(panel_table, axis_height_t, i - 1)
    panel_table <- gtable::gtable_add_grob(panel_table,
                                           rep(axes$x$top, length(panel_pos_h)), t = i, l = panel_pos_h,
                                           clip = "off")
  }

  # ============================================================================
  # Add grobs in sub panel
  # ============================================================================
  if(params$position %in% c("left","top")){
    panel_table$grobs[1] <- list(gTree(children = params$addGrob))
  }else if(params$position %in% c("right","bottom")){
    panel_table$grobs[2] <- list(gTree(children = params$addGrob))
  }

  # ============================================================================
  # Add strips
  # ============================================================================
  strips <- render_strips(
    x = data.frame(name = params$stripTitle),
    labeller = label_value, theme = theme)

  panel_pos_h <- panel_cols(panel_table)$l
  panel_pos_v <- panel_rows(panel_table)$t
  strip_height <- unit(grid::convertHeight(
    grid::grobHeight(strips$x$top[[1]]), "cm", TRUE), "cm")

  for (i in rev(seq_along(panel_pos_v))) {
    panel_table <- gtable::gtable_add_rows(panel_table, strip_height, panel_pos_v[i] - 1)
    if (params$position %in% c("left","right")) {
      panel_table <- gtable::gtable_add_grob(panel_table, strips$x$top,
                                             t = panel_pos_v[i], l = panel_pos_h, clip = "off")
    } else {
      panel_table <- gtable::gtable_add_grob(panel_table, strips$x$top[i],
                                             t = panel_pos_v[i], l = panel_pos_h, clip = "off")
    }
  }


  panel_table
}


# Constructor: shrink is required to govern whether scales are trained on
# Stat-transformed data or not.
#' facet_sub to produce a sub panel around main panel
#'
#' @param position the position of sub panel, one of the "left","right","top" and "bottom".
#' @param shrink shrink
#' @param addGrob the grobs to be added into sub panel, default grid::gList(rectGrob()).
#' @param stripTitle the strip title, default "Main panel" and "Sub panel".
#' @param panelWidth the relative width of main and sub panel, default c(1,1).
#'
#' @export
#' @import ggplot2
#' @import grid
facet_sub <- function(position = c("left","right","top","bottom"),
                      shrink = TRUE,
                      addGrob = grid::gList(rectGrob()),
                      stripTitle = c("Main panel", "Sub panel"),
                      panelWidth = c(1,1)) {
  ggplot2::ggproto(NULL, FacetSub,
                   shrink = shrink,
                   params = list(
                     position = position,
                     addGrob = addGrob,
                     stripTitle = stripTitle,
                     panelWidth = panelWidth
                   )
  )
}


FacetSub <- ggplot2::ggproto("FacetSub", ggplot2::Facet,
                             compute_layout = layout,
                             map_data = mapping,
                             draw_panels = render)


