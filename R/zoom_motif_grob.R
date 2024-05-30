# ==============================================================================
# construct grob
# ==============================================================================
#' construct zoom_motif_grob
#'
#' @param pwmMatrix a pwmMatrix from seq2pwmMatrix function, default NULL.
#' @param rel_height the relative height for each part, default c(0.4,0.2,0.4).
#' @param link_gp_col the color of link, default "white".
#' @param link_gp_fill the fill of link, default "grey90".
#' @param border_col the border color, default "white".
#' @param name name
#' @param gp gp
#' @param vp vp
#'
#' @return a grob
#' @export
#' @import grid
zoom_motif_grob <- function(pwmMatrix = NULL,
                            rel_height = c(0.4,0.2,0.4),
                            link_gp_col = "white",
                            link_gp_fill = "grey90",
                            border_col = "white",
                            name = NULL,gp = NULL,vp = NULL){
  gTree(pwmMatrix = pwmMatrix,
        rel_height = rel_height,
        link_gp_col = link_gp_col,
        link_gp_fill = link_gp_fill,
        border_col = border_col,
        name = name,gp = gp,vp = vp,
        cl = "zoom_motif_grob")
}


#' @export
#' @import grid
makeContent.zoom_motif_grob <- function(x){
  # get relative panel y position
  rel_y = cumsum(x$rel_height)

  # get zoom range
  mid <- ncol(x$pwmMatrix$pwm_mat)/2
  n_col <- ncol(x$pwmMatrix$pwm_mat)
  zoom_rg <- c((mid + x$pwmMatrix$shift[1]):(mid + x$pwmMatrix$shift[2]))

  # ============================================================================
  # first layer
  # ============================================================================
  # grid.newpage()
  # pushViewport(viewport(y = unit(rel_y[1],"npc"),height = unit(rel_height[1],"npc"),
  #                       just = "top",
  #                       xscale = c(1,ncol(pfm)))
  # )
  vp1 <- viewport(y = unit(rel_y[3],"npc"),height = unit(x$rel_height[1],"npc"),
                  just = "top",
                  xscale = c(1,n_col))

  zoom_rect1 <-
    rectGrob(x = 100,
             width = abs(range(zoom_rg)[2] - range(zoom_rg)[1]),
             gp = gpar(fill = x$link_gp_fill,col = x$link_gp_col),
             default.units = "native",
             vp = vp1)
  # grid::grid.draw(JunJunZai::seqLogoGrob2(x = pwmMatrix$pwm_mat))
  logo1 <- JunJunZai::seqLogoGrob2(x = x$pwmMatrix$pwm_mat,vp = vp1)
  panel_rect1 <- rectGrob(gp = gpar(fill = "transparent",col = x$border_col),vp = vp1)
  # popViewport()

  # ============================================================================
  # second layer
  # ============================================================================
  # pushViewport(viewport(y = unit(rel_y[2],"npc"),height = unit(rel_height[2],"npc"),
  #                       just = "top",
  #                       xscale = c(1,ncol(pfm)))
  # )
  vp2 <- viewport(y = unit(rel_y[2],"npc"),height = unit(x$rel_height[2],"npc"),
                  just = "top",
                  xscale = c(1,n_col))

  zoom_rect2 <-
    polygonGrob(x = c(1,range(zoom_rg),n_col),
                y = c(0,1,1,0),
                gp = gpar(fill = x$link_gp_fill,col = x$link_gp_col),
                default.units = "native",
                vp = vp2)
  panel_rect2 <- rectGrob(gp = gpar(fill = "transparent",col = x$border_col),vp = vp2)
  # popViewport()

  # ============================================================================
  # third layer
  # ============================================================================
  # pushViewport(viewport(y = unit(rel_y[3],"npc"),height = unit(rel_height[3],"npc"),
  #                       just = "top",
  #                       xscale = c(1,ncol(pfm)))
  # )
  vp3 <- viewport(y = unit(rel_y[1],"npc"),height = unit(x$rel_height[3],"npc"),
                  just = "top",
                  xscale = c(1,n_col))

  # grid::grid.draw(JunJunZai::seqLogoGrob2(x = pwmMatrix$pwm_mat_sub))
  logo2 <- JunJunZai::seqLogoGrob2(x = x$pwmMatrix$pwm_mat_sub,vp = vp3)
  panel_rect3 <- rectGrob(gp = gpar(fill = "transparent",col = x$border_col),vp = vp3)
  # popViewport()


  # return
  setChildren(x,gList(zoom_rect1,logo1,panel_rect1,
                      zoom_rect2,panel_rect2,
                      logo2,panel_rect3))

}
