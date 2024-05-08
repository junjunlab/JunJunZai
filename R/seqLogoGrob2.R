pfm2ic <- getFromNamespace("pfm2ic","monaLisa")
addLetter <- getFromNamespace("addLetter","monaLisa")
annoSeqlogo <- getFromNamespace("annoSeqlogo","monaLisa")

#' @title Create a simple sequence logo grob(slightly modified with seqLogoGrob).
#'
#' @description Create a simple sequence logo grob (grid-graphics object) for a
#'     transcription factor from a position frequency matrix. The logo drawing
#'     code is a simplified version from \code{\link[seqLogo]{seqLogo}} and for
#'     example can be used to embedd sequence logos within other plots.
#'
#' @param x A \code{\link[TFBSTools]{PFMatrix}} object
#' @param xmax A numeric scalar with the maximal width for the logo
#'     (in base-pairs). A value of \code{NULL} will scale the logo to the
#'     full width of the viewport.
#' @param ymax A numeric scalar with the maximal height for the logo (in bits)
#'     A value of \code{NULL} will scale the logo to the full height of the
#'     viewport.
#' @param xjust A character scalar specifying the horizontal adjustment of the
#'     sequence log withint the viewport; one of \code{"left"}, \code{"center"}
#'     or \code{"right"}.
#' @param ic.scale whether scale to bits or probability, default TRUE.
#'
#' @import grid
#'
#' @return A polygon grob.
seqLogoGrob2 <- function(x, xmax = NULL, ymax = 2.0, ic.scale = TRUE,
                         xjust = c("left", "center", 'right')) {
  stopifnot(is(x, "PFMatrix"))
  stopifnot(is.null(xmax) ||
              (is.numeric(xmax) && length(xmax) == 1L && xmax > 0))
  stopifnot(is.null(ymax) ||
              (is.numeric(ymax) && length(ymax) == 1L && ymax > 0))
  xjust <- match.arg(xjust)

  xm <- TFBSTools::Matrix(x)
  xm <- sweep(xm, MARGIN = 2, colSums(xm), "/")
  xm[is.nan(xm)] <- 0.25

  chars <- c("A", "C", "G", "T")
  letters <- list(x = NULL, y = NULL, id = NULL, fill = NULL)
  npos <- ncol(xm)
  facs <- pfm2ic(xm)
  wt <- 1
  x.pos <- 0

  if (ic.scale) {
    facs <- pfm2ic(xm)
  } else {
    facs <- rep(1, npos)
  }

  for (j in seq_len(npos)) {
    hts <- 0.95 * xm[, j] * facs[j]
    letterOrder <- order(hts)
    y.pos <- 0
    for (i in seq_len(4)) {
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht > 0)
        letters <- addLetter(letters, letter, x.pos, y.pos, ht, wt)
      y.pos <- y.pos + ht + 0.01
    }
    x.pos <- x.pos + wt
  }

  if (is.null(xmax)) {
    xmax <- max(letters$x) # use full width of viewport
  }
  if (is.null(ymax)) {
    ymax <- max(letters$y) # use full height of viewport
  }

  xoff <- switch(xjust,
                 left = 0,
                 center = (xmax - max(letters$x)) / 2,
                 right = (xmax - max(letters$x)))

  x <- unit((letters$x + xoff) / xmax, "npc")
  y <- unit(letters$y / ymax, "npc")

  grid::polygonGrob(x = x, y = y, id = letters$id,
                    name = as.character(ncol(xm)),
                    gp = grid::gpar(fill = letters$fill, col = "transparent"))
}
