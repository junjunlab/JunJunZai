globalVariables(c("Occupancy", "TF", "group", "name", "pos"))

makePFM <- getFromNamespace("makePFM","ggseqlogo")

# ==============================================================================
# footprint plot
# ==============================================================================

#' plot footprint line plot from rgt-hint output
#'
#' @param linePath the rgt-hint output directory, default NULL.
#' @param TFname the transcription factor name to be plotted, default NULL.
#'
#' @return a ggplot
#' @export
#' @importFrom stringr str_detect str_replace
#' @importFrom tidyr pivot_longer
#' @import ggplot2
footprintPlot <- function(linePath = NULL,TFname = NULL){
  fp_file <- list.files(path = linePath,pattern = "txt",full.names = T)
  fp_file_name <- list.files(path = linePath,pattern = "txt")
  fp_file_name <- sapply(strsplit(fp_file_name,split = ".txt"), "[",1)

  # loop load footprinting
  # x = 719
  lapply(seq_along(fp_file),function(x){
    # check names

    var <- stringr::str_detect(fp_file_name[x],pattern = "var")

    if(var){
      tmp_name <- stringr::str_replace(fp_file_name[x],pattern = "_",replacement = "(")
      tmp_name <- paste0(tmp_name,")")
    }else{
      tmp_name <- fp_file_name[x]
    }

    # load data
    ft_tmp <- read.delim(fp_file[x]) %>%
      mutate(pos = -100:99,TF = tmp_name)
    ft <- ft_tmp %>%
      tidyr::pivot_longer(cols = colnames(ft_tmp)[1:(ncol(ft_tmp) - 2)],
                          names_to = "group",
                          values_to = "Occupancy")

  }) %>% do.call("rbind",.) -> ft_df

  # footprinting line plot
  footprint <-
    ggplot(ft_df %>%
             dplyr::filter(TF == TFname),
           aes(x = pos,y = Occupancy,color = group)) +
    geom_line() +
    theme_bw() +
    theme(panel.grid = element_blank(),
          strip.background = element_rect(fill = "grey90"),
          legend.background = element_blank(),
          strip.text = element_text(face = "bold.italic",size = rel(1)),
          axis.text = element_text(colour = "black")) +
    facet_wrap(~TF)

  footprint
}


# ==============================================================================
# get TF sequence
# ==============================================================================


#' extract sequence from bed
#'
#' @param mbpsPath the mbps file from rgt-motifanalysis output directory, default NULL.
#' @param TFname the transcription factor name to be plotted, default NULL.
#' @param genome the BSgenome object for genome sequence, default NULL.
#'
#' @return sequence
#' @export
#' @importFrom rtracklayer import.bed
#' @importFrom IRanges resize
#' @importFrom ChIPpeakAnno getAllPeakSequence
getSeqFrombed <- function(mbpsPath = NULL,TFname = NULL,genome = NULL){
  mbps <- list.files(path = mbpsPath,pattern = "mpbs.bed",full.names = T)
  lapply(mbps,function(x){
    rtracklayer::import.bed(x) %>% subset(.,name == TFname)
  }) %>% do.call("c",.) -> TFtarget

  TFtarget_ed <- IRanges::resize(TFtarget,width = 200,fix = "center")

  # subtract sequence for TF
  mtf_seq <- ChIPpeakAnno::getAllPeakSequence(myPeakList = TFtarget_ed,
                                              upstream = 0,downstream = 0,
                                              genome = genome)

  return(mtf_seq$sequence)
}



# ==============================================================================
# sequence to pwmMatrix
# ==============================================================================


#' create a pwmMatrix
#'
#' @param seqs sequence to make pwmMatrix, default NULL.
#' @param shift a number of shift to zoom the motif, default c(-10,10).
#'
#' @return A list pwmMatrix
#' @export
#' @importFrom TFBSTools PWMatrix
seq2pwmMatrix <- function(seqs = NULL,shift = c(-10,10)){
  pfm <- makePFM(seqs = seqs,seq_type = "dna")

  pwm_mat <- TFBSTools::PWMatrix(ID = "Unknown", name = "Unknown", matrixClass = "Unknown",
                                 strand = "+",
                                 bg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                 tags = list(),
                                 profileMatrix = matrix(pfm[c(1,4,3,2),],
                                                        byrow = FALSE, nrow = 4,
                                                        dimnames=list(c("A", "C", "G", "T"))),
                                 pseudocounts = numeric())

  # get zoom range
  mid <- ncol(pfm)/2
  zoom_rg <- c((mid + shift[1]):(mid + shift[2]))

  pwm_mat_sub <- TFBSTools::PWMatrix(ID = "Unknown", name = "Unknown", matrixClass = "Unknown",
                                     strand = "+",
                                     bg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                     tags = list(),
                                     profileMatrix = matrix(pfm[c(1,4,3,2),zoom_rg],
                                                            byrow = FALSE, nrow = 4,
                                                            dimnames=list(c("A", "C", "G", "T"))),
                                     pseudocounts = numeric())

  return(list(pwm_mat = pwm_mat,
              pwm_mat_sub = pwm_mat_sub,
              shift = shift))
}
