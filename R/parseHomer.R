globalVariables(c(".", "BackgroundPercent", "TargetPercent", "if_else",
                  "log2Ratio", "neglogPvalue"))

# ==============================================================================
# basic parse function
# ==============================================================================
#' parse homer output into meta information
#'
#' @param homerDir the path of homer output data.
#' @param motifIndex the numbers of motif index.
#' @param type the homer output type results to choose, "known" or "novo".
#' @param pattern The pattern to extract de novo match score, default NULL.
#'
#' @return data.frame
#' @export
parseHomer <- function(homerDir = NULL,
                       motifIndex = NULL,
                       type = c("known","novo"),
                       pattern = NULL){
  type <- match.arg(type,c("known","novo"))

  if(type == "known"){
    motif.path <- paste0(homerDir, "knownResults/","known",motifIndex,".motif")
  }else{
    motif.path <- paste0(homerDir, "homerResults/","motif",motifIndex,".motif")
  }

  # x = 1
  lapply(seq_along(motif.path),function(x){
    # summary info
    row1 <- read.table(motif.path[x],nrows = 1,sep = "#")
    info <- unlist(strsplit(row1$V1,split = '\t'))

    # prepare meta info
    if(type == "known"){
      consensus <- sapply(strsplit(info[1],split = ">"),"[",2)
      tf <- sapply(strsplit(info[2],split = "\\("),"[",1)
      experiment <- sapply(strsplit(info[2],split = "\\/"),"[",2)
      matchs <- 0
    }else if(type == "novo"){
      consensus <- sapply(strsplit(info[1],split = ">"),"[",2)
      tf_tmp <- sapply(strsplit(info[2],split = "\\("),"[",1)
      tf <- sapply(strsplit(tf_tmp,split = "\\:|\\/"),"[",2)
      experiment <- sapply(strsplit(info[2],split = "\\/"),"[",2)
      matchs_tmp <- sapply(strsplit(info[2],split = "\\/"),"[",3)

      # check value
      if(is.na(matchs_tmp)){
        matchs_tmp <- sapply(strsplit(info[2],split = "\\/"),"[",2)
      }

      matchs <- sapply(strsplit(matchs_tmp,split = "\\(|\\)"),"[",2) %>% as.numeric()

      if(is.na(matchs)){
        matchs <- sapply(strsplit(info[2],split = pattern[1]),"[",as.numeric(pattern[2])) %>%
          as.numeric()
      }
    }


    pval <- sapply(strsplit(info[6],split = "P:"),"[",2)
    neglogpval <- -as.numeric(info[4])
    qval <- as.numeric(info[5])
    targetper <- sapply(strsplit(info[6],split = "\\(|\\)"),"[",2) %>%
      strsplit(.,split = "%") %>% unlist() %>% as.numeric()
    backgroundper <- sapply(strsplit(info[6],split = "\\(|\\)"),"[",4) %>%
      strsplit(.,split = "%") %>% unlist() %>% as.numeric()


    # merge into dataframe
    info_df <- data.frame(consensus = consensus,
                          TF = tf,
                          experiment = experiment,
                          matchScore = matchs,
                          Pvalue = pval,
                          neglogPvalue = neglogpval,
                          Qvalue = qval,
                          TargetPercent = targetper,
                          BackgroundPercent = backgroundper) %>%
      dplyr::mutate(ratio = TargetPercent/BackgroundPercent,
                    log2Ratio = log2(TargetPercent/BackgroundPercent))

    rownames(info_df) <- tf

    return(info_df)
  }) %>% do.call("rbind",.) -> homerRes

}

# ==============================================================================
# basic parse motif matrix function
# ==============================================================================
#' extract pwm matrix from motif
#'
#' @param homerDir the path of homer output data.
#' @param motifIndex the numbers of motif index.
#' @param type the homer output type results to choose, "known" or "novo".
#'
#' @return pwmList
#' @export
parseHomerMotif <- function(homerDir = NULL,
                            motifIndex = NULL,
                            type = c("known","novo")){
  type <- match.arg(type,c("known","novo"))

  if(type == "known"){
    motif.path <- paste0(homerDir, "knownResults/","known",motifIndex,".motif")
  }else{
    motif.path <- paste0(homerDir, "homerResults/","motif",motifIndex,".motif")
  }

  # load ppm matrix
  lapply(seq_along(motif.path),function(x){
    # summary info
    row1 <- read.table(motif.path[x],nrows = 1,sep = '_')
    info <- unlist(strsplit(row1$V1,split = '\t'))

    if(type == "known"){
      tf <- sapply(strsplit(info[2],split = "\\("),"[",1)
    }else{
      tf_tmp <- sapply(strsplit(info[2],split = "\\("),"[",1)
      tf <- sapply(strsplit(tf_tmp,split = "\\:|\\/"),"[",2)
    }

    # load ppm matrix
    mat <- read.table(motif.path[x],skip = 1) %>% t() %>% as.matrix()

    # matrix2pwm
    pwm <- seqLogo::makePWM(mat)

    # assign name
    attr(pwm,"names") <- tf

    return(pwm)
  }) -> pwmList

  return(pwmList)
}


# ==============================================================================
# basic parse motif PFM matrix function
# ==============================================================================
#' extract pfm matrix from motif
#'
#' @param homerDir the path of homer output data.
#' @param motifIndex the numbers of motif index.
#' @param type the homer output type results to choose, "known" or "novo".
#'
#' @return pfmList
#' @export
preparePFMmat <- function(homerDir = NULL,
                          motifIndex = NULL,
                          type = c("known","novo")){
  type <- match.arg(type,c("known","novo"))

  if(type == "known"){
    motif.path <- paste0(homerDir, "knownResults/","known",motifIndex,".motif")
  }else{
    motif.path <- paste0(homerDir, "homerResults/","motif",motifIndex,".motif")
  }

  # load ppm matrix
  lapply(seq_along(motif.path),function(x){
    # summary info
    row1 <- read.table(motif.path[x],nrows = 1,sep = '_')
    info <- unlist(strsplit(row1$V1,split = '\t'))

    if(type == "known"){
      tf <- sapply(strsplit(info[2],split = "\\("),"[",1)
    }else{
      tf_tmp <- sapply(strsplit(info[2],split = "\\("),"[",1)
      tf <- sapply(strsplit(tf_tmp,split = "\\:|\\/"),"[",2)
    }

    # load pfm matrix
    pfm <- monaLisa::homerToPFMatrixList(filename = motif.path[x])

    return(pfm)
  }) -> pfmList

  # PFMatrix into PFMatrixList
  # x = 1
  lapply(seq_along(pfmList), function(x){
    pfm_ele <- pfmList[[x]]
    return(pfm_ele@listData[[1]])
  }) %>% do.call("PFMatrixList",.) -> pfm_List

  return(pfm_List)
}


# ==============================================================================
# combine pwm matrix and meta info
# ==============================================================================
#' parse homer output and assemble into homerResult object
#'
#' @param homerDir the path of homer output data.
#' @param motifIndex the numbers of motif index.
#' @param novo whether parse de novo motif data, default TRUE.
#' @param known whether parse de known motif data, default TRUE.
#' @param pattern The pattern to extract de novo match score, default NULL.
#'
#' @return homerResult object
#' @export
loadHomerRes <- function(homerDir = NULL,
                         motifIndex = NULL,
                         novo = TRUE,
                         known = TRUE,
                         pattern = NULL){
  # ============================================================================
  # deal with known motif
  # ============================================================================
  # extract known motif
  if(known == TRUE){
    known <- parseHomer(homerDir = homerDir,
                        motifIndex = motifIndex,
                        type = "known")

    known_pwm <- parseHomerMotif(homerDir = homerDir,
                                 motifIndex = motifIndex,
                                 type = "known")

    known_pfm <- preparePFMmat(homerDir = homerDir,
                               motifIndex = motifIndex,
                               type = "known")
  }else{
    known <- data.frame()
    known_pwm <- list()
    known_pfm <- list()
  }

  # ============================================================================
  # deal with novo motif
  # ============================================================================
  # extract novo motif
  if(novo == TRUE){
    novo <- parseHomer(homerDir = homerDir,
                       motifIndex = motifIndex,
                       type = "novo",
                       pattern = pattern)

    novo_pwm <- parseHomerMotif(homerDir = homerDir,
                                motifIndex = motifIndex,
                                type = "novo")

    novo_pfm <- preparePFMmat(homerDir = homerDir,
                              motifIndex = motifIndex,
                              type = "novo")
  }else{
    novo <- data.frame()
    novo_pwm <- list()
    novo_pfm <- list()
  }


  # ============================================================================
  # create homerResult object
  # ============================================================================
  res <- methods::new("homerResult",
                      "known res" = known,
                      "known motif PWM matrix" = known_pwm,
                      "known motif PFM matrix" = list(known_pfm),
                      "novo res" = novo,
                      "novo motif PWM matrix" = novo_pwm,
                      "novo motif PFM matrix" = list(novo_pfm))

  return(res)
}
