globalVariables(c('type', 'ic.scale', 'ymax', 'clustForMotif', 'cluster_rows', 'logoSpace',
                  'logoWidth', 'rightAnnoWidth', 'barWidth', 'col1', 'col2', '...'))

#' Default method
#'
#' @param object homerResult.
#'
#' @import utils
#'
#' @export
setMethod("show",
          signature(object = "homerResult"),
          function(object){
            cat("It's a homerResult. Here showing the known motif content.\n")

            print(head(object@`known res`,3))
          }
)


# ==============================================================================
# gglogo function for plot
# ==============================================================================
#' motif logo plot
#'
#' @param object homerResult.
#' @param type the homer output type results to choose, "known" or "novo".
#' @param ... other args can be passed with ggseqLogo.
#'
#' @export
setGeneric("gglogo",function(object,type = c("known","novo"),...) standardGeneric("gglogo"))



#' method for gglogo
#'
#' @param object homerResult.
#'
#' @import ggseqlogo
#'
#' @return ggplot object
#' @export
setMethod("gglogo",
          signature(object = "homerResult"),
          function(object,type = c("known","novo"),...){
            type <- match.arg(type,c("known","novo"))

            if(type == "known"){
              mtx_list <- object@`known motif PWM matrix`
            }else{
              mtx_list <- object@`novo motif PWM matrix`
            }

            # extract pwm list
            lapply(seq_along(mtx_list),function(x){
              mtx <- mtx_list[[x]]@pwm
            }) -> mtxlist

            lapply(seq_along(mtx_list),function(x){
              tf <- names(mtx_list[[x]])
            }) %>% unlist() -> tf

            # make unique names
            uni_tf <- make.unique(tf)

            # assign names for list
            names(mtxlist) <- uni_tf

            # plot
            ggseqlogo(mtxlist,...)
          }
)



# ==============================================================================
# plotMotifHeatmap function for plot
# ==============================================================================
#' heatmap plot for homerResult
#'
#' @param object homerResult.
#' @param type the homer output type results to choose, "known" or "novo".
#' @param ic.scale whether scale to bits or probability, default TRUE.
#' @param ymax seqlogo height, default 2.
#' @param clustForMotif whether cluster rows by motif similarity, default FALSE.
#' @param cluster_rows whether cluster rows, default TRUE.
#' @param logoSpace the space between motif logo and heatmap, default 0.5.
#' @param logoWidth motif logo width, default 0.8.
#' @param rightAnnoWidth right annotation width, default 0.5.
#' @param barWidth the match score for de novo results annotation width, default 3.
#' @param col1 the color for negLogPvalue heatmap, for example: col1 = colorRamp2(c(-2, 0, 2), c("green", "white", "red")).
#' @param col2 the color for log2Enrichment heatmap.
#' @param ... other args can be passed with ComplexHeatmap.
#'
#' @return Complexheatmap list
#' @export
setGeneric("plotMotifHeatmap",function(object,type = c("known","novo"),
                                       ic.scale = TRUE,ymax = 2,
                                       clustForMotif = FALSE,
                                       cluster_rows = TRUE,
                                       logoSpace = 0.5,
                                       logoWidth = 0.8,
                                       rightAnnoWidth = 0.5,
                                       barWidth = 3,
                                       col1 = NULL,col2 = NULL,
                                       ...) standardGeneric("plotMotifHeatmap"))





#' plotMotifHeatmap method
#'
#' @param object homerResult.
#'
#' @import ComplexHeatmap
#' @importFrom dplyr select filter mutate if_else
#' @import methods
#' @importFrom stats as.dist hclust
#' @importFrom TFBSTools PWMatrix PFMatrixList
#'
#' @export
setMethod("plotMotifHeatmap",
          signature(object = "homerResult"),
          function(object,type = c("known","novo"),
                   ic.scale = TRUE,ymax = 2,
                   clustForMotif = FALSE,
                   cluster_rows = TRUE,
                   logoSpace = 0.5,
                   logoWidth = 0.8,
                   rightAnnoWidth = 0.5,
                   barWidth = 3,
                   col1 = NULL,col2 = NULL,...){
            type <- match.arg(type,c("known","novo"))

            # check type and extract data
            if(type == "known"){
              meta <- object@`known res`
              mtx_list <- object@`known motif PWM matrix`
              pfm_List <-  object@`known motif PFM matrix`
            }else{
              meta <- object@`novo res`
              mtx_list <- object@`novo motif PWM matrix`
              pfm_List <-  object@`novo motif PFM matrix`
            }

            # heatmap mat
            htmat <- meta %>% select(neglogPvalue)
            enrmat <- meta %>% select(log2Ratio) %>%
              mutate(log2Ratio = if_else(is.infinite(log2Ratio),0,log2Ratio))

            # top annotation
            top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "grey85"),
                                                                labels = c("neglogPvalue"),
                                                                labels_gp = gpar(col = "black",
                                                                                 fontsize = 10,
                                                                                 fontface = "bold")))

            # row annotation
            if(type == "known"){
              row_ha = rowAnnotation(tgPer = meta$TargetPercent,
                                     bgPer =  meta$BackgroundPercent,
                                     show_annotation_name = FALSE,
                                     border = TRUE,
                                     simple_anno_size = unit(rightAnnoWidth , "cm"))
            }else{
              row_ha = rowAnnotation(tgPer = meta$TargetPercent,
                                     bgPer =  meta$BackgroundPercent,
                                     matchScore = anno_barplot(meta$matchScore,
                                                               add_numbers = T,
                                                               width = unit(barWidth,"cm")),
                                     show_annotation_name = FALSE,
                                     border = TRUE,
                                     simple_anno_size = unit(rightAnnoWidth , "cm"))
            }

            # motif logo annotation
            # x = 1
            grobL <- lapply(seq_along(mtx_list), function(x){
              pwm_ele <- mtx_list[[x]]
              pwm_mat <- PWMatrix(ID = "Unknown", name = "Unknown", matrixClass = "Unknown",
                                  strand = "+",
                                  bg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                                  tags = list(),
                                  profileMatrix = as.matrix(pwm_ele@pwm),
                                  pseudocounts = numeric())

              logogrob <- seqLogoGrob2(x = pwm_mat,ic.scale = ic.scale, ymax = ymax)

              return(logogrob)
            })


            # logo annotation
            hmSeqlogo <- HeatmapAnnotation(logo = annoSeqlogo(grobL = grobL, which = "row",
                                                              space = unit(logoSpace, "mm"),
                                                              width = unit(logoWidth, "inch")),
                                           show_legend = FALSE,
                                           show_annotation_name = FALSE,
                                           which = "row")

            # whether clust for motif similarity
            if(clustForMotif == TRUE){
              SimMatSel <- monaLisa::motifSimilarity(x = pfm_List[[1]])

              cr <- hclust(as.dist(1 - SimMatSel), method = "average")
            }else{
              cr <- cluster_rows
            }


            # final plot
            htmain <-
              Heatmap(matrix = htmat,
                      name = "negLogP",
                      col = col1,
                      top_annotation = top_annotation,
                      show_column_names = F,
                      row_names_side = "left",
                      show_row_dend = T,
                      cluster_rows = cr,
                      left_annotation = hmSeqlogo,
                      right_annotation = row_ha,
                      border = T,
                      ...
              )

            # ==================================================================
            # add log2enrichment
            # ==================================================================
            enr_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "grey85"),
                                                                labels = c("log2 enrichment"),
                                                                labels_gp = gpar(col = "black",
                                                                                 fontsize = 10,
                                                                                 fontface = "bold")))


            htsub <-
              Heatmap(matrix = enrmat,
                      name = "log2Ratio",
                      col = col2,
                      cluster_rows = F,
                      show_column_names = F,
                      show_row_names = F,
                      top_annotation = enr_annotation,
                      border = T,
                      ...)

            # output
            htlist <- htmain + htsub
            ComplexHeatmap::draw(object = htlist)

          }
)
