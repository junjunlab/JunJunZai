#' create homerResult object class
#'
#' @slot known res data.frame. known meta info.
#' @slot known motif PWM matrix list. motif PWM matrix.
#' @slot known motif PFM matrix list. motif PFM matrix.
#' @slot novo res data.frame. novo meta info.
#' @slot novo motif PWM matrix list. motif PWM matrix.
#' @slot novo motif PFM matrix list. motif PFM matrix.
#'
#' @return homerResult object
#' @export
homerResult <- setClass("homerResult",
                        slots = list("known res" = "data.frame",
                                     "known motif PWM matrix" = "list",
                                     "known motif PFM matrix" = "list",
                                     "novo res" = "data.frame",
                                     "novo motif PWM matrix" = "list",
                                     "novo motif PFM matrix" = "list"))
