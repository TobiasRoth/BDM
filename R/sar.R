#' Species area relationship
#'
#' The function can be used to standardize measures of species richness that
#' where sampled from sites that differ in size. Note that all measures with
#' zero species will be removed without furhter warning.
#'
#' @param SR Vector with measured number of species (species richness).
#' @param A Vector with same length as SR that contains the size of area within
#'   which SR were sampled.
#' @param stand Sample area for which species richness should be predicted.
#' @param minA Remove all samples that with sampler area smaler than 'minA'.
#' @return Predicted species richness for a sample area with stand area.
#' @export
#'
#' @examples
#' sar(SR = rpois(20, 60), A = rnorm(20 ,100, 25), stand = 200, minA = 20)
sar <- function(SR, A, stand, minA) {
  d <- tibble(SR = SR, A = A) %>% 
    filter(!is.na(SR) & !is.na(A)) %>% 
    filter(A >= minA) %>% 
    filter(SR > 0)
  koef <- coef(lm(log(SR) ~ log(A), data = d))
  SR_stand <- (exp(koef[1]) * stand^koef[2]) + (d$SR - (exp(koef[1]) * d$A^koef[2]))
  as.integer(round(SR_stand, 0))
}

