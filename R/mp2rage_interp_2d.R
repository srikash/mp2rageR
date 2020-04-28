#' @title mp2rage_interp_2d
#'
#' @description
#' \code{mp2rage_interp_2d} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of MATLAB code made publicly available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#'
#' For methodological details, see \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.1910150117}{MPRAGE paper} and \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE paper}.
#'
#' @author Sriranga Kashyap
#'
#' @export
#' @return Data as a numerical array
#' @importFrom akima bicubic
#' @examples
#' c_int <- mp2rage_interp_2d(a, b, c, a_int, b_int)
mp2rage_interp_2d <- function(a, b, c, a_int, b_int) {
  c_temp <- bicubic(
    x = a,
    y = b,
    z = t(c),
    x0 = a_int,
    y0 = b_int
  )

  c_int <- array(data = c_temp$z, dim = dim(a_int))

  return(c_int)
}
