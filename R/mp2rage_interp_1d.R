#' @title mp2rage_interp_1d
#'
#' @description
#' \code{mp2rage_interp_1d} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of MATLAB code made publicly available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#'
#' For methodological details, see \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.1910150117}{MPRAGE paper} and \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE paper}.
#'
#' @author Sriranga Kashyap
#'
#' @export
#' @return Data as a numerical array
#' @importFrom signal interp1
#' @importFrom stats spline
#' @examples
#' bi <- mp2rage_interp_1d(a, b, a_int, is_linear = TRUE)
mp2rage_interp_1d <- function(a, b, a_int, is_linear = TRUE) {
  if (isTRUE(is_linear)) {
    b_int <- interp1(
      x = a,
      y = c(b),
      xi = as.vector(a_int),
      method = "linear",
      extrap = TRUE
    )
    return(b_int)
  }
  else {
    b_temp <- spline(
      x = a,
      y = c(b),
      xout = as.vector(a_int),
      method = "fmm"
    )

    b_int <- b_temp$y

    return(b_int)
  }
}
