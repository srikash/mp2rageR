#' @title mp2rage_b1_sa2rage_lookuptable
#'
#' @description
#' \code{mp2rage_b1_sa2rage_lookuptable} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{MATLAB code} made publicly available by José P. Marques.
#'
#' For methodological details, see \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069294}{Sa2RAGE correction}.
#'
#' Please cite the linked papers if you used these methods in your work.
#'
#' @author Sriranga Kashyap
#'
#' @param mprage_tr (required) \emph{TR} of the Sa2RAGE in s
#' @param flash_tr (required) \emph{TR} of the FLASH readout in s
#' @param inv_times_a_b (required) \eqn{c(TI1,TI2)}
#' @param flip_angle_a_b_deg (required) \emph{c(\eqn{\alpha}1,\eqn{\alpha}2)}
#' @param num_z_slices (required) calculate as: \eqn{base_resolution * c(phase_partial_fourier-0.5,0.5)/phase_ipat+c(ref_lines/2,ref_lines/2) * (1-1/phase_ipat)}
#' @param sequence_type (has default) \emph{"normal"}, set to \emph{NULL} for water excitation
#' @param t1_vector (has default) \emph{1.5}, provide input vector as necessary
#' @param b0 (has default) \emph{7.0} can be \emph{3.0}
#' @param m0 (has default) \emph{1.0}
#' @param inversion_efficiency (has default) \emph{0.96} is the efficiency of the Siemens MP2RAGE inversion pulse
#' @param n_images (has default) \emph{2}, change if required
#'
#' @export
#'
#' @return List of numerical arrays (in order): intensity, t1_vector, intensity_before_comb
#'
#' @importFrom pracma and
#'
#' @examples
#'
#' sa2rage_params <- list(
#'                    mprage_tr = 2.4,
#'                    flash_tr = 2.2e-3,
#'                    inv_times_a_b = c(58e-3,1800e-3),
#'                    flip_angle_a_b_deg = c(4,10),
#'                    num_z_slices = c(22,38)
#'                    )
#' \dontrun{
#' out_list <- mp2rage_b1_sa2rage_lookuptable(
#'                    sa2rage_params$mprage_tr,
#'                    sa2rage_params$flash_tr,
#'                    sa2rage_params$inv_times_a_b,
#'                    sa2rage_params$flip_angle_a_b_deg,
#'                    sa2rage_params$num_z_slices,
#'                    sequence_type = "normal",
#'                    t1_vector = 1.5,
#'                    b0 = 7,
#'                    m0 = 1,
#'                    inversion_efficiency = 0.96,
#'                    n_images = 2
#'                    )}

mp2rage_b1_sa2rage_lookuptable <-
  function(mprage_tr,
           flash_tr,
           inv_times_a_b,
           flip_angle_a_b_deg,
           num_z_slices,
           sequence_type = "normal",
           t1_vector = 1.5,
           b0 = 7,
           m0 = 1,
           inversion_efficiency = 0.96,
           n_images = 2) {
    inv_times_a <- inv_times_a_b[1]
    inv_times_b <- inv_times_a_b[2]

    flip_angle_a <- flip_angle_a_b_deg[1]
    flip_angle_b <- flip_angle_a_b_deg[2]

    b1_vector <- seq(from = 0.005, to = 2.5, by = 0.005)

    if (length(num_z_slices) == 2) {
      num_z_bef <- num_z_slices[1]
      num_z_aft <- num_z_slices[2]
      num_z_slices_2 <- num_z_slices
      num_z_slices <- sum(num_z_slices)
    } else if (length(num_z_slices) == 1) {
      num_z_bef <- num_z_slices / 2
      num_z_aft <- num_z_slices / 2
      num_z_slices_2 <- num_z_slices
    }

    signal <- matrix(data = NA,
                     nrow = length(b1_vector),
                     ncol = 2)

    for (j in seq_along(b1_vector)) {
      if (and(and(
        diff(inv_times_a_b) >= (num_z_slices * flash_tr),
        inv_times_a >= (num_z_bef * flash_tr)
      ),
      (inv_times_b <= (mprage_tr - num_z_aft * flash_tr))))
      {

        signal[j, 1:2] <- mp2rage_bloch_func(
          mprage_tr,
          flash_tr,
          inv_times_a_b,
          b1_vector[j] * c(flip_angle_a, flip_angle_b),
          num_z_slices_2,
          sequence_type,
          t1_vector,
          b0,
          m0,
          -cos(b1_vector[j] * pi / 2),
          n_images
        )
      }
      else {
        signal[j, 1:2] <- 0
      }
    }

    # Look up table
    intensity_out <- (Re(signal[, 1]) / Re(signal[, 2]))

    return(list(b1_vector=b1_vector, intensity=intensity_out, signal=signal))
  }
