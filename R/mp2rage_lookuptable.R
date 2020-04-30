#' @title mp2rage_lookuptable
#'
#' @description
#' \code{mp2rage_lookuptable} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{MATLAB code} made publicly available by José P. Marques.
#'
#' For methodological details, see \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE}.
#'
#' Please cite the aforementioned paper if you used this method in your work.
#'
#' @author Sriranga Kashyap
#'
#' @param mprage_tr (required) \emph{TR} of the MP2RAGE in s
#' @param flash_tr (required) \emph{TR} of the FLASH readout in s
#' @param inv_times_a_b (required) \eqn{c(TI1,TI2)}
#' @param flip_angle_a_b_deg (required) \emph{c(\eqn{\alpha}1,\eqn{\alpha}2)}
#' @param num_z_slices (required) \eqn{slices_per_slab * c(slice_partial_fourier-0.5,0.5)}
#' @param sequence_type (has default) \emph{"normal"}, set to \emph{NULL} for water excitation
#' @param t1_vector (has default) \emph{NULL}, provide input vector as necessary
#' @param b0 (has default) \emph{7.0} can be \emph{3.0}
#' @param m0 (has default) \emph{1.0}
#' @param inversion_efficiency (has default) \emph{0.96} is the efficiency of the Siemens MP2RAGE inversion pulse
#' @param n_images (has default) \emph{2}, change if required
#' @param all_data (has default) \emph{0}
#'
#' @export
#' @return List of numerical arrays (in order): intensity, t1_vector, intensity_before_comb
#' @importFrom pracma isempty and
#' @examples
#' \dontrun{
#' mp2rage_params <- list(
#'                    mprage_tr = 5.0,
#'                    flash_tr = 6.9e-3,
#'                    inv_times_a_b = c(900e-3,2750e-3),
#'                    flip_angle_a_b_deg = c(5,3),
#'                    num_z_slices = c(120,120)
#'                    )
#'
#' out_list <- mp2rage_lookuptable(
#'                    mp2rage_params$mprage_tr,
#'                    mp2rage_params$flash_tr,
#'                    mp2rage_params$inv_times_a_b,
#'                    mp2rage_params$flip_angle_a_b_deg,
#'                    mp2rage_params$num_z_slices,
#'                    sequence_type = "normal",
#'                    t1_vector = NULL,
#'                    b0 = 7,
#'                    m0 = 1,
#'                    inversion_efficiency = 0.96,
#'                    n_images = 2,
#'                    all_data = 0
#'                    )}

mp2rage_lookuptable <-
  function(mprage_tr,
           flash_tr,
           inv_times_a_b,
           flip_angle_a_b_deg,
           num_z_slices,
           sequence_type = "normal",
           t1_vector = NULL,
           b0 = 7,
           m0 = 1,
           inversion_efficiency = 0.96,
           n_images = 2,
           all_data = 0) {
    inv_times_a <- inv_times_a_b[1]
    inv_times_b <- inv_times_a_b[2]

    flip_angle_a <- flip_angle_a_b_deg[1]
    flip_angle_b <- flip_angle_a_b_deg[2]

    if (is.null(t1_vector)) {
      t1_vector <- seq(from = 0.05, to = 5.0, by = 0.05)
    }

    if (length(num_z_slices) == 2) {
      num_z_bef <- num_z_slices[1]
      num_z_aft <- num_z_slices[2]
      num_z_slices_2 <- sum(num_z_slices)
    } else if (length(num_z_slices) == 1) {
      num_z_bef <- num_z_slices / 2
      num_z_aft <- num_z_slices / 2
      num_z_slices_2 <- num_z_slices
    }

    signal <- matrix(data = NA,
                     nrow = length(t1_vector),
                     ncol = 2)

    for (j in seq_along(t1_vector)) {
      if (and(and(
        diff(inv_times_a_b) >=
        (num_z_bef * flash_tr +
         num_z_aft * flash_tr),
        (inv_times_a >=
         num_z_bef * flash_tr)
      ), (inv_times_b <=
          (mprage_tr - num_z_aft * flash_tr)))) {
        signal[j, 1:2] <- mp2rage_bloch_func(
          mprage_tr,
          flash_tr,
          inv_times_a_b,
          c(flip_angle_a, flip_angle_b),
          num_z_slices_2,
          sequence_type,
          t1_vector[j],
          b0,
          m0,
          inversion_efficiency = 1,
          n_images
          )
      }
      else {
        signal[j, 1:2] <- 0
      }
    }

    # Look up table
    intensity <- (Re(signal[, 1] * Conj(signal[, 2])) /
                    (abs(signal[, 1]) ^ 2 + abs(signal[, 2]) ^ 2))

    if (all_data == 0) {
      min_index=which.max(intensity)
      max_index=which.min(intensity)
      intensity <- intensity[min_index:max_index]

      t1_vector <-
        t1_vector[min_index:max_index]

      intensity_before_comb <-
        signal[min_index:max_index, ]

      intensity[c(1, length(intensity))] <- c(0.5,-0.5)

      return(
        list(
          intensity = intensity,
          t1_vector = t1_vector,
          intensity_before_comb = intensity_before_comb
        )
      )

    } else {
      intensity_before_comb <- signal

      return(
        list(
          intensity = intensity,
          t1_vector = t1_vector,
          intensity_before_comb = intensity_before_comb
        )
      )
    }

  }
