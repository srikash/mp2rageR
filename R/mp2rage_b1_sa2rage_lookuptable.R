#' @title mp2rage_b1_sa2rage_lookuptable
#'
#' @description
#' \code{mp2rage_b1_sa2rage_lookuptable} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of MaTLab code made publicly
#' available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#'
#' For methodological details, see \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE paper}.
#'
#' Please cite the aforementioned paper if you used this method in your work.
#'
#' @author Sriranga Kashyap
#'
#' @export
#' @return List of numerical arrays (in order): intensity, t1_vector, intensity_before_comb
#' @importFrom pracma and
#' @examples
#' list_of_b1_vector_intensity_signal <- mp2rage_b1_sa2rage_lookuptable(mprage_tr,
#'                                                                      flash_tr,
#'                                                                      inv_times_a_b,
#'                                                                      flip_angle_a_b_deg,
#'                                                                      num_z_slices,
#'                                                                      sequence_type = "normal",
#'                                                                      t1_vector = 1.5,
#'                                                                      b0 = 7,
#'                                                                      m0 = 1,
#'                                                                      inversion_efficiency = 0.96,
#'                                                                      n_images = 2)
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
