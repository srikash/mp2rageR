#' @title mp2rage_lookuptable
#'
#' @description
#' \code{mp2rage_lookuptable} is part of the \pkg{mp2rageR} package.
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
#' @importFrom pracma isempty and
#' @examples
#' intensity_t1_vector_intensity_bef_comb <- mp2rage_lookuptable(mprage_tr,
#' inv_times_a_b,
#' flip_angle_a_b_deg,
#' num_z_slices,
#' flash_tr,
#' sequence_type,
#' n_images = 2,
#' b0 = 7.0,
#' m0 = 1.0,
#' inversion_efficiency = 0.96,
#' all_data = 0,
#' t1_vector = NULL)

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
