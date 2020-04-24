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
#' @importFrom pracma isempty
#' @examples
#' list_of_intensity_t1_vector_intensity_before_comb <- mp2rage_b1_sa2rage_lookuptable(mprage_tr,
#'                                                                          inv_times_a_b,
#'                                                                          flip_angle_a_b_deg,
#'                                                                          num_z_slices,
#'                                                                          flash_tr,
#'                                                                          sequence_type,
#'                                                                          ...)
mp2rage_b1_sa2rage_lookuptable <-
  function(mprage_tr,
           inv_times_a_b,
           flip_angle_a_b_deg,
           num_z_slices,
           flash_tr,
           sequence_type,
           n_images = 2,
           b0 = 7,
           m0 = 1,
           inversion_efficiency = 0.96,
           all_data = 0,
           t1_vector = NULL) {
    inv_times_a <- inv_times_a_b[1]
    inv_times_b <- inv_times_a_b[2]

    flip_angle_a <- flip_angle_a_b_deg[1]
    flip_angle_b <- flip_angle_a_b_deg[2]

    if (is.null(t1_vector)) {
      t1_vector <- seq(from = 0.05, to = 5.0, by = 0.05)

      flash_tr <- as.array(flash_tr)
    }

    if (length(flash_tr) == 1) {
      flash_tr <- rep(flash_tr, times = n_images)
      num_z_slices <- as.array(num_z_slices)
    }

    if (length(num_z_slices) == 2) {
      num_z_bef <- num_z_slices[1]
      num_z_aft <- num_z_slices[2]
      num_z_slices_2 <- sum(num_z_slices)
    } else if (length(num_z_slices) == 1) {
      num_z_bef <- num_z_aft <- num_z_slices / 2
      num_z_slices_2 <- num_z_slices
    }

    signal <- as.array(0, c(length(t1_vector), 2))

    for (j in seq_along(t1_vector)) {
      if (and(and(
        diff(inv_times_a_b) >=
        (num_z_bef * flash_tr[2] +
         num_z_aft * flash_tr[1]),
        (inv_times_a >=
         num_z_bef * flash_tr[1])
      ), (inv_times_b <=
          (mprage_tr - num_z_aft * flash_tr[2])))) {
        signal[j,] <- mprage_func(
          mprage_tr,
          inv_times_a_b,
          num_z_slices_2,
          flash_tr,
          c(flip_angle_a, flip_angle_b),
          sequence_type,
          t1_vector[j],
          n_images,
          b0,
          m0,
          inversion_efficiency
        )
      }
      else {
        signal[j,] <- 0
      }
    }

    # Look up table
    intensity <- drop(Re(signal[, 1] * Conj(signal[, 2])) /
                        (abs(signal[, 1]) ^ 2 + abs(signal[, 2]) ^ 2))

    if (all_data == 0) {
      intensity <- intensity[which.min(intensity):which.max(intensity)]

      t1_vector <-
        t1_vector[which.min(intensity):which.max(intensity)]

      intensity_before_comb <-
        signal[which.min(intensity):which.max(intensity)]

      intensity[c(1, tail(intensity, 1))] <- c(0.5,-0.5)

    } else {
      intensity = NULL

      t1_vector = NULL

      intensity_before_comb <- signal
    }
    return(list(intensity, t1_vector, intensity_before_comb))
  }
