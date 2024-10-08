#' @title mp2rage_bloch_func
#'
#' @description
#' \code{mp2rage_bloch_func} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{MATLAB code} made publicly available by José P. Marques.
#'
#' For methodological details, see \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.1910150117}{MPRAGE} and \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE}.
#'
#' Please cite the linked papers if you used these methods in your work.
#'
#' @author Sriranga Kashyap
#'
#' @param mprage_tr (required) \emph{TR} of the MP2RAGE in s
#' @param flash_tr (required) \emph{TR} of the FLASH readout in s
#' @param inv_times_a_b (required) \eqn{c(TI1,TI2)}
#' @param flip_angle_a_b_deg (required) \emph{c(\eqn{\alpha}1,\eqn{\alpha}2)}
#' @param num_z_slices (required) calculate as: \eqn{slices_per_slab * c(slice_partial_fourier-0.5,0.5)}
#' @param sequence_type (has default) \emph{"normal"}, set to \emph{NULL} for water excitation
#' @param t1_vector (has default) \emph{NULL}, provide input vector as necessary
#' @param b0 (has default) \emph{7.0} can be \emph{3.0}
#' @param m0 (has default) \emph{1.0}
#' @param inversion_efficiency (has default) \emph{0.96} is the efficiency of the Siemens MP2RAGE inversion pulse
#' @param n_images (has default) \emph{2}, change if required
#'
#' @export
#'
#' @return Data as a numerical array
#'
#' @importFrom pracma strcmp isempty zeros
#'
#' @examples
#' mp2rage_params <- list(
#'                    mprage_tr = 5.0,
#'                    flash_tr = 6.9e-3,
#'                    inv_times_a_b = c(900e-3,2750e-3),
#'                    flip_angle_a_b_deg = c(5,3),
#'                    num_z_slices = c(120,120)
#'                    )
#'
#' \dontrun{
#' signal <- mp2rage_bloch_func(
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
#'                    n_images = 2
#'                    )}

mp2rage_bloch_func <-
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
           n_images = 2)
  {
    # fat sat check
    if (isTRUE(sequence_type == "normal")) {
      normal_sequence_type <- TRUE
      water_excitation <- FALSE
    } else {
      normal_sequence_type <- FALSE
      water_excitation <- TRUE
    }

    num_z_slices <- as.array(num_z_slices)
    inv_times_a_b <- as.array(inv_times_a_b)
    flash_tr <- as.array(flash_tr)
    flip_angle_a_b_deg <- as.array(flip_angle_a_b_deg)

    fat_water_cs <- 3.3 #ppm
    gamma <- 42.576 #MHz/T
    pulse_space <- 1 / 2 / (fat_water_cs * b0 * gamma)

    #conversion from degrees to radians
    flip_angle_a_b_rad <-
      flip_angle_a_b_deg / 180 * pi

    if (n_images != length(flip_angle_a_b_rad)) {
      flip_angle_a_b_rad <- rep(flip_angle_a_b_rad, times = n_images)
    }
    # if (n_images != length(flash_tr)) {
    #   flash_tr <- rep(flash_tr, times = n_images)
    # }

    if (length(num_z_slices) == 2) {
      num_z_bef <- num_z_slices[1]
      num_z_aft <- num_z_slices[2]
      num_z_slices <- sum(num_z_slices)
    } else if (length(num_z_slices) == 1) {
      num_z_bef <- num_z_slices / 2
      num_z_aft <- num_z_slices / 2
    }
    ## calculating the relevant timing and associated values
    if (isTRUE(water_excitation)) {
      e_1 <- exp(-flash_tr / t1_vector)
      e_1a <- exp(-pulse_space / t1_vector)
      e_2a <-
        exp(-pulse_space / 0.06) #60 ms is an estimation of the T2*
      e_1b <- exp(-(flash_tr - pulse_space) / t1_vector)

      ta <- num_z_slices * flash_tr
      ta_bef <- num_z_bef * flash_tr
      ta_aft <- num_z_aft * flash_tr

      td[1] <- inv_times_a_b[1] - ta_bef
      td[n_images + 1] = mprage_tr - inv_times_a_b[n_images] - ta_aft

      e_td[1] <- exp(-td[1] / t1_vector)
      e_td[n_images + 1] = exp(-td[n_images + 1] / t1_vector)

      cos_alpha_e1 <- matrix()
      one_minus_e1 <- matrix()
      sin_alpha <- matrix()

      if (n_images > 1)
      {
        for (k in seq(2, n_images)) {
          td[k] <- inv_times_a_b[k] - inv_times_a_b[k - 1] - ta
          e_td[k] <- exp(-td[k] / t1_vector)
        }
      }

      for (k in seq(n_images))
      {
        cos_alpha_e1[k] <- (cos(flip_angle_a_b_rad[k] / 2)) ^ 2 *
          (e_1a * e_1b) -
          (sin(flip_angle_a_b_rad[k] / 2)) ^ 2 *
          (e_2a * e_1b)

        one_minus_e1[k] <- (1 - e_1a) *
          cos(flip_angle_a_b_rad[k] / 2) *
          e_1b +
          (1 - e_1b)

        sin_alpha[k] <- sin(flip_angle_a_b_rad[k] / 2) *
          cos(flip_angle_a_b_rad[k] / 2) *
          (e_1a + e_2a)
      }
    }

    if (isTRUE(normal_sequence_type))
    {
      e_1 <- exp(-flash_tr / t1_vector)
      ta <- num_z_slices * flash_tr
      ta_bef <- num_z_bef * flash_tr
      ta_aft <- num_z_aft * flash_tr

      td <- zeros(1, n_images + 1)
      td[1] <- inv_times_a_b[1] - ta_bef[1]

      e_td <- zeros(1, n_images + 1)
      e_td[1] <- exp(-td[1] / t1_vector)

      td[n_images + 1] <- mprage_tr -
        inv_times_a_b[n_images] -
        ta_aft

      e_td[n_images + 1] <- exp(-td[n_images + 1] / t1_vector)

      cos_alpha_e1 <- matrix()
      one_minus_e1 <- matrix()
      sin_alpha <- matrix()

      if (n_images > 1) {
        for (k in seq(2, n_images)) {
          td[k] = inv_times_a_b[k] - inv_times_a_b[k - 1] - ta
          e_td[k] = exp(-td[k] / t1_vector)
        }
      }
      for (k in seq(1, n_images)) {
        cos_alpha_e1[k] <- cos(flip_angle_a_b_rad[k]) * e_1
        one_minus_e1[k] <- 1 - e_1
        sin_alpha[k] <- sin(flip_angle_a_b_rad[k])

      }
    }
    # Steady state calculations
    mz_steady_state_numerator <- m0 * (1 - e_td[1])

    for (k in seq(n_images)) {
      # term relative to acquisition
      mz_steady_state_numerator <-
        mz_steady_state_numerator * cos_alpha_e1[k] ^ num_z_slices +
        m0 * (1 - e_1) * (1 - (cos_alpha_e1[k]) ^ num_z_slices) / (1 - cos_alpha_e1[k])

      # term relative to relaxation after
      mz_steady_state_numerator <-
        mz_steady_state_numerator * e_td[k + 1] + m0 * (1 - e_td[k + 1])
    }

    mz_steady_state_denominator <-
      1 + inversion_efficiency * ((prod(cos_alpha_e1)) ^ (num_z_slices)) * prod(e_td)

    mz_steady_state <-
      mz_steady_state_numerator / mz_steady_state_denominator

    # allocate
    signal <- vector()

    m <- 1
    temp <- (-inversion_efficiency *
               mz_steady_state *
               e_td[m] +
               m0 *
               (1 - e_td[m])) *
      ((cos_alpha_e1[m]) ^ (num_z_bef)) +
      m0 *
      (1 - e_1[m]) *
      (1 - ((cos_alpha_e1[m]) ^ (num_z_bef))) /
      (1 - (cos_alpha_e1[m]))

    signal[m] <- sin_alpha[m] * (temp)

    if (n_images > 1) {
      for (m in seq(2, n_images))
      {
        temp <- temp *
          ((cos_alpha_e1[m - 1]) ^ (num_z_aft)) +
          m0 *
          (1 - e_1) *
          (1 - ((cos_alpha_e1[m - 1]) ^ (num_z_aft))) /
          (1 - (cos_alpha_e1[m - 1]))

        temp <- (temp * e_td[m] +
                   m0 *
                   (1 - e_td[m])) *
          ((cos_alpha_e1[m]) ^ (num_z_bef)) +
          m0 *
          (1 - e_1) *
          (1 - ((cos_alpha_e1[m]) ^ (num_z_bef))) /
          (1 - (cos_alpha_e1[m]))

        signal[m] <- sin_alpha[m] * (temp)
      }
    }

    return(signal)
  }
