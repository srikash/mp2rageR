#' @title mprage_func
#'
#' @description
#' \code{mprage_func} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of MATLAB code made publicly available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#'
#' For methodological details, see \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.1910150117}{MPRAGE paper} and \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE paper}.
#'
#' @author Sriranga Kashyap
#'
#' @export
#' @return Data as a numerical array
#' @importFrom pracma strcmp isempty linspace
#' @examples
#' Signal <- mprage_func(mprage_tr,
#'                       inv_times_a_b,
#'                       num_z_slices,
#'                       flash_tr,
#'                       flip_angle_a_b_deg,
#'                       sequence_type="normal",
#'                       t1s,
#'                       ...)

mprage_func <-
  function(mprage_tr,
           inv_times_a_b,
           num_z_slices,
           flash_tr,
           flip_angle_a_b_deg,
           sequence_type="normal",
           t1s,
           n_images=2,
           b0=7,
           m0=1,
           inversion_efficiency=0.96)
    {
    if (isTRUE(sequence_type=="normal")) {
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
    pulse_space <- 1/2/(fat_water_cs*b0*gamma)

    flip_angle_a_b_rad <-
      flip_angle_a_b_deg / 180 * pi #conversion from degrees to radians

    if (n_images != length(flip_angle_a_b_rad)) {
      flip_angle_a_b_rad <- rep(flip_angle_a_b_rad, times = n_images)
    }
    if (n_images != length(flash_tr)) {
      flash_tr <- rep(flash_tr, times = n_images)
    }

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
      e_1 <- exp(-flash_tr / T1s)
      e_1a <- exp(-pulse_space / t1s)
      e_2a <- exp(-pulse_space / 0.06) #60 ms is an estimation of the T2*
      e_1b <- exp(-(flash_tr - pulse_space) / t1s)

      ta <- num_z_slices * flash_tr
      ta_bef <- num_z_bef * flash_tr
      ta_aft <- num_z_aft * flash_tr

      td[1] <- inv_times_a_b[1] - ta_bef
      e_td[1] <- exp(-td[1] / t1s)
      td[nimages + 1] = mprage_tr - inv_times_a_b[n_images] - ta_aft
      e_td[nimages + 1] = exp(-td[nimages + 1] / t1s)

      if (n_images > 1)
        {
        for (k in seq(2,n_images)) {
          td[k] <- inv_times_a_b[k] - inv_times_a_b[k - 1] - ta
          e_td[k] <- exp(-td[k] / t1s)
        }
      }

      for (k in seq(n_images))
        {
          cos_alpha_e1[k] <- (cos(flip_angle_a_b_rad[k] / 2)) ^ 2 *
            (e_1a * e_1b) -
            (sin(flip_angle_a_b_rad[k] / 2)) ^ 2 *
              (e_2a * e_1b)
          one_minus_e1[k] <- (1 - E_1A) *
            cos(flip_angle_a_b_rad[k] / 2) *
            e_1b +
            (1 - e_1b)
          sin_alpha[k] <- sin(fliprad[k] / 2) *
            cos(flip_angle_a_b_rad[k] / 2) *
            (e_1a + e_2a)
          }
    }

    if (isTRUE(normal_sequence_type))
    {
      e_1 <- exp(-flash_tr / t1s)
      ta <- num_z_slices * flash_tr
      ta_bef <- num_z_bef * flash_tr
      ta_aft <- num_z_aft * flash_tr

      td <- zeros(1, n_images+1 )
      e_td <- zeros(1, n_images+1 )

      td[1] <- inv_times_a_b[1] - ta_bef[1]
      e_td[1] <- exp(-td[1] / t1s)

      td[n_images] <- mp2rage_tr -
        inv_times_a_b[n_images - 1] -
        ta_aft[-1]
      e_td[nimages] <- exp(-td[n_images] / t1s)

      if (num_z_slices > 1) {
        td[1:nimages] = inv_times_a_b[1, ] - inv_times_a_b[ ,-1] -(ta_aft[ ,-1] + ta_bef[1, ])
        e_td[1:nimages] = exp(-td[1:n_images] / t1s)
      }
      cos_alpha_e1 <- cos(flip_angle_a_b_rad) * e_1
      one_minus_e1 <- 1 - e_1
      sin_alpha <- sin(flip_angle_a_b_rad)
    }

    mz_steady_state <- 1 / (1 + inversion_efficiency * (prod(cos_alpha_e1))^(num_z_slices) * prod(e_td))
    mz_steady_state_numerator <- m0 * (1 - e_td[1])

    for (i in seq(n_images)) {
      mz_steady_state_numerator <- mz_steady_state_numerator * cos_alpha_e1[i]^num_z_slices +
        m0 * (1 - e_1[i]) * (1 - (cos_alpha_e1[i])^nZslices) / (1 - cos_alpha_e1[i])
      mz_steady_state_numerator <- mz_steady_state_numerator * e_td[i + 1] + m0 * (1 - e_td[i + 1])
    }

    mz_steady_state <- mz_steady_state * mz_steady_state_numerator

    signal <- zeros(n_images)

    m <- 0
    temp <- (-inversion_efficiency *
      mz_steady_state *
      e_td[m] +
      m0 *
        (1-e_td[m])) *
      (cos_alpha_e1[m]) ^ (num_z_bef) +
      m0 *
        (1 - e_1[m]) *
        (1 - (cos_alpha_e1[m]) ^ (num_z_bef)) /
        (1-(cos_alpha_e1[m]))

    signal[1] <- sin_alpha[m] * temp

    for (m in seq_along(nimages))
      {
      temp <- temp *
        (cos_alpha_e1[m - 1])^(num_z_aft) +
        M0 *
          (1 - e_1[m - 1]) *
          (1 - (cos_alpha_e1[m - 1])^(num_z_aft)) /
          (1 - (cos_alpha_e1[m - 1]))

      temp <- (temp * e_td[m] +
        M0 *
          (1 - e_td[m])) *
        (cos_alpha_e1[m])^(num_z_bef) +
        M0 *
          (1 - e_1[m]) *
          (1 - (cos_alpha_e1[m])^(num_z_bef)) /
          (1 - (cos_alpha_e1[m]))

      signal[m] <- sin_alpha[m] * temp
    }

    return(signal)
}
