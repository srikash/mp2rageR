#' @title mp2rage_estimate_t1_m0
#'
#' @description
#' \code{mp2rage_estimate_t1_m0} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{MATLAB code} made publicly available by José P. Marques.
#'
#' For methodological details, see \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE}.
#'
#' Please cite the linked papers if you used these methods in your work.
#'
#' @author Sriranga Kashyap
#'
#' @param in_uni_data (required) path to the magnitude MP2RAGE UNI NIfTI file
#' @param in_inv2_data (has default) path to the magnitude MP2RAGE INV2 NIfTI file, not required for only T1 estimation
#' @param param_list_mp2rage (required) list of parameters of MP2RAGE acquisition
#' @export
#' @return List of numerical arrays (not NIfTI)
#' @importFrom signal interp1
#' @examples
#' \dontrun{
#' list_t1_m0_maps <- mp2rage_estimate_t1_m0(
#'                          in_uni_data,
#'                          in_inv2_data,
#'                          param_list_mp2rage)}
mp2rage_estimate_t1_m0 <- function(in_uni_data,
                                   in_inv2_data = NULL,
                                   param_list_mp2rage) {
  list_of_intensity_t1_vector_intensity_before_comb <-
    mp2rage_lookuptable(
      param_list_mp2rage$mprage_tr,
      param_list_mp2rage$flash_tr,
      param_list_mp2rage$inv_times_a_b,
      param_list_mp2rage$flip_angle_a_b_deg,
      param_list_mp2rage$num_z_slices,
      sequence_type =
        "normal",
      b0 = 7,
      m0 = 1,
      inversion_efficiency = 0.96,
      n_images = 2,
      all_data = 0
    )

  # checks if it is a standard MP2RAGE image, from -0.5 to 0.5
  if (max(abs(as.vector(in_uni_data))) > 1) {
    t1_vector = signal::interp1(
      list_of_intensity_t1_vector_intensity_before_comb$intensity,
      list_of_intensity_t1_vector_intensity_before_comb$t1_vector,-0.5 + 1 / 4095 * as.vector(in_uni_data),
      method = "linear",
      extrap = TRUE
    )

  } else {
    t1_vector = signal::interp1(
      list_of_intensity_t1_vector_intensity_before_comb$intensity,
      list_of_intensity_t1_vector_intensity_before_comb$t1_vector,
      as.vector(in_uni_data),
      method = "linear",
      extrap = TRUE
    )
  }

  t1_vector[is.na(t1_vector)] <- 0
  t1_vector[t1_vector < 0] <- 0
  t1_vector[t1_vector > 5] <- 0

  t1_map <- array(data = t1_vector * 1000, dim = dim(in_uni_data))

  if (is.null(in_inv2_data)) {
    return(list(t1_map = t1_map))
  } else {
    intensity_before_comb <-
      signal::interp1(
        list_of_intensity_t1_vector_intensity_before_comb$t1_vector,
        list_of_intensity_t1_vector_intensity_before_comb$intensity_before_comb[, 2],
        t1_vector
      )

    m0_map = array(
      data = (as.vector(in_inv2_data) / intensity_before_comb) / 10,
      dim = dim(in_inv2_data)
    )

    return(list(t1_map = t1_map, m0_map = m0_map))
  }

}
