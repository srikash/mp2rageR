#' #' @title mprage_func
#' #'
#' #' @description
#' #' \code{mprage_func} is part of the \pkg{mp2rageR} package.
#' #'
#' #' \pkg{mp2rageR} is the R implementation of MATLAB code made publicly available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#' #'
#' #' For methodological details, see \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.1910150117}{MPRAGE paper} and \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE paper}.
#' #'
#' #' @author Sriranga Kashyap
#' #'
#' #' @export
#' #' @return Data as a numerical array
#' #' @importFrom pracma strcmp isempty linspace
#' #' @examples
#' #' Signal <- mprage_func(mprage_tr,
#' #'                       inv_times_a_b,
#' #'                       num_z_slices,
#' #'                       flash_tr,
#' #'                       flip_angle_a_b_deg,
#' #'                       sequence_type="normal",
#' #'                       t1s,
#' #'                       ...)
#'
#' mp2rage_func <-
#'   function(mprage_tr = NULL,
#'            inv_times_a_b = NULL,
#'            flip_angle_a_b_deg = NULL,
#'            num_z_slices = NULL,
#'            flash_tr = NULL,
#'            sequence_type = "normal",
#'            inversion_efficiency = 0.96,
#'            b0 = 7,
#'            inv1_mag = NULL,
#'            inv1_phs = NULL,
#'            inv2_mag = NULL,
#'            inv2_phs = NULL,
#'            b1_map = NULL)
#'     {
#'       if (!is.null(b1_map)) {
#'             nii_b1 = readnii(b1_map)
#'             data_b1 = b1_map@.Data
#'         }
#'       if (!is.null(inv1_phs)) {
#'             nii_inv1_phs = readnii(inv1_phs)
#'             data_inv1_phs = inv1_phs@.Data
#'         }
#'       if (!is.null(inv2_phs)) {
#'             nii_inv2_phs = readnii(inv2_phs)
#'             data_inv2_phs = inv2_phs@.Data
#'         }
#'
#'   }
