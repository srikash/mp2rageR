#' @title mp2rage_robust_combination
#'
#' @description
#' \code{mp2rage_robust_combination} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of MATLAB code made publicly available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#'
#' For methodological details, see \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE sequence} and \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099676}{MP2RAGE background supression algorithm}.
#'
#' Please cite the aforementioned papers if you used these methods in your work.
#'
#' @author Sriranga Kashyap
#'
#' @param inv1 (required) path to the magnitude MP2RAGE INV1 NIfTI file
#' @param inv2 (required) path to the magnitude MP2RAGE INV2 NIfTI file
#' @param uni (required) path to the magnitude M?mp2P2RAGE UNI NIfTI file
#' @param out_uni (optional)  path to background suppressed MP2RAGE uni NIfTI file
#' @param regularisation (optional) a scalar (1-10) noise supression strength
#' @export
#' @return Data as a numerical array (not NIfTI, unless out_uni is specified)
#' @importFrom neurobase readnii writenii
#' @examples
#' mp2rage_robust <- mp2rage_robust_combination("sub-01_inv1.nii.gz","sub-01_inv2.nii.gz","sub-01_uni.nii.gz","sub-01_out_uni.nii.gz",regularisation = 10)
mp2rage_robust_combination <-
  function(in_inv1,
           in_inv2,
           in_uni,
           out_uni = NULL,
           regularisation = 1) {
    # Functions used for robust combination
    mp2rage_robust_func <- function (inv1 , inv2 , beta) {
      return((Conj(inv1) * inv2 - beta) / (inv1 ^ 2 + inv2 ^ 2 + 2 * beta))
    }

    pos_rootsquares <- function (a , b , c)
    {
      return((-b + sqrt(b ^ 2 - (4 * a * c))) / (2 * a))
    }

    neg_rootsquares <- function (a , b , c)
    {
      return((-b - sqrt(b ^ 2 - (4 * a * c))) / (2 * a))
    }

    rescale_intensities <-
      function(data) {
        if (min(data) >= 0.0 &
            max(data) >= 0.51) {
          return((data - max(data) / 2) / max(data))
        } else {
          integer_format <- 0
          return(integer_format)
        }
      }

    # Load data from NIfTI filenames
    nii_inv1 <- readnii(in_inv1)
    nii_inv2 <- readnii(in_inv2)
    nii_uni <- readnii(in_uni)

    data_inv1 <- nii_inv1@.Data
    data_inv2 <- nii_inv2@.Data
    data_uni <- nii_uni@.Data

    # Rescale uni between -0.5 to 0.5

    data_uni_rescaled <- rescale_intensities(data_uni)

    # Fix inv1 polarity
    data_inv1 <- data_inv1 * 1.0
    data_inv1_corrected <- sign(data_uni_rescaled) * data_inv1

    # because the MP2RAGE's inv1 and inv2 is a sum of squares data, while the
    # MP2RAGE image is a phase sensitive coil combination.. some more maths has to
    # be performed to get a better inv1 estimate which here is done by assuming
    # both inv2 is closer to a real phase sensitive combination

    a <- data_uni_rescaled
    b <- data_inv2 * 1 # simple conversion to double
    c <- (data_inv2 ^ 2) * data_uni_rescaled

    data_inv1_pos <- pos_rootsquares(-a, b,-c)
    data_inv1_pos[is.na(data_inv1_pos)] <- 0 # Remove NaNs

    data_inv1_neg <- neg_rootsquares(-a, b,-c)
    data_inv1_neg[is.na(data_inv1_neg)] <- 0 # Remove NaNs

    data_inv1_final <- data_inv1_corrected

    data_inv1_final[abs(data_inv1_corrected - data_inv1_pos) > abs(data_inv1_corrected - data_inv1_neg)] <-
      data_inv1_neg[abs(data_inv1_corrected - data_inv1_pos) > abs(data_inv1_corrected - data_inv1_neg)]

    data_inv1_final[abs(data_inv1_corrected - data_inv1_pos) <= abs(data_inv1_corrected - data_inv1_neg)] <-
      data_inv1_pos[abs(data_inv1_corrected - data_inv1_pos) <= abs(data_inv1_corrected - data_inv1_neg)]

    # Lambda calculation
    # usually the multiplicative factor shouldn't be greater then 10
    # but that is not the ase when the image is bias field corrected.
    # in that case the noise estimated at the edge of the image may
    # not be a good measure to use

    noise_level <-
      regularisation * mean(data_inv2[1:dim(data_inv2)[1], (dim(data_inv2)[2] -
                                                              10):dim(data_inv2)[2], (dim(data_inv2)[3] - 10):dim(data_inv2)[3]])

    mp2rage_robust <-
      mp2rage_robust_func(data_inv1_final, data_inv2, noise_level ^ 2)

    # Check for output or return calculations
    if (is.null(out_uni)) {
      return(mp2rage_robust)
    } else {
      if (exists("integer_format")) {
        data_uni_clean <- mp2rage_robust
      } else
      {
        data_uni_clean <- round(4095 * (mp2rage_robust + 0.5))
      }
      nii_uni@.Data <- (data_uni_clean)
      writenii(nii_uni, out_uni)
    }
  }
