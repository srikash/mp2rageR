#' @title robust_combination
#'
#' @description
#' \code{robust_combination} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of MATLAB code made publicly available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#'
#' For methodological details, see \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE sequence} and \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099676}{MP2RAGE background supression algorithm}.
#'
#' Please cite the aforementioned papers if you used these methods in your work.
#'
#' @author Sriranga Kashyap
#'
#' @param inv1 path to MP2RAGE inv1 NIfTI file (required)
#' @param inv2 path to MP2RAGE inv1 NIfTI file (required)
#' @param uni path to MP2RAGE uni NIfTI file (required)
#' @param uni_out (optional)  path to background suppressed MP2RAGE uni NIfTI file
#' @param regularisation (optional) a scalar (1-10) for background noise supression strength
#' @export
#' @return Data as a numerical array (not NIfTI, unless uni_out is specified)
#' @importFrom neurobase readnii writenii
#' @examples
#' mp2rage_robust_phase_sensitive <- robustCombination("sub-01_inv1.nii.gz","sub-01_inv2.nii.gz","sub-01_uni.nii.gz","sub-01_uni_out.nii.gz",regularisation = 10)
robust_combination <-
  function(inv1,
           inv2,
           uni,
           uni_out = NULL,
           regularisation = 1) {
    # Functions used for robustCombination
    mp2rage_robust_func <- function (inv1 , inv2 , beta) {
      # MP2RAGErobustfunc=@(inv1,inv2,beta)(conj(inv1).*inv2-beta)./(inv1.^2+inv2.^2+2*beta);
      return((Conj(inv1) * inv2 - beta) / (inv1 ^ 2 + inv2 ^ 2 + 2 * beta))
    }

    pos_rootsquares <- function (a , b , c)
    {
      # rootsquares_pos  =@(a,b,c)(-b+sqrt(b.^2 -4 *a.*c))./(2*a);
      return((-b + sqrt(b ^ 2 - (4 * a * c))) / (2 * a))
    }

    neg_rootsquares <- function (a , b , c)
    {
      # rootsquares_neg  =@(a,b,c)(-b-sqrt(b.^2 -4 *a.*c))./(2*a);
      return((-b - sqrt(b ^ 2 - (4 * a * c))) / (2 * a))
    }

    rescale_intensities <-
      function(data) {
        if (min(data) >= 0.0 &
            max(data) >= 0.51) {
          # converts MP2RAGE to -0.5 to 0.5 scale - assumes that it is getting only positive values
          return((data - max(data) / 2) / max(data))
        } else {
          integer_format <- 0
          return(integer_format)
        }
      }

    # Load data from NIfTI filenames
    nii_inv1 <- readnii(inv1)
    nii_inv2 <- readnii(inv2)
    nii_uni <- readnii(uni)

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

    #inv1pos=rootsquares_pos(-MP2RAGEimg.img,inv2img.img,-inv2img.img.^2.*MP2RAGEimg.img);
    data_inv1_pos <- pos_rootsquares(-a, b,-c)
    data_inv1_pos[is.na(data_inv1_pos)] <- 0 # Remove NaNs

    #inv1neg=rootsquares_neg(-MP2RAGEimg.img,inv2img.img,-inv2img.img.^2.*MP2RAGEimg.img);
    data_inv1_neg <- neg_rootsquares(-a, b,-c)
    data_inv1_neg[is.na(data_inv1_neg)] <- 0 # Remove NaNs

    #inv1final=inv1img.img;
    data_inv1_final <- data_inv1_corrected

    #inv1final(abs(inv1img.img-inv1pos)> abs(inv1img.img-inv1neg))=inv1neg(abs(inv1img.img-inv1pos)>abs(inv1img.img-inv1neg));
    data_inv1_final[abs(data_inv1_corrected - data_inv1_pos) > abs(data_inv1_corrected - data_inv1_neg)] <- data_inv1_neg[abs(data_inv1_corrected - data_inv1_pos) > abs(data_inv1_corrected - data_inv1_neg)]

    #inv1final(abs(inv1img.img-inv1pos)<=abs(inv1img.img-inv1neg))=inv1pos(abs(inv1img.img-inv1pos)<=abs(inv1img.img-inv1neg));
    data_inv1_final[abs(data_inv1_corrected - data_inv1_pos) <= abs(data_inv1_corrected - data_inv1_neg)] <- data_inv1_pos[abs(data_inv1_corrected - data_inv1_pos) <= abs(data_inv1_corrected - data_inv1_neg)]

    # Lambda calculation
    # usually the multiplicative factor shouldn't be greater then 10
    # but that is not the ase when the image is bias field corrected.
    # in that case the noise estimated at the edge of the image may
    # not be a good measure to use

    # noiselevel=multiplyingFactor*mean(mean(mean(inv2img.img(1:end,end-10:end,end-10:end))));
    noiseLevel <- regularisation * mean(data_inv2[1:dim(data_inv2)[1], (dim(data_inv2)[2] -
                                                                         10):dim(data_inv2)[2], (dim(data_inv2)[3] - 10):dim(data_inv2)[3]])

    # MP2RAGEimgRobustScanner=MP2RAGErobustfunc(inv1img.img,inv2img.img,noiselevel.^2);
    # MP2RAGEimgRobustPhaseSensitive=MP2RAGErobustfunc(inv1final,inv2img.img,noiselevel.^2);
    mp2rage_robust_phase_sensitive <- mp2rage_robust_func(data_inv1_final, data_inv2, noiseLevel ^ 2)

    # Check for output or return calculations
    if (is.null(uni_out)) {
      return(mp2rage_robust_phase_sensitive)
    } else {
      if (exists("integer_format")) {
        data_uni_clean <- mp2rage_robust_phase_sensitive
      } else
      {
        data_uni_clean <- round(4095 * (mp2rage_robust_phase_sensitive + 0.5))
      }
      nii_uni@.Data <- (data_uni_clean)
      writenii(nii_uni, uni_out)
    }
  }
