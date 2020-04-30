#' @title mp2rage_register_sa2rage
#'
#' @description
#' \code{mp2rage_register_sa2rage} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{MATLAB code} made publicly available by José P. Marques.
#'
#' Uses \href{https://fsl.fmrib.ox.ac.uk/fsl/fslwiki}{fsl} via \href{https://cran.r-project.org/web/packages/fslr/index.html}{fslr}.
#'
#' Please cite the linked papers if you used these methods in your work.
#'
#' @author Sriranga Kashyap
#'
#' @param in_sa2rage_inv2 (required) path to the magnitude Sa2RAGE INV2 NIfTI file
#' @param in_sa2rage_b1_map (required) path to the magnitude Sa2RAGE B1 map NIfTI file
#' @param in_mp2rage_inv2 (required) path to the magnitude MP2RAGE INV2 NIfTI file
#' @param out_reg_b1_map (required) path to the registered Sa2RAGE B1 map NIfTI file
#'
#' @export
#'
#' @return Saves Sa2RAGE B1 map NIfTI file registered and resampled into MP2RAGE space.
#'
#' @importFrom fslr flirt flirt_apply
#'
#' @examples
#' \dontrun{
#' mp2rage_register_sa2rage(
#'        in_sa2rage_inv2 = "sub-01_Sa2RAGE_INV2.nii.gz",
#'        in_sa2rage_b1_map = "sub-01_Sa2RAGE_B1map.nii.gz",
#'        in_mp2rage_inv2 = "sub-01_MP2RAGE_INV2.nii.gz",
#'        out_reg_b1_map = "sub-01_Sa2RAGE_B1map_reg-to-mp2rage.nii.gz")}
mp2rage_register_sa2rage <-
  function(in_sa2rage_inv2,
           in_sa2rage_b1_map,
           in_mp2rage_inv2,
           out_reg_b1_map) {
    # Step 1: Register Sa2RAGE INV2 image to MP2RAGE INV2 image
    flirt(
      infile = in_sa2rage_inv2,
      reffile = in_mp2rage_inv2,
      omat = "fsl_sa2rage-reg-to-mp2rage.mat",
      opts = "-nosearch"
    )
    # Step 2: Apply transform to B1 map
    flirt_apply(
      infile = in_sa2rage_b1_map,
      reffile = in_mp2rage_inv2,
      initmat = "fsl_sa2rage-reg-to-mp2rage.mat",
      outfile = out_reg_b1_map,
      retimg = FALSE
    )
  }
