#' @title mp2rage_register_sa2rage2
#'
#' @description
#' \code{mp2rage_register_sa2rage2} is part of the \pkg{mp2rageR} package.
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
#' @importFrom RNifti readNifti writeNifti
#' @importFrom RNiftyReg niftyreg forward applyTransform
#'
#' @examples
#' \dontrun{
#' mp2rage_register_sa2rage22(
#'        in_sa2rage_inv2 = "sub-01_Sa2RAGE_INV2.nii.gz",
#'        in_sa2rage_b1_map = "sub-01_Sa2RAGE_B1map.nii.gz",
#'        in_mp2rage_inv2 = "sub-01_MP2RAGE_INV2.nii.gz",
#'        out_reg_b1_map = "sub-01_Sa2RAGE_B1map_reg-to-mp2rage.nii.gz")}

mp2rage_register_sa2rage2 <-
  function(in_sa2rage_inv2,
           in_sa2rage_b1_map,
           in_mp2rage_inv2,
           out_reg_b1_map) {
    # Step 1: Register Sa2RAGE INV2 image to MP2RAGE INV2 image
    nii_sa2rage_inv2 <- RNifti::readNifti(file = in_sa2rage_inv2)
    nii_mp2rage_inv2 <- RNifti::readNifti(file = in_mp2rage_inv2)
    reg_out <- RNiftyReg::niftyreg(
      source = nii_sa2rage_inv2,
      target = nii_mp2rage_inv2,
      scope = "rigid",
      interpolation = 1,
      estimateOnly = TRUE,
      threads = 4,
      precision = "single"
    )
    nii_sa2rage_b1 <- RNifti::readNifti(in_sa2rage_b1_map)
    nii_sa2rage_b1_linreg <- RNiftyReg::applyTransform(
      transform = RNiftyReg::forward(reg_out),
      x = nii_sa2rage_b1,
      interpolation = 1)
    RNifti::writeNifti(image = nii_sa2rage_b1_linreg,
                       file = out_reg_b1_map)
  }
