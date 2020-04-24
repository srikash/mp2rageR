#' @title mp2rage_b1_correction
#'
#' @description
#' \code{mp2rage_b1_correction} is part of the \pkg{mp2rageR} package.
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
mp2rage_b1_correction <-
  function(param_list_mp2rage,
           in_b1_map = NULL,
           is_sa2rage = TRUE)
  {
    # """ This function corrects the bias-field corrected T1-weighted image (`t1w_uni`-attribute)
    #       and the quantitative T1 map (`t1map`-attribute) for B1 inhomogenieties using a B1 field map.
    #       (see Marques and Gruetter, 2013).
    #       It assumes that the B1 field map is either a ratio of the real and intended
    #       flip angle (range of approximately 0 - 2) *or* the percentage of the real
    #       vs intended flip angle (range of approximately 0 - 200).
    #
    #       If the B1 map has a different resolution, it is resampled to the resolution
    #       of INV1 and INV2. *This function assumes your MP2RAGE images and the B1 map
    #       are in the same space*.
    #
    #       If the B1 map is not immediately acquired after the MP2RAGE sequence,
    #       you should register the (magnitude image corresponding to) the B1 map to
    #       INV1/INV2 first.
    #
    #       The corrected images are stored in the `t1w_uni_b1_corrected` and the
    #       `t1_b1_corrected`-attributes as well as returned as a tuple
    #
    #
    #       Args:
    #           B1 (filename): B1 field map, either as a ratio or as a percentage. If
    #                          set to None, use self.B1 (set when the MP2RAGE class was
    #                          initialized).
    #           check_B1_range (bool): whether the fuction should check whether the range
    #                                  of the B1 fieldmap makes sense (centered at 1, range of
    #                                  roughly 0-2 or 0-200).
    #
    #       Returns:
    #           (tuple): tuple containing:
    #
    #               t1w_uni_b1_corrected: A T1-weighted image corrected for B1 inhomogenieties
    #               t1map_b1_corrected: A quantiative T1-weighted image corrected for B1 inhomogenieties
    #
    #
    #       """

    if (is.null(file_b1_map)) {
      stop("B1 map is not specified")
    } else {
      print(">>>> Loading B1 map")
      nii_b1 = readnii(file_b1_map)
      data_b1 = nii_b1@.Data
    }

    if (isTRUE(sa2rage)) {
      data_b1_deg = data_b1 / 100 # in degrees
    }

    # creates a lookup table of MP2RAGE intensities as a function of B1 and T1
    b1_vector = seq(from = 0.005, to = 2.5, by = 0.005)
    t1_vector =  seq(from = 0.05, to = 5.0, by = 0.05)

    # mp2rage_matrix = zeros(length(b1_vector), length(t1_vector))
    #
    # list_of_intensity_t1_vector <-
    #   list(intensity, t1_vector, intensity_before_comb)
    # for n in seq_along(b1_vector) {
    #   for b1 in seq_along(b1_vector) {
    #     flip_angle_corr = b1 * as.array(flip_angle_a_b)
    #
    #     list_of_intensity_t1_vector = mp2rage_lookuptable(
    #       mprage_tr,
    #       inv_times_a_b,
    #       flip_angle_corr,
    #       num_z_slices,
    #       flash_tr,
    #       sequence_type,
    #       n_images = 2,
    #       inversion_efficiency = 0.96,
    #       b0 = 7,
    #       all_data = 0
    #     )
    #
    #     temp = interp1(t1_vector, intensity, method = "linear")
    #     mp2rage_matrix[n,] = temp[t1_vector]
    #
    #     #     return MP2RAGEmatrix
    #
    #     # make the matrix  MP2RAGEMatrix into T1_matrix(B1,ratio)
    #     npoints = 40
    #
    #     MP2RAGE_vector = np.linspace(-0.5, 0.5, npoints)
    #
    #
    #     T1matrix = np.zeros((len(B1_vector), npoints))
    #
    #     for k, b1val in enumerate(B1_vector):if np.isnan(MP2RAGEmatrix[k,:]).any():signal = MP2RAGEmatrix[k,:].copy()
    #     signal[np.isnan(signal)] = np.linspace(-0.5,-1, np.isnan(signal).sum())
    #
    #     f = interpolate.interp1d(signal,
    #                              T1_vector,
    #                              bounds_error = False,
    #                              fill_value = 'extrapolate')#fill_value='extrapolate')
    #
    #     T1matrix[k,:] = f(MP2RAGE_vector)
    #
    #     else:signal = MP2RAGEmatrix[k,:]
    #     f = interpolate.interp1d(
    #       np.sort(MP2RAGEmatrix[k,:]),
    #       T1_vector[np.argsort(MP2RAGEmatrix[k,:])],
    #       'cubic',
    #       bounds_error = False,
    #       fill_value = 'extrapolate'
    #     )
    #     T1matrix[k, :] = f(MP2RAGE_vector)
    #
    #
    #     # *** Create correted T1 map ***
    #     # Make interpolation function that gives T1, given B1 and T1w signal
    #     f = interpolate.RectBivariateSpline(B1_vector,
    #                                         MP2RAGE_vector,
    #                                         T1matrix,
    #                                         kx = 1,
    #                                         ky = 1)
    #
    #
    #     x = B1.get_data()
    #
    #     # Rescale T1w signal to [-.5, .5]
    #     y = self.t1w_uni.get_data() / 4095 - .5
    #
    #     # Precache corrected T1 map
    #     t1c = np.zeros_like(x)
    #
    #     # Make a mask that excludes non-interesting voxels
    #     mask = (x != 0) & (y != 0) & ~ np.isnan(y)
    #
    #     # Interpolate T1-corrected map
    #     t1c[mask] = f(x[mask], y[mask], grid = False)
    #     self.t1map_b1_corrected = nb.Nifti1Image(t1c * 1000, self.t1map.affine)
    #
    #     # *** Create corrected T1-weighted image ***
    #     Intensity, T1vector,   _   = MP2RAGE_lookuptable(
    #       self.MPRAGE_tr,
    #       self.invtimesAB,
    #       self.flipangleABdegree,
    #       self.nZslices,
    #       self.FLASH_tr,
    #       self.sequence,
    #       nimages = 2,
    #       inversion_efficiency = self.inversion_efficiency,
    #       B0 = self.B0,
    #       all_data = 0
    #     )
    #
    #     f = interpolate.interp1d(T1vector,
    #                              Intensity,
    #                              bounds_error = False,
    #                              fill_value = -0.5)
    #     t1w_uni_corrected = (f(t1c) + .5) * 4095
    #     self.t1w_uni_b1_corrected = nb.Nifti1Image(t1w_uni_corrected, self.t1w_uni.affine)
    #
    #     return self.t1map_b1_corrected, self.t1w_uni_b1_corrected
  }
