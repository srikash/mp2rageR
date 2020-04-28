#' @title mp2rage_b1_correction_short
#'
#' @description
#' \code{mp2rage_b1_correction_short} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of MATLAB code made publicly available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#'
#' For methodological details, see \href{https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.1910150117}{MPRAGE paper} and \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE paper}.
#'
#' @author Sriranga Kashyap
#'
#' @export
#' @return Data as a numerical array
#' @importFrom pracma strcmp isempty eps
#' @importFrom stats spline
#' @importFrom signal interp1
#' @importFrom neurobase readnii writenii
#' @examples
#' list_of_b1_corr_t1_corr_uni_corr <- mp2rage_b1_correction_short(in_b1_map,
#' is_sa2rage = TRUE,
#' in_uni,
#' param_list_sa2rage,
#' param_list_mp2rage,
#' inversion_efficiency = 0.96)
mp2rage_b1_correction_short <-
  function(in_b1_map,
           is_sa2rage = TRUE,
           # if FALSE, assumes TFL
           in_uni,
           param_list_sa2rage,
           param_list_mp2rage,
           inversion_efficiency = 0.96)
  {
    if (isTRUE(is_sa2rage)) {
      if (is.null(in_b1_map)) {
        stop("b1_map is not specified")
      } else if (is.null(param_list_sa2rage)) {
        stop("param_list_sa2rage is not specified")
      } else {
        nii_b1 = readnii(in_b1_map)
        data_b1 = nii_b1@.Data
        nii_uni = readnii(in_uni)
        data_uni = nii_uni@.Data
      }
    }
    else {
      if (is.null(in_b1_map)) {
        stop("b1_map is not specified")
      } else {
        print("is_sa2rage = FALSE")
        nii_b1 = readnii(in_b1_map)
        data_b1 = nii_b1@.Data

        nii_uni = readnii(in_uni)
        data_uni = nii_uni@.Data
      }
    }

    # sanity check to see how B1 sensitive your sequence was
    b1_vector <- seq(from = 0.6, to = 1.4, by = 0.2)
    for (j in seq_along(b1_vector)) {
      list_of_intensity_t1_vector_intensity_before_comb <-
        mp2rage_lookuptable(
          param_list_mp2rage$mprage_tr,
          param_list_mp2rage$flash_tr,
          param_list_mp2rage$inv_times_a_b,
          b1_vector[j] * param_list_mp2rage$flip_angle_a_b_deg,
          param_list_mp2rage$num_z_slices,
          sequence_type =
            "normal",
          b0 = 7,
          m0 = 1,
          inversion_efficiency = 0.96,
          n_images = 2,
          all_data = 0
        )
    }

    mp2rage_amp <-
      list_of_intensity_t1_vector_intensity_before_comb$intensity
    t1_vector_initial <-
      list_of_intensity_t1_vector_intensity_before_comb$t1_vector
    intensity_before_comb <-
      list_of_intensity_t1_vector_intensity_before_comb$intensity_before_comb
    rm(list_of_intensity_t1_vector_intensity_before_comb)

    if (b0 == 7.0) {
      t1_wm <- 1.1
      t1_gm <- 1.85
      t1_csf <- 3.9
      b1_vector <- seq(from = 0.6, to = 1.4, by = 0.2)
    } else if (b0 == 3) {
      t1_wm <- 0.85
      t1_gm <- 1.35
      t1_csf <- 2.8
      b1_vector <- seq(from = 0.8, to = 1.2, by = 0.1)
    } else {
      t1_wm <- 1.1
      t1_gm <- 1.85
      t1_csf <- 3.9
      b1_vector <- seq(from = 0.6, to = 1.4, by = 0.2)
    }
    # gcf=figure(3);
    # set(gcf,'Color',[1 1 1]);
    # hold off
    # for B1=0.6:0.2:1.4
    # [MP2RAGEamp T1vector IntensityBeforeComb]=
    #   MP2RAGE_lookuptable(
    #     2,
    #     MP2RAGE.TR,
    #     MP2RAGE.TIs,
    #     B1*MP2RAGE.FlipDegrees,
    #     MP2RAGE.NZslices,
    #     MP2RAGE.TRFLASH,
    #     'normal');
    # plot(MP2RAGEamp,T1vector,'color',[0.5 0.5 0.5]*B1,'Linewidth',2)
    # hold on
    # end
    # legend('B1=0.6','B1=0.8','B1=1','B1=1.2','B1=1.4','B1=1.6')

    # definition of range of B1s and T1s and creation of MP2RAGE and Sa2RAGE lookupvector to make sure the input data for the rest of the code is the Sa2RAGEimg and the MP2RAGEimg
    b1_vector <- seq(from = 0.005, to = 1.9, by = 0.05)
    t1_vector <- seq(from = 0.5, to = 5.2, by = 0.05)

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


    # rescale uni image between -0.5 and 0.5
    data_uni_scaled <- data_uni / 4095 - 0.5

    if (isTRUE(is_sa2rage)) {
      list_of_b1_vector_intensity_signal <-
        mp2rage_b1_sa2rage_lookuptable(
          param_list_sa2rage$mprage_tr,
          param_list_sa2rage$flash_tr,
          param_list_sa2rage$inv_times_a_b,
          param_list_sa2rage$flip_angle_a_b_deg,
          param_list_sa2rage$num_z_slices,
          sequence_type = "normal",
          t1_vector = 1.5,
          b0 = 7,
          m0 = 1,
          inversion_efficiency = 0.96,
          n_images = 2
        )

      # create sa2rageimg from b1
      data_b1_scaled <- data_b1 / 1000

      data_sa2rage_vector <- mp2rage_interp_1d(list_of_b1_vector_intensity_signal$b1_vector,
                                               list_of_b1_vector_intensity_signal$intensity,
                                               data_b1_scaled)

      data_sa2rage <-
        array(data_sa2rage_vector, dim(data_b1_scaled))
    }

    # create a lookup table of MP2RAGE intensities as a function of B1 and T1
    mp2rage_matrix <- matrix(
      data = NA,
      nrow = length(b1_vector),
      ncol = length(t1_vector)
    )

    for (k in seq_along(b1_vector)) {
      list_of_intensity_t1_vector = mp2rage_lookuptable(
        param_list_mp2rage$mprage_tr,
        param_list_mp2rage$flash_tr,
        param_list_mp2rage$inv_times_a_b,
        b1_vector[k] * param_list_mp2rage$flip_angle_a_b_deg,
        param_list_mp2rage$num_z_slices,
        sequence_type =
          "normal",
        b0 = 7,
        m0 = 1,
        inversion_efficiency = 0.96,
        n_images = 2,
        all_data = 0
      )

      mp2rage_matrix[k, ] <- mp2rage_interp_1d(list_of_intensity_t1_vector$t1_vector,
                                               list_of_intensity_t1_vector$intensity,
                                               t1_vector)
    }

    mp2rage_matrix[mp2rage_matrix == 0] = NA

    if (isTRUE(is_sa2rage)) {
      # create a lookup table of Sa2RAGE intensities as a function of B1 and T1
      sa2rage_matrix <- matrix(
        data = NA,
        nrow = length(t1_vector),
        ncol = length(b1_vector)
      )

      for (k in seq_along(t1_vector)) {
        list_of_b1_vector_intensity = mp2rage_b1_sa2rage_lookuptable(
          param_list_sa2rage$mprage_tr,
          param_list_sa2rage$flash_tr,
          param_list_sa2rage$inv_times_a_b,
          param_list_sa2rage$flip_angle_a_b_deg,
          param_list_sa2rage$num_z_slices,
          sequence_type =
            "normal",
          t1_vector = t1_vector[k],
          b0 = 7,
          m0 = 1,
          inversion_efficiency = 0.96,
          n_images = 2
        )

        sa2rage_matrix[k, ] <-  signal::interp1(
          x = list_of_b1_vector_intensity$b1_vector,
          y = c(list_of_b1_vector_intensity$intensity),
          xi = b1_vector,
          method = "linear",
          extrap = T
        )
      }
      #imagesc(C = sa2rage_matrix, col = jet.colors(1000))

      ## make the sa2rage_matrix into b1_matrix(t1,ratio)
      sa2rage_vector <-
        seq(min(as.vector(sa2rage_matrix)), max(as.vector(sa2rage_matrix)), length.out = 40)

      b1_matrix <- matrix(
        data = NA,
        nrow = length(t1_vector),
        ncol = length(sa2rage_vector)
      )

      for (k in seq_along(t1_vector)) {
        if (any(is.na(sa2rage_matrix[k,])) == TRUE)
        {
          b1_matrix[k, ] = 0
        }
        else
        {
          temp <- spline(
            x = sa2rage_matrix[k, ],
            y = c(b1_vector),
            xout = sa2rage_vector,
            method = "fmm"
          )

          b1_matrix[k,] <- temp$y

          rm(temp)

        }
      }
      #imagesc(C = b1_matrix, col = jet.colors(1000))
    }

    ## make the mp2rage_matrix into t1_matrix(b1,ratio)
    mp2rage_vector <-
      seq(from = -0.5,
          to = 0.5,
          length.out = 40)

    t1_matrix <- matrix(
      data = NA,
      nrow = length(b1_vector),
      ncol = length(mp2rage_vector)
    )

    for (k in seq_along(b1_vector)) {
      if (isTRUE(any(is.na(mp2rage_matrix[k, ]))))
      {
        temp <- mp2rage_matrix[k,]
        temp[is.na(temp)] <-
          seq(
            from = -0.5 - eps(x = 1),
            to = -1,
            length.out = length(which(is.na(temp) == TRUE))
          )
        temp2 <- signal::interp1(temp,
                                 c(t1_vector),
                                 mp2rage_vector,
                                 method = "linear",
                                 extrap = T)
        t1_matrix[k,] <- temp2
        rm(temp)
        rm(temp2)
      }
      else
      {
        temp <- spline(
          x = mp2rage_matrix[k,],
          y = c(t1_vector),
          xout = mp2rage_vector,
          method = "fmm"
        )

        t1_matrix[k, ] <- temp$y

        rm(temp)
      }
    }
    #imagesc(C = t1_matrix, col = jet.colors(1000))

    ## iterative correction of T1 and B1 estimates
    t1_temp_img <-
      array(data = 1.5, dim = dim(data_uni)) #t1_temp.img

    if (isTRUE(is_sa2rage)) {
      b1_temp_img <- data_sa2rage #b1_temp.img
      data_sa2rage[which(is.na(data_sa2rage))] <- -0.5 #sa2rage.img

      for (k in seq(3))
      {
        print(paste(">>>> Iteration", as.character(k)))
        #### B1
        b1_temp_img_mid = b1_temp_img[, dim(b1_temp_img)[2] / 2,]
        b1_temp_interp <- akima::bicubic(
          x = sa2rage_vector,
          #40
          y = c(t1_vector),
          # 95
          z = t(b1_matrix),
          # dimensions different in MATLAB code
          # 95x40
          x0 = data_sa2rage,
          y0 = t1_temp_img
        )

        # output is vector, needs to reshaped to image size
        data_b1_corr <-
          array(data = b1_temp_interp$z, dim = dim(data_sa2rage))
        data_b1_corr[which(is.na(data_b1_corr))] = 2

        ## Rescale B1
        data_b1_corr <- (data_b1_corr) * 1000

        b1_temp_img_mid2 = data_b1_corr[, dim(data_b1_corr)[2] / 2,]

        #image((b1_temp_img_mid2 - b1_temp_img_mid) * 1000, col = gray.colors(1000))

        #### T1
        t1_temp_img_mid = t1_temp_img[, dim(t1_temp_img)[2] / 2,]
        t1_temp_interp <- akima::bicubic(
          x = mp2rage_vector,
          #40
          y = c(b1_vector),
          # 95
          z = t(t1_matrix),
          # 95x40
          x0 = data_uni_scaled,
          y0 = data_b1_corr
        )
        # output is vector, needs to reshaped to image size
        data_t1_corr <-
          array(data = t1_temp_interp$z, dim = dim(data_uni))
        data_t1_corr[which(is.na(data_t1_corr))] = 4
        ## Rescale T1
        data_t1_corr <- (data_t1_corr) * 1000

        t1_temp_img_mid2 = data_t1_corr[, dim(data_t1_corr)[2] / 2,]

        #image((t1_temp_img_mid2 - t1_temp_img_mid) * 1000, col = gray.colors(1000))


      }
    } else {
      #### T1
      t1_temp_img_mid = t1_temp_img[, dim(t1_temp_img)[2] / 2, ]
      t1_temp_interp <- akima::bicubic(
        x = mp2rage_vector,
        #40
        y = c(b1_vector),
        # 95
        z = t(t1_matrix),
        # 95x40
        x0 = data_uni_scaled,
        y0 = data_b1
      )
      # output is vector, needs to reshaped to image size
      data_t1_corr <-
        array(data = t1_temp_interp$z, dim = dim(data_uni))

      data_t1_corr[which(is.na(data_t1_corr))] = 4

      ## Rescale T1
      data_t1_corr <- (data_t1_corr) * 1000

      t1_temp_img_mid2 = data_t1_corr[, dim(data_t1_corr)[2] / 2, ]

      #image((t1_temp_img_mid2 - t1_temp_img_mid) * 1000, col = gray.colors(1000))
    }


    ## Make the B1 corrected UNI image and put B1 and T1 in ms scale
    list_of_intensity_t1_vector = mp2rage_lookuptable(
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
    nii_uni_corr <- nii_uni
    data_uni_corr <- data_uni_scaled

    data_uni_corr <-
      array(
        data = signal::interp1(
          x = list_of_intensity_t1_vector$t1_vector,
          y = c(list_of_intensity_t1_vector$intensity),
          xi = as.vector(data_t1_corr),
          method = "linear",
          extrap = T
        ),
        dim = size(data_t1_corr)
      )

    data_uni_corr[which(is.na(data_uni_corr))] <- -0.5
    data_uni_corr_rescaled <-
      round(4095 * (data_uni_corr + 0.5))


    if (isTRUE(is_sa2rage)) {
      return(
        list(
          b1_corr = data_b1_corr,
          t1_corr = data_t1_corr,
          uni_corr = data_uni_corr_rescaled
        )
      )
    } else {
      return(list(
        b1_corr = data_b1,
        t1_corr = data_t1_corr,
        uni_corr = data_uni_corr_rescaled
      ))
    }

  }
