#' @title mp2rage_plot_properties
#'
#' @description
#' \code{mp2rage_plot_properties} is part of the \pkg{mp2rageR} package.
#'
#' \pkg{mp2rageR} is the R implementation of MaTLab code made publicly
#' available by \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{José P. Marques}.
#'
#' For methodological details, see \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE paper}.
#'
#' Please cite the aforementioned paper if you used this method in your work.
#'
#' @author Sriranga Kashyap
#'
#' @export
#' @return List of numerical arrays (in order): intensity, t1_vector, intensity_before_comb
#' @importFrom pracma isempty
#' @examples
#' list_of_intensity_t1_vector_intensity_before_comb <- mp2rage_plot_properties(mprage_tr,
#'                                                                          inv_times_a_b,
#'                                                                          flip_angle_a_b_deg,
#'                                                                          num_z_slices,
#'                                                                          flash_tr,
#'                                                                          sequence_type,
#'                                                                          n_images = 2,
#'                                                                          b0 = 7.0,
#'                                                                          m0 = 1.0,
#'                                                                          inversion_efficiency = 0.96,
#'                                                                          all_data = 0,
#'                                                                          t1_vector = NULL)

mp2rage_plot_properties <- function(mprage_tr,
                                 inv_times_a_b,
                                 flip_angle_a_b_deg,
                                 num_z_slices,
                                 flash_tr,
                                 sequence_type,
                                 n_images = 2,
                                 b0 = 7,
                                 m0 = 1,
                                 inversion_efficiency = 0.96,
                                 all_data = 1,
                                 t1_vector = NULL) {

# define signal and noise function as in the paper
signal_res <- function (x1 , x2) {
  return(((x1 * x2) / (x2 ^ 2 + x1 ^ 2)))
}

noise_res <- function (x1 , x2) {
  return(sqrt(((x2 ^ 2 - x1 ^ 2) ^ 2) / ((x2 ^ 2 + x1 ^ 2) ^ 3)))
}

b0 <-7.0 # remove later as it will be part of function input

if (b0 == 7.0) {
  t1_wm <- 1.1
  t1_gm <- 1.85
  t1_csf <- 3.9
  b1_range <- seq(from = 0.6, to = 1.4, by = 0.2)
} else if (b0 == 3) {
  t1_wm <- 0.85
  t1_gm <- 1.35
  t1_csf <- 2.8
  b1_range <- seq(from = 0.8, to = 1.2, by = 0.1)
} else {
  t1_wm <- 1.1
  t1_gm <- 1.85
  t1_csf <- 3.9
  b1_range <- seq(from = 0.6, to = 1.4, by = 0.2)
}
mp2rage_amp = list()

for (k in seq_along(b1_range)) {
  list_of_intensity_t1_vector_intensity_before_comb <- mp2rage_lookuptable(mprage_tr,
      inv_times_a_b,
      b1_range[k] * flip_angle_a_b_deg,
      num_z_slices,
      flash_tr,
      sequence_type,
      n_images = 2,
      b0 = 7,
      m0 = 1,
      inversion_efficiency = 0.96,
      all_data = 1,
      t1_vector = NULL)

  mp2rage_amp[[k]] <-list_of_intensity_t1_vector_intensity_before_comb$intensity
  # plot(mp2rage_amp, t1_vector, type = 'l') #, 'color', [0.5 0.5 0.5] * b1_range[k], 'Linewidth', 2)
  # par(new=T)

  pos_wm = which.min(abs(t1_wm - t1_vector))
  pos_gm = which.min(abs(t1_gm - t1_vector))
  pos_csf = which.min(abs(t1_csf - t1_vector))

  signal <-
    signal_res(intensity_before_comb[c(pos_wm, pos_gm, pos_csf), 1],
               intensity_before_comb[c(pos_wm, pos_gm, pos_csf), 2])

  noise <-
    noise_res(intensity_before_comb[c(pos_wm, pos_gm, pos_csf), 1],
              intensity_before_comb[c(pos_wm, pos_gm, pos_csf), 2])

#  contrast[k] <- num2str(1000 * sum((signal[2:length(signal)] - signal[1:(length(signal) - 1)]) / sqrt(noise[2:length(noise)] ^ 2 + noise[1:(length(noise) - 1)]) ^ 2)) / sqrt(mp2rage_tr))

# legendcell[k] = ['B1= ', num2str(b1_range[k])]

}

# # legend('B1=0.6', 'B1=0.8', 'B1=1', 'B1=1.2', 'B1=1.4')
# legend(legendcell)
# # examples of T1 values at 3T
# plot([-0.5 0.5],
#      [T1CSF T1CSF
#       T1GM T1GM
#       T1WM T1WM]','Linewidth',2)
# text(0.35,T1WM,'White Matter')
# text(0.35,T1GM,'Grey Matter')
# # text(-0.3,(T1CSF+T1GM)/2,['Contrast over B1 range = ',num2str(round(1000*Contrast))])
# text(-0.3,(T1CSF+T1GM)/2,['Contrast over B1 range' ])
# text(0,(T1CSF+T1GM)/2,Contrast)
#
#
# ylabel('T1');
# xlabel('MP2RAGE');
}


## bangu the great's sexy graph ;)
library(ggplot2)
testdata = data.frame(mp2rage_amp, t1_vector=t1_vector)
colnames(testdata)[1:5] = paste("b1_", num2str(b1_range), sep="")
ggplot(data=testdata, aes(y=t1_vector))+
  geom_line(aes(x = testdata[,1]), color = "darkred") +
  geom_smooth()
  geom_line(aes(x = testdata[,2]), color = "darkred") +
  geom_line(aes(x = testdata[,3]), color = "darkred") +
  geom_line(aes(x = testdata[,4]), color = "darkred") +
  geom_line(aes(x = testdata[,5]), color = "darkred")
