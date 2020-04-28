param_list_mp2rage <- list(mprage_tr=5.0,
flash_tr=6.9e-3,
inv_times_a_b=c(900e-3,2750e-3),
flip_angle_a_b_deg=c(5,3),
num_z_slices=c(120,120))

param_list_sa2rage <- list(mprage_tr=2.4,
                           flash_tr=2.2e-3,
                           inv_times_a_b=c(58e-3,1800e-3),
                           flip_angle_a_b_deg=c(4,10),
                           num_z_slices=c(22,38))

other_params <- list(sequence_type="normal",
b0 = 7,
m0 = 1,
inversion_efficiency = 0.96,
n_images = 2,
all_data = 1)
