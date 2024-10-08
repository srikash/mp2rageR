% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mp2rage_b1_sa2rage_lookuptable.R
\name{mp2rage_b1_sa2rage_lookuptable}
\alias{mp2rage_b1_sa2rage_lookuptable}
\title{mp2rage_b1_sa2rage_lookuptable}
\usage{
mp2rage_b1_sa2rage_lookuptable(
  mprage_tr,
  flash_tr,
  inv_times_a_b,
  flip_angle_a_b_deg,
  num_z_slices,
  sequence_type = "normal",
  t1_vector = 1.5,
  b0 = 7,
  m0 = 1,
  inversion_efficiency = 0.96,
  n_images = 2
)
}
\arguments{
\item{mprage_tr}{(required) \emph{TR} of the Sa2RAGE in s}

\item{flash_tr}{(required) \emph{TR} of the FLASH readout in s}

\item{inv_times_a_b}{(required) \eqn{c(TI1,TI2)}}

\item{flip_angle_a_b_deg}{(required) \emph{c(\eqn{\alpha}1,\eqn{\alpha}2)}}

\item{num_z_slices}{(required) calculate as: \eqn{base_resolution * c(phase_partial_fourier-0.5,0.5)/phase_ipat+c(ref_lines/2,ref_lines/2) * (1-1/phase_ipat)}}

\item{sequence_type}{(has default) \emph{"normal"}, set to \emph{NULL} for water excitation}

\item{t1_vector}{(has default) \emph{1.5}, provide input vector as necessary}

\item{b0}{(has default) \emph{7.0} can be \emph{3.0}}

\item{m0}{(has default) \emph{1.0}}

\item{inversion_efficiency}{(has default) \emph{0.96} is the efficiency of the Siemens MP2RAGE inversion pulse}

\item{n_images}{(has default) \emph{2}, change if required}
}
\value{
List of numerical arrays (in order): intensity, t1_vector, intensity_before_comb
}
\description{
\code{mp2rage_b1_sa2rage_lookuptable} is part of the \pkg{mp2rageR} package.

\pkg{mp2rageR} is the R implementation of \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{MATLAB code} made publicly available by José P. Marques.

For methodological details, see \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0069294}{Sa2RAGE correction}.

Please cite the linked papers if you used these methods in your work.
}
\examples{

sa2rage_params <- list(
                   mprage_tr = 2.4,
                   flash_tr = 2.2e-3,
                   inv_times_a_b = c(58e-3,1800e-3),
                   flip_angle_a_b_deg = c(4,10),
                   num_z_slices = c(22,38)
                   )
\dontrun{
out_list <- mp2rage_b1_sa2rage_lookuptable(
                   sa2rage_params$mprage_tr,
                   sa2rage_params$flash_tr,
                   sa2rage_params$inv_times_a_b,
                   sa2rage_params$flip_angle_a_b_deg,
                   sa2rage_params$num_z_slices,
                   sequence_type = "normal",
                   t1_vector = 1.5,
                   b0 = 7,
                   m0 = 1,
                   inversion_efficiency = 0.96,
                   n_images = 2
                   )}
}
\author{
Sriranga Kashyap
}
