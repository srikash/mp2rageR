% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mp2rage_robust_combination.R
\name{mp2rage_robust_combination}
\alias{mp2rage_robust_combination}
\title{mp2rage_robust_combination}
\usage{
mp2rage_robust_combination(
  in_inv1,
  in_inv2,
  in_uni,
  out_uni = NULL,
  regularisation = 1
)
}
\arguments{
\item{in_inv1}{(required) path to the magnitude MP2RAGE INV1 NIfTI file}

\item{in_inv2}{(required) path to the magnitude MP2RAGE INV2 NIfTI file}

\item{in_uni}{(required) path to the magnitude MP2RAGE UNI NIfTI file}

\item{out_uni}{(optional)  path to background suppressed MP2RAGE UNI NIfTI file}

\item{regularisation}{(optional) a scalar (1-10) noise supression strength}
}
\value{
Data as a numerical array (not NIfTI, unless out_uni is specified)
}
\description{
\code{mp2rage_robust_combination} is part of the \pkg{mp2rageR} package.

\pkg{mp2rageR} is the R implementation of \href{https://github.com/JosePMarques/MP2RAGE-related-scripts}{MATLAB code} made publicly available by José P. Marques.

For methodological details, see \href{https://www.sciencedirect.com/science/article/abs/pii/S1053811909010738}{MP2RAGE} and \href{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0099676}{MP2RAGE background supression}.

Please cite the linked papers if you used these methods in your work.
}
\examples{
\dontrun{
mp2rage_robust <- mp2rage_robust_combination(
                            "sub-01_inv1.nii.gz",
                            "sub-01_inv2.nii.gz",
                            "sub-01_uni.nii.gz",
                            "sub-01_out_uni.nii.gz",
                             regularisation = 10)}
}
\author{
Sriranga Kashyap
}
