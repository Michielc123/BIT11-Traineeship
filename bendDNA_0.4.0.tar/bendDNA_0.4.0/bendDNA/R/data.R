#' Scales for trinucleotide bendability coefficients.
#'
#' Available scales are "con", "conrigid", "dnase", "dnaserigid", "nuc", "nucrigid".
#'
#' @format A data.table with 64 rows and 7 columns:
#' \itemize{
#'  \item{\strong{V1}} : {all 64 possible trinucleotides}
#'  \item{\strong{con, conrigid}} : {consensus DNA bendability scales}
#'  \item{\strong{dnase, dnaserigid}} : {DNA bendability scales derived from DNase I experiments}
#'  \item{\strong{nuc, nucrigid}} : {DNA bendability scales derived from nucleosome positioning experiments}
#' }
#'
#' @details DNA bendability scale based on DNase I digestion experiments
#' was calculated by Brukner et al. in their 1995 paper. All 3D structures of
#' DNA-DNase I complexes show the DNA to be bent. Sequences which are flexible
#' (or already bent) are therefore more acessible to DNase I, and more likely
#' to be cleaved by it. The scale was calculated from cleavage frequencies at
#' different trinucleotides.
#'
#' Data used to calculate the nucleosome-based bendability scale was published
#' by Satchewll, Drew and Travers in 1986. Nucleosome binding requires DNA to
#' be wrapped around the core proteins in a very tight fashion. Hence, frequency
#' with which individual trinucleotides appear in such sharp bends can be correlated
#' with their bending propensity.
#'
#' Gabrielian and Pongor introduced a consensus bendability scale - an average between
#' DNase I-based and nucleosome-based parameters.
#' @source Brukner I., Sanchez R., Suck D. and Pongor S.: Sequence dependent
#' bending propensity of DNA as revealed by DNase I: Parameters for trinucleotides,
#' Embo J. 14 (1995): 1812-1818.
#'
#' Satchwell S.C., Drew H.R. and Travers A.A.: Sequence periodicities in chicken
#' nucleosome core DNA, J. Mol. Biol. 191 (1986): 639-659.
#'
#' Gabrielian A. and Pongor S.: Correlation of intrinsic DNA curvature with DNA
#' property periodicity, Febs. Lett. 393 (1996): 65-68.
#'
"scales"
