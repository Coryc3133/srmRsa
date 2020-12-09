#' Check Social Relations Model
#'
#' This function concatenates and prints the string srm model. If you're familiar with
#' lavaan syntax, this can be used to make sure the SRM is specified correctly in the
#' event of unexpected results.
#' @param fitted_srm_rsa The fitted SRM RSA object from which you want the
#' printed SRM (in lavaan syntax).
#' @keywords social relations model
#' @return A printed, concatenated, string that corresponds to the SRM passed to lavaan.
#' @export

check_srm <- function(fitted_srm_rsa){
  return(cat(fitted_srm_rsa$srm_model))
}

#' Check Response Surface Analysis Paths
#'
#' This function concatenates and prints the string RSA Paths. If you're familiar with
#' lavaan syntax, this can be used to make sure the RSA is specified correctly in the event
#' of unexpected results.
#' @param fitted_srm_rsa The fitted SRM RSA object from which you want the
#' printed RSA Paths (in lavaan syntax).
#' @keywords RSA
#' @return A printed, concatenated, string that corresponds to the RSA Paths
#' that get combined with the SRM and passed to lavaan.
#' @export
check_rsa_paths <- function(fitted_srm_rsa){
  return(cat(fitted_srm_rsa$rsa_paths))
}

#' Check Social Relations Model Response Surface Analysis
#'
#' This function concatenates and prints the full SRM RSA model. If you're familiar with
#' lavaan syntax, this can be used to make sure the model is specified correctly in the event
#' of unexpected results.
#' @param fitted_srm_rsa The fitted SRM RSA object from which you want the
#' printed RSA Paths (in lavaan syntax).
#' @keywords RSA
#' @return A printed, concatenated, string that corresponds to the SRM RSA
#' that is passed to lavaan.
#' @export
check_srm_rsa <- function(fitted_srm_rsa){
  cat(fitted_srm_rsa$srm_rsa_model)
}

#' RSA (vs. NULL) Model Comparison Table
#'
#' This function extracts just the model comparison table from a fitted SRM
#' RSA object.
#' @param fitted_srm_rsa The fitted SRM RSA object from which you want the
#' model comparison table.
#' @param caption The caption for the table. Defaults to generic label.
#' @keywords RSA
#' @return a wrapper for knitr::kable that extracts and tables
#' the model comparison part of the SRM RSA output.
#'  Contains AIC, BIC, Chi-squared, difference in
#' model Chi-Squared values, and the corresponding df and p value
#' for the comparison.
#' @export
rsa_modelcomp_tbl <- function(fitted_srm_rsa,
                              caption =  "Comparing RSA to NULL (no Interaction or Polynomials) Model"){
  fitted_srm_rsa$srm_rsa_model_comp %>%
    knitr::kable(caption = caption,
                 digits = c(0, 2, 2, 2, 2, 0, 3))
}

#' RSA SRM Parameter Table
#'
#' This function extracts the SRM components, the slopes from the regression
#' used to define the response surface, the surface parameters (a1 through a5),
#' and the intercepts and variances of x and y variables from RSA.
#' @param fitted_srm_rsa The fitted SRM RSA object from which you want the
#' parameter table.
#' @param caption The caption for the table. Defaults to generic label.
#' @keywords RSA
#' @return a wrapper for knitr::kable that extracts and tables
#' the parameters from the SRMRSA model.
#' @export
#' @import knitr broom dplyr tidyr
srm_rsa_params_tbl <- function(fitted_srm_rsa,
                               caption = "SRM RSA Parameters"){
  design <- fitted_srm_rsa$model_info$design
  ratings <- fitted_srm_rsa$model_info$ratings
  if(ratings == "identical" &&
     design == "reciprocal" |
     ratings == "identical" &&
     design == "psxts"){
    tbl <- fitted_srm_rsa$srm_rsa_fit %>%
      broom::tidy(conf.int = TRUE) %>%
      dplyr::filter(label != "") %>%
      dplyr::distinct(label, .keep_all = TRUE) %>%
      dplyr::mutate(label = forcats::fct_relevel(label,
                                                 "targ_var",
                                                 "perc_var",
                                                 "rel_var",
                                                 "pt_cov",
                                                 "rel_cov",
                                                 "b1",
                                                 "b2",
                                                 "b3",
                                                 "b4",
                                                 "b5",
                                                 "a1",
                                                 "a2",
                                                 "a3",
                                                 "a4",
                                                 "a5",
                                                 "p10",
                                                 "p11",
                                                 "xy_var",
                                                 "xy_sq_var",
                                                 "xy_intx_var",
                                                 "xy_int",
                                                 "xy_sq_int",
                                                 "xy_intx_int"))
  }
  if(ratings == "different" &&
     design == "reciprocal" |
     ratings == "different" &&
     design == "psxts" |
     design %in% c("pxp", "pxps", "pxts")){
    tbl <- fitted_srm_rsa$srm_rsa_fit %>%
      broom::tidy(conf.int = TRUE) %>%
      dplyr::filter(label != "") %>%
      dplyr::distinct(label, .keep_all = TRUE) %>%
      dplyr::mutate(label = forcats::fct_relevel(label,
                                                 "targ_var",
                                                 "perc_var",
                                                 "rel_var",
                                                 "pt_cov",
                                                 "rel_cov",
                                                 "b1",
                                                 "b2",
                                                 "b3",
                                                 "b4",
                                                 "b5",
                                                 "a1",
                                                 "a2",
                                                 "a3",
                                                 "a4",
                                                 "a5",
                                                 "p10",
                                                 "p11",
                                                 "x_var",
                                                 "x_sq_var",
                                                 "y_var",
                                                 "y_sq_var",
                                                 "xy_intx_var",
                                                 "x_int",
                                                 "y_int",
                                                 "x_sq_int",
                                                 "y_sq_int",
                                                 "xy_intx_int"))
  }
  tbl %>%
    dplyr::arrange(label) %>%
    dplyr::select(label, estimate, p.value, conf.low, conf.high) %>%
    knitr::kable(caption = caption,
                 digits = c(NA, 2, 3, 2, 2))
}

#' SRM RSA Surface Plot
#'
#' This function extracts the beta coefficients and generates a
#' 3-d surface plot using the plotRSA() function from the RSA library
#' @param fitted_srm_rsa The fitted SRM RSA object from which you want the
#' surface plot.
#' @param ...  Optional additional arguments passed directly to \link[RSA]{plotRSA}
#' For example, it can be used to set X, Y, or Z labels, axis-limits, etc.
#' @keywords RSA
#' @return a wrapper for RSA::plotRSA() that generates 3-d surface plot based
#' on the SRMRSA model output.
#' @export
#' @import knitr broom dplyr tidyr RSA
srm_rsa_surface_plot <- function(fitted_srm_rsa, ...){
  betas <- fitted_srm_rsa$srm_rsa_fit %>%
    broom::tidy() %>%
    dplyr::filter(label != "") %>%
    dplyr::distinct(label, .keep_all = TRUE) %>%
    dplyr::filter(str_detect(label, "b")) %>%
    dplyr::pull(estimate)
RSA::plotRSA(x = betas[1], y = betas[2], x2 = betas[3], xy = betas[4], y2 = betas[5], ...)
}
