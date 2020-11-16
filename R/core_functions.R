#' Build Social Relations Model Function
#'
#' This function generates a string object corresponding to a social relations model
#' in lavaan syntax. It requires a dataframe that contains columns for percever, target,
#' and group ID variables and a rating variable. It does not usually need to be called directly,
#' but is used in other functions.
#' @param data The dataframe. It must contain columns for percever, target,
#' and group ID variables and a rating variable. It should be in long format such
#' that each row is a rating by a perceiver for a given target in a given group.
#' @param perceiver_id A quoted string with the name of
#' the column containing perceiver IDs.
#' Perceiver IDs should be recycled across groups
#' (i.e., each group should have perceiver 1 to i
#' where i is the number of participants per group).
#' It can either be a number of a character string.
#' @param target_id A quoted string with the name of
#' the column containing target IDs.
#' target IDs should be recycled across groups
#' (i.e., each group should have perceiver 1 to i
#' where i is the number of participants per group).
#' It can either be a number of a character string.
#' @param group_id A quoted string with the name of
#' the column containing group IDs.It can either
#' be a number of a character string.
#' @param rating A quoted string with the name of
#' the column that contains ratings.
#' @keywords social relations model
#' @export
build_srm <-
  function(data,
           perceiver_id,
           target_id,
           group_id,
           rating) {
    p <- unique(data[,perceiver_id])
    t <- unique(data[,target_id])
    vars <- paste(c(rating), data[,perceiver_id], data[,target_id],  sep = "_")
    # Target Effects
    tfx <- ""
    for (i in t) {
      targ_vars <- vars[stringr::str_detect(vars,
                                            paste0("\\w*_\\d*_", i, "$"))]
      targ_vars <- unique(targ_vars)
      targ_vars <- targ_vars[!stringr::str_detect(targ_vars,
                                                  paste0("\\w*_", i, "_", i, "$"))]
      targ_vars <- stringr::str_flatten(targ_vars, " + 1*")
      tfx <- paste(tfx,
                   paste0("target",
                          i,
                          " =~  + 1*",
                          targ_vars),
                   sep = "\n \n")
    }
    # Perceiver Effects
    pfx <- ""
    for (i in p) {
      perc_vars <- vars[stringr::str_detect(vars,
                                            paste0("\\w*_", i, "_", "\\d"))]
      perc_vars <- unique(perc_vars)
      perc_vars <- perc_vars[!stringr::str_detect(perc_vars,
                                                  paste0("\\w*_", i, "_", i, "$"))]
      perc_vars <- stringr::str_flatten(perc_vars, " + 1*")
      pfx <- paste(pfx,
                   paste0("perceiver",
                          i,
                          " =~  + 1*",
                          perc_vars),
                   sep = "\n \n")
    }
    # group Effects
    gfx <- ""
    all_vars <- unique(vars)
    # remove self-reports
    for (i in 1:max(c(p, t))) {
      if (i == 1){
        all_vars_no_self <- all_vars[!stringr::str_detect(all_vars,
                                                          paste0("\\w*_", i, "_", i, "$"))]
      }
      else{all_vars_no_self <- all_vars_no_self[!stringr::str_detect(all_vars_no_self,
                                                                     paste0("\\w*_", i, "_", i, "$"))]}
    }
    grp_vars <- stringr::str_flatten(all_vars_no_self, " + 1*")
    gfx <- paste(gfx,
                 paste0("group",
                        " =~  + 1*",
                        grp_vars),
                 sep = "\n \n")
    # Relationship Variance
    dyads <- expand.grid(p = unique(data[, perceiver_id]),
                         t = unique(data[, target_id]))
    dyads <- dyads[which(dyads$t != dyads$p),]
    # create vector of relationship effects
    rfx <- stringr::str_flatten(
      paste0(
        paste(rating, dyads$p, dyads$t, sep = "_"),
        "~~ rel_var*",
        paste(rating, dyads$p, dyads$t, sep = "_")),
      "\n\n"
    )
    # Relationship (dyad) covariances
    rcov <- NULL
    for(i in 1:nrow(dyads)){
      rcov <-     c(rcov,
                    paste(paste(rating, dyads$p[i], dyads$t[i], sep = "_"),
                          "~~rel_cov*",
                          paste(rating, dyads$t[i], dyads$p[i], sep = "_")))
      dyads <- dyads[which(dyads$t != dyads$p[i] | dyads$p != dyads$t[i]),]
    }
    rcov <- rcov[-which(stringr::str_detect(rcov, "NA"))]
    rcov <- stringr::str_flatten(rcov, collapse = "\n\n")
    # intercepts
    ## relationship intercepts set to zero
    rel_ints <- stringr::str_flatten(paste(all_vars_no_self, "~0"), "\n\n")
    ## target intercepts set to zero
    targ_ints <- NULL
    for(i in t){
      targ_ints <- c(targ_ints,
                     paste0("target", i, "~0"))
    }
    targ_ints <- stringr::str_flatten(targ_ints, "\n\n")
    ## perceiver intercepts set to zero
    perc_ints <- NULL
    for(i in t){
      perc_ints <- c(perc_ints,
                     paste0("perceiver", i, "~0"))
    }
    perc_ints <- stringr::str_flatten(perc_ints, "\n\n")
    # Variances
    ## Target Variance
    targ_var <- NULL
    for(i in t){
      targ_var <- c(targ_var,
                    paste0("target", i, " ~~ ", "targ_var*target", i))
    }
    targ_var <- stringr::str_flatten(targ_var, "\n\n")
    ## Perceiver Variance
    perc_var <- NULL
    for(i in p){
      perc_var <- c(perc_var,
                    paste0("perceiver", i, " ~~ ", "perc_var*perceiver", i))
    }
    perc_var <- stringr::str_flatten(perc_var, "\n\n")
    # Perceiver-target Covariances
    pt_cov <- NULL
    for (i in 1:max(c(p, t))) {
      pt_cov <- c(pt_cov,
                  paste0("perceiver", i, " ~~ ", "pt_cov*target", i))
    }
    pt_cov <- stringr::str_flatten(pt_cov, "\n\n")
    # perceiver-group covariances (zero-out)
    pg_cov <- NULL
    for (i in p) {
      pg_cov <- c(pg_cov,
                  paste0("group", " ~~ ", "0*perceiver", i))
    }
    pg_cov <- stringr::str_flatten(pg_cov, "\n\n")
    # target-group covariances (zero-out)
    tg_cov <- NULL
    for (i in t) {
      tg_cov <- c(tg_cov,
                  paste0("group", " ~~ ", "0*target", i))
    }
    tg_cov <- stringr::str_flatten(tg_cov, "\n\n")
    # zero out actor with actor
    # NEED TO REMOVE DUPLICATE COVS
    targ0cov <- NULL
    t_t <- expand.grid(t1 = t, t2 = t)
    t_t <- t_t[which(t_t$t1 != t_t$t2),]
    for(i in 1:nrow(t_t)){
      targ0cov <- c(targ0cov,
                    paste0(paste0("target", t_t$t1[i]),
                           "~~0*",
                           paste0("target", t_t$t2[i])))
      t_t <- t_t[which(t_t$t2 != t_t$t1[i] | t_t$t1 != t_t$t2[i]),]
    }
    targ0cov <- targ0cov[-which(stringr::str_detect(targ0cov, "NA"))]
    targ0cov <- stringr::str_flatten(targ0cov, collapse = "\n\n")
    # zero out partner with partner
    # NEED TO REMOVE DUPLICATE COVS
    perc0cov <- NULL
    p_p <- expand.grid(p1 = p, p2 = p)
    p_p <- p_p[which(p_p$p1 != p_p$p2),]
    for(i in 1:nrow(p_p)){
      perc0cov <- c(perc0cov,
                    paste0(paste0("perceiver", p_p$p1[i]),
                           "~~0*",
                           paste0("perceiver", p_p$p2[i])))
      p_p <- p_p[which(p_p$p2 != p_p$p1[i] | p_p$p1 != p_p$p2[i]),]
    }
    perc0cov <- perc0cov[-which(stringr::str_detect(perc0cov, "NA"))]
    perc0cov <- stringr::str_flatten(perc0cov, collapse = "\n\n")
    # zero out actor with partner
    t_p <- expand.grid(tfx = t,
                       pfx = p)
    t_p <- t_p[which(t_p$tfx != t_p$pfx),]
    t_p$t <- paste0("target", t_p$tfx)
    t_p$p <- paste0("perceiver", t_p$pfx)
    t_p <- paste(t_p$t, "~~0*", t_p$p)
    t_p <- stringr::str_flatten(t_p, "\n\n")
    model_str <- paste(pfx, tfx, gfx, rfx,
                       rcov, rel_ints, targ_ints, perc_ints,
                       "group ~ 1", targ_var, perc_var,
                       pt_cov, pg_cov, tg_cov, t_p,
                       targ0cov, perc0cov,
                       sep = "\n\n")
    return(model_str)
  }

#' Build Response Surface Analysis Paths
#'
#' This function generates a string object corresponding
#' to the response surface analysis paths. It requires
#' a dataframe that contains columns for percever, target,
#' and group ID variables and rating variables for the X, Y,
#' Z variables for a response surface analysis (Z ~ X * Y).
#' It does not usually need to be called directly,
#' but is used in other functions.
#' @param data The dataframe. It must contain columns for percever, target,
#' and group ID variables and X, Y, and Z rating variables. Note that X and Y
#' can be the same variable. It should be in long format such
#' that each row is a rating by a perceiver for a given target in a given group.
#' @param perceiver_id A quoted string with the name of
#' the column containing perceiver IDs.
#' Perceiver IDs should be recycled across groups
#' (i.e., each group should have perceiver 1 to i
#' where i is the number of participants per group).
#' It can either be a number of a character string.
#' @param target_id A quoted string with the name of
#' the column containing target IDs.
#' target IDs should be recycled across groups
#' (i.e., each group should have perceiver 1 to i
#' where i is the number of participants per group).
#' It can either be a number of a character string.
#' @param group_id A quoted string with the name of
#' the column containing group IDs.It can either
#' be a number of a character string.
#' @param rating_x A quoted string with the name of
#' the column that contains ratings for the
#' x variable in the RSA.
#' @param rating_y A quoted string with the name of
#' the column that contains ratings for the
#' y variable in the RSA. Note that this can be the same
#' variable as x as long as design is not pxp.
#' @param rating_z A quoted string with the name of
#' the column that contains ratings for the
#' z variable (the outcome/DV) in the RSA.
#' @param design A quoted string specifying the design
#' of the RSA. Valid entries include:
#' \describe{
#' \item{reciprocal}{X and Y are reciprocal
#' ratings for each dyad; this can be on the same variable (e.g.,
#' A(B) Liking & B(A) Liking) or on different variables (e.g.,
#' A(B) Liking & B(A) Meta-Liking)}
#' \item{pxp}{X and Y are two ratings from the same
#' perceiver rating the same target (e.g., A(B) Liking
#' and A(B) Meta-Liking). These have to be different variables.}
#' \item{pxps}{X is a perceiver's rating of a target
#' and y is the perceivers' self-report.
#' This can be on the same or different variables.}
#' \item{pxts}{X is a perceiver's rating of a target
#' and y is the targets' self-report.
#' This can be on the same or different variables.}
#' \item{psxts}{X is a perceiver's self-report
#' and y is the targets' self-report.
#' This can be on the same or different variables.}}
#' @keywords social relations model, response surface analysis
#' @export
build_rsa_paths <- function(data,
                            perceiver_id,
                            target_id,
                            group_id,
                            rating_x,
                            rating_y,
                            rating_z,
                            design = NULL){

  p <- unique(data[,perceiver_id])
  t <- unique(data[,target_id])
  # vectorised function to order and combine values for dyad id
  f = function(x,y) paste(sort(c(x, y)), collapse="_")
  f = Vectorize(f)
  # Check if ratings are the same & what design is specified
  if(rating_x == rating_y &&
     is.null(design)){
    stop("X and Y variable are the same and no design is specified. Please specify the design using the design argument.")
  }
  # X Y identical variables, Reciprocal Perception (1_2) on X * Perception (2_1) on X
  if(rating_x == rating_y &&
     design == "reciprocal"){
    message("X and Y variable are identical and reciprocal (1_2 X 2_1) design specified")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(dyad_id = factor(f(p, t))) %>%
      dplyr::mutate(z = paste(rating_z, p, t, sep = "_"),
                   b1 = paste(rating_x, p, t, sep = "_"),
                   b2 = paste(rating_x, t, p, sep = "_"),
                   b3 = paste(rating_x, "sq", p, t, sep = "_"),
                   b4 = paste(rating_x, dyad_id, "intx", sep = "_"),
                   b5 = paste(rating_x, "sq", t, p, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }

    # Variances
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~~ ", "xy_var*" ,coef_var_mat[i,"b1"],
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~~ ", "xy_sq_var*", coef_var_mat[i,"b3"])
    }
    # variance for interaction term
    unique_xy <- unique(coef_var_mat[,"b4"])
    for(i in 1:nrow(unique_xy)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xy[i,"b4"], " ~~ ", "xy_int_var*" , unique_xy[i,"b4"])
    }
    # intercepts
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~ ", "xy_int*" ,   1,
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~ ", "xy_sq_int*", 1)
    }
    # intercepts forinteraction term
    for(i in 1:nrow(unique_xy)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xy[i,"b4"], " ~ ", "xy_int_int*" ,   1)
    }
  }
  # X Y different variables, Reciprocal Perception (1_2) on X * Perception (2_1) on Y
  if(rating_x != rating_y &&
     design == "reciprocal"){
    message("X and Y variable are different variables and pxp (perception X perception) design specified")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(z  = paste(rating_z, p, t, sep = "_"),
                    b1 = paste(rating_x, p, t, sep = "_"),
                    b2 = paste(rating_y, t, p, sep = "_"),
                    b3 = paste(rating_x, "sq", p, t, sep = "_"),
                    b4 = paste0(rating_x,"_", p, "_", t, "x", rating_y, t, "_", p),
                    b5 = paste(rating_y, "sq", t, p, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }
    # Variances
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~~ ", "x_var*" ,coef_var_mat[i,"b1"],
                             "\n\n",
                             coef_var_mat[i,"b2"], " ~~ ", "y_var*" ,coef_var_mat[i,"b2"],
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~~ ", "x_sq_var*", coef_var_mat[i,"b3"],
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~~ ", "xy_int_var*", coef_var_mat[i,"b4"],
                             "\n\n",
                             coef_var_mat[i,"b5"], " ~~ ", "y_sq_var*" , coef_var_mat[i,"b5"])
    }
    # intercepts
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~ ", "x_int*" ,   1,
                             "\n\n",
                             coef_var_mat[i,"b2"], " ~ ", "y_int*" , 1,
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~ ", "x_sq_int*", 1,
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~ ", "xy_int_int*",1,
                             "\n\n",
                             coef_var_mat[i,"b5"], " ~ ", "y_sq_int*" , 1)
    }
  }
  # X Y different variables, Perception (1_2) on X * Perception (1_2) on Y
  if(rating_x != rating_y &&
     design == "pxp"){
    message("X and Y variable are different variables and pxp (perception X perception) design specified")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(z  = paste(rating_z, p, t, sep = "_"),
                    b1 = paste(rating_x, p, t, sep = "_"),
                    b2 = paste(rating_y, p, t, sep = "_"),
                    b3 = paste(rating_x, "sq", p, t, sep = "_"),
                    b4 = paste0(rating_x,"_", p, "_", t, "x", rating_y, p, "_", t),
                    b5 = paste(rating_y, "sq", p, t, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }
    # Variances
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~~ ", "x_var*" ,coef_var_mat[i,"b1"],
                             "\n\n",
                             coef_var_mat[i,"b2"], " ~~ ", "y_var*" ,coef_var_mat[i,"b2"],
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~~ ", "x_sq_var*", coef_var_mat[i,"b3"],
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~~ ", "xy_int_var*", coef_var_mat[i,"b4"],
                             "\n\n",
                             coef_var_mat[i,"b5"], " ~~ ", "y_sq_var*" , coef_var_mat[i,"b5"])
    }
    # intercepts
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~ ", "x_int*" ,   1,
                             "\n\n",
                             coef_var_mat[i,"b2"], " ~ ", "y_int*" , 1,
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~ ", "x_sq_int*", 1,
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~ ", "xy_int_int*",1,
                             "\n\n",
                             coef_var_mat[i,"b5"], " ~ ", "y_sq_int*" , 1)
    }
  }
  # X Y identical; Perception (1_2) on X * perceiver Self-Report on X (1_1)
  if(rating_x == rating_y &&
     design == "pxps"){
    message("X and Y variable are identical and pxps (perception X perceiver self-perception) design specified;
            NOTE: this assumes that self-reports are coded as where perceiver_id and target_id are the same.")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(z = paste(rating_z, p, t, sep = "_"),
                    b1 = paste(rating_x, p, t, sep = "_"),
                    b2 = paste(rating_x, p, p, sep = "_"),
                    b3 = paste(rating_x, "sq", p, t, sep = "_"),
                    b4 = paste0(rating_x,"_", p, "_", t, "x", p, "_", p),
                    b5 = paste(rating_x, "sq", p, p, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }
    # Variances
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~~ ", "x_var*" ,coef_var_mat[i,"b1"],
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~~ ", "x_sq_var*", coef_var_mat[i,"b3"],
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~~ ", "xy_int_var*", coef_var_mat[i,"b4"])
    }
    # Variances for self-reports
    # special because they repeat across dyads
    unique_ys <- unique(coef_var_mat[,"b2"])
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~~ ", "y_var*" , unique_ys[i,"b2"])
    }
    # Variances for self-report squared terms
    # special because they repeat across dyads
    unique_ysq <- unique(coef_var_mat[,"b5"])
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~~ ", "y_sq_var*" , unique_ysq[i,"b5"])
    }
    # intercepts
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~ ", "x_int*" ,   1,
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~ ", "x_sq_int*", 1,
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~ ", "xy_int_int*",1)
    }
    # intercetps for self; special bc they repeat
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~ ", "y_int*" ,   1)
    }

    # intercetps for self squared term; special bc they repeat
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~ ", "y_sq_int*" ,   1)
    }
  }
  # X Y different; Perception (1_2) on X * perceiver Self-Report on Y (1_1)
  if(rating_x != rating_y &&
     design == "pxps"){
    message("X and Y variable are different and pxps (perception X perceiver self-perception) design specified;
            NOTE: this assumes that self-reports are coded as cases where perceiver_id and target_id are the same.")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(z  = paste(rating_z, p, t, sep = "_"),
                    b1 = paste(rating_x, p, t, sep = "_"),
                    b2 = paste(rating_y, p, p, sep = "_"),
                    b3 = paste(rating_x, "sq", p, t, sep = "_"),
                    b4 = paste0(rating_x,"_", p, "_", t, "x", rating_y, p, "_", p),
                    b5 = paste(rating_y, "sq", p, p, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }
    # Variances
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~~ ", "x_var*" ,coef_var_mat[i,"b1"],
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~~ ", "x_sq_var*", coef_var_mat[i,"b3"],
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~~ ", "xy_int_var*", coef_var_mat[i,"b4"])
    }
    # Variances for self-reports
    # special because they repeat across dyads
    unique_ys <- unique(coef_var_mat[,"b2"])
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~~ ", "y_var*" , unique_ys[i,"b2"])
    }
    # Variances for self-report squared terms
    # special because they repeat across dyads
    unique_ysq <- unique(coef_var_mat[,"b5"])
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~~ ", "y_sq_var*" , unique_ysq[i,"b5"])
    }
    # intercepts
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~ ", "x_int*" ,   1,
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~ ", "x_sq_int*", 1,
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~ ", "xy_int_int*",1)
    }
    # intercetps for self; special bc they repeat
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~ ", "y_int*" ,   1)
    }

    # intercetps for self squared term; special bc they repeat
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~ ", "y_sq_int*" ,   1)
    }
  }
  # X Y Identical Variable; Perception (1_2) on X * target Self-Report on X (2_2)
  if(rating_x == rating_y &&
     design == "pxts"){
    message("X and Y variable are identical and pxts (perception X target self-perception) design specified;
            NOTE: this assumes that self-reports are coded as where perceiver_id and target_id are the same.")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(z  = paste(rating_z, p, t, sep = "_"),
                    b1 = paste(rating_x, p, t, sep = "_"),
                    b2 = paste(rating_x, t, t, sep = "_"),
                    b3 = paste(rating_x, "sq", p, t, sep = "_"),
                    b4 = paste0(rating_x,"_", p, "_", t, "x", t, "_", t),
                    b5 = paste(rating_x, "sq", t, t, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }
    # Variances
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~~ ", "x_var*" ,coef_var_mat[i,"b1"],
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~~ ", "x_sq_var*", coef_var_mat[i,"b3"],
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~~ ", "xy_int_var*", coef_var_mat[i,"b4"])
    }
    # Variances for self-reports
    # special because they repeat across dyads
    unique_ys <- unique(coef_var_mat[,"b2"])
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~~ ", "y_var*" , unique_ys[i,"b2"])
    }
    # Variances for self-report squared terms
    # special because they repeat across dyads
    unique_ysq <- unique(coef_var_mat[,"b5"])
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~~ ", "y_sq_var*" , unique_ysq[i,"b5"])
    }
    # intercepts
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~ ", "x_int*" ,   1,
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~ ", "x_sq_int*", 1,
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~ ", "xy_int_int*",1)
    }
    # intercetps for self; special bc they repeat
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~ ", "y_int*" ,   1)
    }

    # intercetps for self squared term; special bc they repeat
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~ ", "y_sq_int*" ,   1)
    }
  }
  # X Y are different; Perception (1_2) on X * target Self-Report on Y (2_2)
  if(rating_x != rating_y &&
     design == "pxts"){
    message("X and Y variable are different variables and pxts (perception X target self-perception) design specified;
            NOTE: this assumes that self-reports are coded as where perceiver_id and target_id are the same.")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(z  = paste(rating_z, p, t, sep = "_"),
                    b1 = paste(rating_x, p, t, sep = "_"),
                    b2 = paste(rating_y, t, t, sep = "_"),
                    b3 = paste(rating_x, "sq", p, t, sep = "_"),
                    b4 = paste0(rating_x,"_", p, "_", t, "x", rating_y, t, "_", t),
                    b5 = paste(rating_y, "sq", t, t, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }
    # Variances
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~~ ", "x_var*" ,coef_var_mat[i,"b1"],
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~~ ", "x_sq_var*", coef_var_mat[i,"b3"],
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~~ ", "xy_int_var*", coef_var_mat[i,"b4"])
    }
    # Variances for self-reports
    # special because they repeat across dyads
    unique_ys <- unique(coef_var_mat[,"b2"])
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~~ ", "y_var*" , unique_ys[i,"b2"])
    }
    # Variances for self-report squared terms
    # special because they repeat across dyads
    unique_ysq <- unique(coef_var_mat[,"b5"])
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~~ ", "y_sq_var*" , unique_ysq[i,"b5"])
    }
    # intercepts
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"b1"], " ~ ", "x_int*" ,   1,
                             "\n\n",
                             coef_var_mat[i,"b3"], " ~ ", "x_sq_int*", 1,
                             "\n\n",
                             coef_var_mat[i,"b4"], " ~ ", "xy_int_int*",1)
    }
    # intercetps for self; special bc they repeat
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~ ", "y_int*" ,   1)
    }

    # intercetps for self squared term; special bc they repeat
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~ ", "y_sq_int*" ,   1)
    }
  }
  # X Y Identical Variable; perceiver self-report (1_1) on X * target Self-Report on X (2_2)
  if(rating_x == rating_y &&
     design == "psxts"){
    message("X and Y variable are identical and psxts (perceiver self-report X target self-report) design specified;
            NOTE: this assumes that self-reports are coded as where perceiver_id and target_id are the same.")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(dyad_id = factor(f(p, t))) %>%
      dplyr::mutate(z  = paste(rating_z, p, t, sep = "_"),
                    b1 = paste(rating_x, p, p, sep = "_"),
                    b2 = paste(rating_x, t, t, sep = "_"),
                    b3 = paste(rating_x, "sq", p, p, sep = "_"),
                    b4 = paste(rating_x, dyad_id, "intx", sep = "_"),
                    b5 = paste(rating_x, "sq", t, t, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }
    # Variances
    # Variances for perceivers' self-reports
    unique_xs <- unique(coef_var_mat[,"b1"])
    for(i in 1:nrow(unique_xs)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xs[i,"b1"], " ~~ ", "x_var*" , unique_xs[i,"b1"])
    }
    # Variances for targets' self-reports
    unique_ys <- anti_join(unique(coef_var_mat[,"b2"]),unique(coef_var_mat[,"b1"]),
                           by = c("b2" = "b1"))
    if(nrow(unique_ys > 0)){
      for(i in 1:nrow(unique_ys)){
        coef_var_str <- paste0(coef_var_str, "\n\n",
                               unique_ys[i,"b2"], " ~~ ", "y_var*" , unique_ys[i,"b2"])
      }
    }
    # Variances for targets' self-reports
    unique_xsq <- unique(coef_var_mat[,"b3"])
    for(i in 1:nrow(unique_xsq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xsq[i,"b3"], " ~~ ", "x_sq_var*" , unique_xsq[i,"b3"])
    }
    # variance for interaction term
    unique_xy <- unique(coef_var_mat[,"b4"])
    for(i in 1:nrow(unique_xy)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xy[i,"b4"], " ~~ ", "xy_int_var*" , unique_xy[i,"b4"])
    }
    # Variances for targets' self-report squared terms
    unique_ysq <- anti_join(unique(coef_var_mat[,"b5"]),unique(coef_var_mat[,"b3"]),
                           by = c("b5" = "b3"))
    if(nrow(unique_ysq > 0)){
      for(i in 1:nrow(unique_ysq)){
        coef_var_str <- paste0(coef_var_str, "\n\n",
                               unique_ysq[i,"b5"], " ~~ ", "y_sq_var*" , unique_ysq[i,"b5"])
      }
    }
    # intercepts
    # intercepts for targets' & perceivers' self
    for(i in 1:nrow(unique_xs)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xs[i,"b1"], " ~ ", "x_int*" ,   1)
    }
    unique_ys <- anti_join(unique(coef_var_mat[,"b2"]),unique(coef_var_mat[,"b1"]),
                           by = c("b2" = "b1"))
    if(nrow(unique_ys > 0)){
      for(i in 1:nrow(unique_ys)){
        coef_var_str <- paste0(coef_var_str, "\n\n",
                               unique_ys[i,"b2"], " ~~ ", "y_int*" , 1)
      }
    }
    # intercepts for perceivers' self Squared term
    for(i in 1:nrow(unique_xsq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xsq[i,"b3"], " ~ ", "x_sq_int*" ,   1)
    }
    # intercetps for perceivers' X Targets' self interaction term
    for(i in 1:nrow(unique_xy)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xy[i,"b4"], " ~ ", "xy_int_int*" ,   1)
    }
    unique_ysq <- anti_join(unique(coef_var_mat[,"b5"]),unique(coef_var_mat[,"b3"]),
                           by = c("b5" = "b3"))
    if(nrow(unique_ysq > 0)){
      for(i in 1:nrow(unique_ysq)){
        coef_var_str <- paste0(coef_var_str, "\n\n",
                               unique_ysq[i,"b5"], " ~~ ", "y_sq_int*" , 1)
      }
    }
  }
  # X Y are different; perceiver self-report (1_1) on X * target Self-Report on Y (2_2)
  if(rating_x != rating_y &&
     design == "psxts"){
    message("X and Y variable are different variables and psxts (perceiver self-report X target self-perception) design specified;
            NOTE: this assumes that self-reports are coded as where perceiver_id and target_id are the same.")
    # get p_t matrix
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    coef_var_mat <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(z  = paste(rating_z, p, t, sep = "_"),
                    b1 = paste(rating_x, p, p, sep = "_"),
                    b2 = paste(rating_y, t, t, sep = "_"),
                    b3 = paste(rating_x, "sq", p, p, sep = "_"),
                    b4 = paste0(rating_x,"_", p, "_", p, "x", rating_y, t, "_", t),
                    b5 = paste(rating_y, "sq", t, t, sep = "_"))
    coef_var_str <- ""
    # regression paths
    for(i in 1:nrow(coef_var_mat)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             coef_var_mat[i,"z"], " ~ ",
                             "b1*", coef_var_mat[i,"b1"], " + ",
                             "b2*", coef_var_mat[i,"b2"], " + ",
                             "b3*", coef_var_mat[i,"b3"], " + ",
                             "b4*", coef_var_mat[i, "b4"], " + ",
                             "b5*", coef_var_mat[i, "b5"], "\n\n")
    }
    # Variances
    # Variances for targets' self-reports
    unique_xs <- unique(coef_var_mat[,"b1"])
    for(i in 1:nrow(unique_xs)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xs[i,"b1"], " ~~ ", "x_var*" , unique_xs[i,"b1"])
    }
    # Variances for targets' self-reports
    unique_ys <- unique(coef_var_mat[,"b2"])
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~~ ", "y_var*" , unique_ys[i,"b2"])
    }
    # Variances for targets' self-reports
    unique_xsq <- unique(coef_var_mat[,"b3"])
    for(i in 1:nrow(unique_xsq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xsq[i,"b3"], " ~~ ", "x_sq_var*" , unique_xsq[i,"b3"])
    }
    # variance for interaction term
    unique_xy <- unique(coef_var_mat[,"b4"])
    for(i in 1:nrow(unique_xy)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xy[i,"b4"], " ~~ ", "xy_int_var*" , unique_xy[i,"b4"])
    }
    # Variances for targets' self-report squared terms
    unique_ysq <- unique(coef_var_mat[,"b5"])
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~~ ", "y_sq_var*" , unique_ysq[i,"b5"])
    }
    # intercepts
    # intercepts for perceivers' self
    for(i in 1:nrow(unique_xs)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xs[i,"b1"], " ~ ", "x_int*" ,   1)
    }
    # intercepts for targets' self
    for(i in 1:nrow(unique_ys)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ys[i,"b2"], " ~ ", "y_int*" ,   1)
    }
    # intercepts for perceivers' self Squared term
    for(i in 1:nrow(unique_xsq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xsq[i,"b3"], " ~ ", "x_sq_int*" ,   1)
    }
    # intercetps for perceivers' X Targets' self interaction term
    for(i in 1:nrow(unique_xy)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_xy[i,"b4"], " ~ ", "xy_int_int*" ,   1)
    }
    # intercepts for targets' self squared term
    for(i in 1:nrow(unique_ysq)){
      coef_var_str <- paste0(coef_var_str, "\n\n",
                             unique_ysq[i,"b5"], " ~ ", "y_sq_int*" ,   1)
    }
  }
  # Check if design is not available
  if(!is.null(design) &&
     !(design %in% c("reciprocal", "pxp", "pxps", "pxts", "psxts"))){
       stop("design is mis-specified or undefined; please choose one of the defined options:
            reciprocal, pxp, pxps, pxts")
     }
  # Specify surface params
  coef_var_str <- paste(coef_var_str,
                        "a1 := b1 + b2
                         a2 := b3 + b4 + b5
                         a3 := b1 - b2
                         a4 := b3 - b4 + b5
                         a5 := b3 - b5", sep = "\n\n")
  # specify return
  return(coef_var_str)
}

#' Fit Response Surface Analysis Paths
#'
#' This function fits the social relations models response surface analysis.
#' It requires a dataframe that contains columns for percever, target,
#' and group ID variables and rating variables for the X, Y,
#' Z variables for a response surface analysis (Z ~ X * Y).
#' The Z variable will be treated as the rating in the social relations model,
#' and the response surface analysis will be conducted on the relationship effect
#' of that rating (i.e., after partialling out target, perceiver, and group effects).
#'
#' This function will mean-center the X and Y variables for you before it runs the analyses.
#' It also creates the cross-products for you based on the rating variables and design specified.
#' It will then fit 3 models:
#' \describe{
#' \item{Social Relations Model}{This is a basic Social Relations model
#' Estimating Target, Perceiver, Relationship effects,
#' and the Target-Perceiver and relationship Covariances}
#' \item{SRM RSA Null Model}{This is the SRM RSA model where slopes for
#' the interation and polynomial terms are set to zero,
#' to test whether and RSA is appropriate.}
#' \item{SRM RSA}{This is the full SRM RSA model with b1 through b5 and
#' surface parameters a1 through a5}}
#' For each model, it returns a string for the model (as lavaan syntax),
#' fitted lavaan models, and a model comparison table comparing model 2 to 3.
#'
#' @param data The dataframe. It must contain columns for percever, target,
#' and group ID variables and X, Y, and Z rating variables. Note that X and Y
#' can be the same variable. It should be in long format such
#' that each row is a rating by a perceiver for a given target in a given group.
#' @param perceiver_id A quoted string with the name of
#' the column containing perceiver IDs.
#' Perceiver IDs should be recycled across groups
#' (i.e., each group should have perceiver 1 to i
#' where i is the number of participants per group).
#' Perceiver and Target id should match such that
#' perceiver i is target i and vice versa.
#' It can either be a number of a character string.
#' @param target_id A quoted string with the name of
#' the column containing target IDs.
#' target IDs should be recycled across groups
#' (i.e., each group should have target 1 to i
#' where i is the number of participants per group).
#' Perceiver and Target id should match such that
#' perceiver i is target i and vice versa.
#' It can either be a number of a character string.
#' @param group_id A quoted string with the name of
#' the column containing group IDs.It can either
#' be a number of a character string.
#' @param rating_x A quoted string with the name of
#' the column that contains ratings for the
#' x variable in the RSA.
#' @param rating_y A quoted string with the name of
#' the column that contains ratings for the
#' y variable in the RSA. Note that this can be the same
#' variable as x as long as design is not pxp.
#' @param rating_z A quoted string with the name of
#' the column that contains ratings for the
#' z variable (the outcome/DV) in the RSA.
#' @param design A quoted string specifying the design
#' of the RSA. Valid entries include:
#' \describe{
#' \item{reciprocal}{X and Y are reciprocal
#' ratings for each dyad; this can be on the same variable (e.g.,
#' A(B) Liking & B(A) Liking) or on different variables (e.g.,
#' A(B) Liking & B(A) Meta-Liking)}
#' \item{pxp}{X and Y are two ratings from the same
#' perceiver rating the same target (e.g., A(B) Liking
#' and A(B) Meta-Liking). These have to be different variables.}
#' \item{pxps}{X is a perceiver's rating of a target
#' and y is the perceivers' self-report.
#' This can be on the same or different variables.}
#' \item{pxts}{X is a perceiver's rating of a target
#' and y is the targets' self-report.
#' This can be on the same or different variables.}
#' \item{psxts}{X is a perceiver's self-report
#' and y is the targets' self-report.
#' This can be on the same or different variables.}}
#' @param ... Optional additional arguments passed directly to
#' \link[lavaan]{lavaan} as it fits each model. For example, it can be used
#' to specify bootstrapped or robust standard errors. Note this will affect
#' all three of the fitted models.
#' @keywords social relations model, response surface analysis
#' @export
#' @import RSA lavaan rlang stringr dplyr tidyr
#' @examples
#' #NEED TO SIMULATE DATA
#' @return The function returns a list containing the following elements:
#' \describe{
#' \item{srm_model}{string. SRM model in lavaan syntax treating z as rating.}
#' \item{rsa_paths}{string. regression paths, intercepts, and variances for rsa in
#' lavaan syntax.}
#' \item{srm_rsa_model}{string. SRM RSA model in lavaan syntax. }
#' \item{srm_fit}{fitted lavaan model. Contains basic SRM (on Z variable).}
#' \item{srm_rsa_null_fit}{fitted lavaan model. Contains SRM RSA Null model
#' where interaction and polynomial are set to zero. Preimarily
#' used for model comparison.}
#' \item{srm_rsa_fit}{fitted lavaan model. Contains full SRM RSA
#' Model including surface parameters.}
#' \item{srm_rsa_model_comp}{a tibble. Contains the model comparison results
#' from comparing the SRM RSA Null model to the full SRM RSA model.}
#' \item{wide_df}{a tibble. Contains the wide version of the data
#' created by the fit_srm_rsa function. It should be a row for each group, and a column
#' for every rating in the format rating_perceiverid_targetid;
#' squared terms are in the format rating_sq_perceiverid_targetid;
#' interaction terms when x and y are the same variable
#' (with different perceiver-target combos) are in the format
#' rating_perceiverid_targetid_x_perceiverid_targetid; interaction
#' terms when x and y are different variables are in the format
#' ratingx_perceiverid_targetid_x_ratingy_perceiverid_targetid.}}

fit_srm_rsa <- function(data,
                        perceiver_id,
                        target_id,
                        group_id,
                        rating_x,
                        rating_y,
                        rating_z,
                        design = NULL,
                        ...){
  p <- unique(data[,perceiver_id])
  t <- unique(data[,target_id])

  # subset data to have just the vars of interest
  data <- as.data.frame(data)
  data <- data[, c(perceiver_id,
                   target_id,
                   group_id,
                   rating_x,
                   rating_y,
                   rating_z)]
  # check Mean Centerin and Mean-center if necessary
  ## Mean Center Function
  mean_center <- function(vec){
    mc_vec <- vec - mean(vec, na.rm = TRUE)
    return(mc_vec)

  }
  if(round(mean(data[,rating_x], na.rm = TRUE), 3) != 0){
    data[,rating_x]<- mean_center(data[,rating_x])
    message("X was not mean-centered; grand mean-centering X")
  }
  if(rating_x == rating_y){message("Y is the same variable as X; it was thus also grand mean centered.")}
  if(round(mean(data[,rating_y], na.rm = TRUE), 3) != 0){
    data[,rating_y]<- mean_center(data[,rating_y])
    message("Y was not mean-centered; grand mean-centering Y")
  }
  # calculate Sq terms
  ## Calculate Sq Term Function
  sq_term <- function(vec){
    sq_vec <- vec^2
    return(sq_vec)
  }
  ## calculate sq term for x
  data[, paste(rating_x, "sq", sep = "_")] <- sq_term(data[,rating_x])
  ## calculate for y if its a different variable
  if(rating_x != rating_y){
    data[, paste(rating_y, "sq", sep = "_")] <- sq_term(data[,rating_y])
  }
  # spread data
  wide_data <- data %>%
    tidyr::gather(., var, rating,
                  -match(group_id, names(.)),
                  -match(perceiver_id, names(.)),
                  -match(target_id, names(.))) %>%
    tidyr::unite(.,
                 var, var,
                 match(perceiver_id, names(.)),
                 match(target_id, names(.))) %>%
    tidyr::spread(var, rating)
  # Calculate Cross-Products
  ## Reciprocal Perception X Perception 1_2 X 2_1 on same variable
  # vectorised function to order and combine values for dyad id
  f = function(x,y) paste(sort(c(x, y)), collapse="_")
  f = Vectorize(f)

  if(rating_x == rating_y &&
     design == "reciprocal"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    cross_prods <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(dyad_id = factor(f(p, t))) %>%
      dplyr::mutate(b4 = paste(rating_x, dyad_id, "intx", sep = "_")) %>%
      distinct(dyad_id, .keep_all = TRUE)
    for(i in 1:nrow(cross_prods)){
      cross_prod <- paste(rating_x, cross_prods$dyad_id[i], "intx", sep = "_")
      cross_x <-  paste0(rating_x,"_", cross_prods$p[i], "_", cross_prods$t[i])
      cross_y <- paste0(rating_x,"_", cross_prods$t[i], "_", cross_prods$p[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }

  ## Reciprocal Perception X Perception 1_2 X 2_1 on different XY variables
  if(rating_x != rating_y &&
     design == "reciprocal"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    for(i in 1:nrow(p_t)){
      cross_prod <- paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i], "x", rating_y, p_t$t[i], "_", p_t$p[i])
      cross_x <-  paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i])
      cross_y <- paste0(rating_y,"_", p_t$t[i], "_", p_t$p[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }
  ## Perception X Perception 1_2 X 1_2 on different XY variables
  if(rating_x != rating_y &&
     design == "pxp"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    for(i in 1:nrow(p_t)){
      cross_prod <- paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i], "x", rating_y, p_t$p[i], "_", p_t$t[i])
      cross_x <-  paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i])
      cross_y <- paste0(rating_y,"_", p_t$p[i], "_", p_t$t[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }
  # Perception X Perceiver Self-Perception 1_2 X 1_1 on Same Variable
  if(rating_x == rating_y &&
     design == "pxps"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    for(i in 1:nrow(p_t)){
      cross_prod <- paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i], "x", p_t$p[i], "_", p_t$p[i])
      cross_x <-  paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i])
      cross_y <- paste0(rating_x,"_", p_t$p[i], "_", p_t$p[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }

  # Perception X Perceiver Self-Perception 1_2 X 1_1 on Different XY Variable
  if(rating_x != rating_y &&
     design == "pxps"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    for(i in 1:nrow(p_t)){
      cross_prod <- paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i], "x", rating_y, p_t$p[i], "_", p_t$p[i])
      cross_x <-  paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i])
      cross_y <- paste0(rating_y,"_", p_t$p[i], "_", p_t$p[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }

  # Perception X Target Self-Perception 1_2 X 2_2 on Same Variable
  if(rating_x == rating_y &&
     design == "pxts"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    for(i in 1:nrow(p_t)){
      cross_prod <- paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i], "x", p_t$t[i], "_", p_t$t[i])
      cross_x <-  paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i])
      cross_y <- paste0(rating_x,"_", p_t$t[i], "_", p_t$t[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }

  # Perception X Target Self-Perception 1_2 X 2_2 on Different XY Variable
  if(rating_x != rating_y &&
     design == "pxts"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    for(i in 1:nrow(p_t)){
      cross_prod <- paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i], "x", rating_y, p_t$t[i], "_", p_t$t[i])
      cross_x <-  paste0(rating_x,"_", p_t$p[i], "_", p_t$t[i])
      cross_y <- paste0(rating_y,"_", p_t$t[i], "_", p_t$t[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }

  # Perceiver self-perception X Target Self-Perception 1_1 X 2_2 on Same Variable
  if(rating_x == rating_y &&
     design == "psxts"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    cross_prods <- p_t %>%
      dplyr::as_tibble() %>%
      dplyr::mutate(dyad_id = factor(f(p, t))) %>%
      dplyr::mutate(b4 = paste(rating_x, dyad_id, "int", sep = "_")) %>%
      distinct(dyad_id, .keep_all = TRUE)
    for(i in 1:nrow(cross_prods)){
      cross_prod <- paste(rating_x, cross_prods$dyad_id[i], "intx", sep = "_")
      cross_x <-  paste0(rating_x,"_", cross_prods$p[i], "_", cross_prods$t[i])
      cross_y <- paste0(rating_x,"_", cross_prods$t[i], "_", cross_prods$p[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }

  # Perceiver self-perception X Target Self-Perception 1_2 X 2_2 on Different XY Variable
  if(rating_x != rating_y &&
     design == "psxts"){
    p_t <- expand.grid(p = p, t = t)
    p_t <- p_t[p_t$p != p_t$t,]
    row.names(p_t) <- 1:nrow(p_t)
    for(i in 1:nrow(p_t)){
      cross_prod <- paste0(rating_x,"_", p_t$p[i], "_", p_t$p[i], "x", rating_y, p_t$t[i], "_", p_t$t[i])
      cross_x <-  paste0(rating_x,"_", p_t$p[i], "_", p_t$p[i])
      cross_y <- paste0(rating_y,"_", p_t$t[i], "_", p_t$t[i])
      wide_data[, cross_prod] <-  wide_data[, cross_x] * wide_data[, cross_y]
    }
  }


  srm <- build_srm(data = data,
                   perceiver_id = perceiver_id,
                   target_id = target_id,
                   group_id = group_id,
                   rating = rating_z)
  rsa_paths <- build_rsa_paths(data = data,
                               perceiver_id = perceiver_id,
                               target_id = target_id,
                               group_id = group_id,
                               rating_x = rating_x,
                               rating_y = rating_y,
                               rating_z = rating_z,
                               design = design)
  srm_rsa_model <- paste(srm, rsa_paths, sep = "\n\n")
  srm_rsa_null_model <- stringr::str_replace_all(srm_rsa_model, "b3\\*|b4\\*|b5\\*", "0*")
  srm_rsa_null_model <- stringr::str_remove_all(srm_rsa_null_model, "a1.*|a2.*|a3.*|a4.*|a5.*|")
  # fit models
  basic_srm_fit <- lavaan::lavaan(srm,
                                  data = wide_data,
                                  missing = "fiml",
                                  ...)
  srm_rsa_null_fit <- lavaan::lavaan(srm_rsa_null_model,
                                     data = wide_data,
                                     missing = "fiml",
                                     ...)
  srm_rsa_fit <- lavaan::lavaan(srm_rsa_model,
                                data = wide_data,
                                missing = "fiml",
                                ...)
  srm_rsa_model_comp <- lavaan::anova(srm_rsa_null_fit,
                                      srm_rsa_fit) %>%
    dplyr::as_tibble()
  wide_data <- dplyr::as_tibble(wide_data)
  return(list(
    srm_model = srm,
    rsa_paths = rsa_paths,
    srm_rsa_model = srm_rsa_model,
    srm_fit =  basic_srm_fit,
    srm_rsa_null_fit = srm_rsa_null_fit,
    srm_rsa_fit = srm_rsa_fit,
    srm_rsa_model_comp = srm_rsa_model_comp,
    wide_df = wide_data))
}
