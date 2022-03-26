


measure_A_con <- function(SCM, x, y) {
  mat_coef <- SCM$coef
  noises   <- SCM$noises
  sing_tot <- calc_sing_tot(SCM)
  size <- length(SCM$noises)
  
  tot_y_x_ypar_noises <- mat_coef[y, x, drop = F] %*% sing_tot[x, , drop = F]
  
  var1 <- sum(tot_y_x_ypar_noises^2 * noises)
  var2 <- calc_var(SCM, y)
  return(var1 / var2)
}

measure_A_max <- function(SCM, x, y) {
  var1 <- calc_var(SCM, y = y, x = x)
  var2 <- calc_var(SCM, y = y, x = NULL)
  return(var1 / var2)
}

measure_A_rel <- function(SCM, x, y) {
  SCM2 <- SCM
  SCM2$coef[,x] <- 0
  SCM2$dag[,x] <- 0
  SCM2$noises[x] <- 0
  
  var1 <- calc_var(SCM2, y)
  var2 <- calc_var(SCM, y)
  return((var2 - var1) / var2)
}

measure_A_min <- function(SCM, x, y) {
  w <- find_K_set(SCM, x, y)
  if (any(w == 0)) {
    var1 <- calc_var(SCM, y = y, x = x)
  } else {
    var1 <- calc_var(SCM, y = y, x = x, w = w)
  }
  
  var2 <- calc_var(SCM, y)
  return(var1 / var2)
}

measure_A_adj <- function(SCM, x, y) {
  w <- find_K_set(SCM, x, y)
  
  if (any(w == 0)) {
    var1 <- calc_var(SCM, y = y, x = x)
    var2 <- calc_var(SCM, y = y)
  } else {
    SCM2 <- SCM
    SCM2$dag[, w] <- 0
    SCM2$coef[, w] <- 0
    return(measure_A_max(SCM2, x, y))
    # var1 <- calc_var(SCM, y = y, x = x, w = w)
    # var2 <- calc_var(SCM, y = y, x = y, w = w)
  }
  if(var1 / var2 >= 1) {
    ment(var1)
    ment(var2)
  }

  return(var1 / var2)
}

