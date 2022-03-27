#
#
#
#
# This file contains function to calculate total effects.
#
#
#
#
#

# Function to calculate all singular total effect between all variables.
# Output: Matrix with entries equal to the columns total effect on the row.
#         The diagonal is set to 1's.
calc_sing_tot <- function(SCM) {
  mat_coef <- SCM$coef
  n        <- NROW(mat_coef)
  result   <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:i) {
      if (j < i) {
        result[i, j] <- sum(mat_coef[i, j:(i - 1)] * result[j:(i - 1), j])
      }
      if (j == i) {
        result[i, j] <- 1
      }
    }
  }
  result
}

# Function to calculate variance. See next function
# Output: If only y: variance of y
#         If x is given: variance of the effect of x on y
#         If w is given: variance of the effect of x on y after removal of w

calc_var_coef <- function(SCM, y, x = NULL, w = NULL) {
  mat_coef <- SCM$coef
  noises   <- SCM$noises
  size     <- NROW(mat_coef)
  sing_tot <- calc_sing_tot(SCM)
  
  if (!is.null(x)) {        
    n <- length(x)
    tmp1 <- sing_tot[x, x, drop = F] - diag(n)
    tmp2 <- diag(n)
    if (n > 1) {
      for (i in 1:(n - 1)) {
        tmp2 <- diag(n) - tmp2 %*% tmp1
      }
    }
    totef_y_x <- sing_tot[y, x] %*% tmp2
    if (!is.null(w)) {
      m <- length(w)
      tmp3 <- sing_tot[w, w, drop = F] - diag(m)
      tmp4 <- diag(m)
      if (m > 1) {
        for (i in 1:(m - 1)) {
          tmp4 <- diag(m) - tmp4 %*% tmp3
        }
      }
      totef_x_w <- sing_tot[x, w] %*% tmp4
      return(totef_y_x %*% (sing_tot[x,] - totef_x_w %*% sing_tot[w, ]))
    }
    return(totef_y_x %*% sing_tot[x,])
  }
  return(sing_tot[y, , drop = F])
}

# Function. Wrapper of calc_var_coef to output actual variance.
calc_var <- function(SCM, y, x = NULL, w = NULL) {
  noises <- SCM$noises
  return(sum(calc_var_coef(SCM, y = y, x = x, w = w)^2 * noises))
}

# Function to calculate simple latent projection of SCM over W. See thesis
# Checks first if the projection is simple.
# Outputs: Simple latent projection of the SCM.
#          For simplicity we swap coefficients and noises on L to 0's insead of returning a smaller model
find_lat_proj <- function(SCM, W) {
  mat_coef <- SCM$coef
  mat_dag  <- SCM$dag
  noises   <- SCM$noises
  size     <- length(noises)
  
  if (length(W) == size) {
    return(SCM)
  }
  
  is_simple <- is_simple_lat_proj(SCM, W)
  if (!is_simple) {
    return(c(TRUE, FALSE))
  }
  
  SCM2 <- SCM
  SCM2$coef[, W] <- 0
  SCM2$dag[, W]  <- 0
  sing_tot <- calc_sing_tot(SCM2)
  des_mat  <- find_des_mat(SCM2$dag)

  noises[W]      <- noises[W]          + sing_tot[W, -W, drop = F]^2 %*% noises[-W]
  mat_dag[W, W]  <- pmin(SCM$dag[W, W] + des_mat[W, -W, drop = F]    %*% mat_dag[-W, W, drop = F], 1)
  mat_coef[W, W] <- SCM$coef[W, W]     + sing_tot[W, -W, drop = F]   %*% mat_coef[-W, W, drop = F]
  
  noises[-W]     <- 0
  mat_dag[-W, ]  <- 0
  mat_dag[, -W]  <- 0
  mat_coef[-W, ] <- 0
  mat_coef[, -W] <- 0
  
  return(list(
    dag    = mat_dag,
    coef   = mat_coef,
    noises = noises
  ))
}
