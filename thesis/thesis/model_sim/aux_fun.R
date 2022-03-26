
# mode version

find_des_mat <- function(M) {
  X <- M 
  for (i in 1:(NROW(M) - 1)) {
    X <- X %*% M + M
  }
  return(pmin(X, 1))
}

# Function: Is y a descendant of x? 
# Output: T or F

is_des <- function(M, y, x) {
  n <- ncol(M)
  if (n == 1) {
    return(FALSE)
  }
  X <- M
  for (i in 1:(n - 1)) {
    tmp <- X[y, x]
    if (tmp == 1) {
      return(TRUE)
    }
    X <- X %*% M
  }
  return(FALSE)
}



find_K_set = function(SCM, x, y) {
  mat_coef <- SCM$coef
  mat_dag  <- SCM$dag
  size <- length(SCM$noises)

  des_mat1 <- find_des_mat(mat_dag)
  m2 <- mat_dag
  m2[,x] <- 0
  des_mat2 <- find_des_mat(m2)

  vertices <- (1:size)[-c(x, y)]
  
  if (length(vertices) == 0) {
    return(0)
  }
  
  result <- rep(0, size)
  for (i in vertices) {
    for (j in x) {
      tmp <-  all(c(des_mat2[j, i], des_mat2[y, i], des_mat1[y, j]) == 1)
      if (tmp) {
        result[i] <- 1
        break
      }
    }
  }
  
  result <- (1:size)[as.logical(result)]
  if (length(result) == 0) {
    return(0)
  } else {
    return(result)
  }
}

path_walker <- function(dag, x, y, w = NULL, type = 0, 
                        xlast = NULL, silent = T, des_mat = NULL, 
                        printpath = F, path = NULL) {
  size <- NROW(dag)
  
  if (is.null(des_mat)) {
    des_mat <- find_des_mat(dag)
  }
  if (is.null(w)) {
    w <- c()
  }
  
  path <- c(path, x)
  xpar <- (1:size)[as.logical(dag[x, ])]
  xchi <- (1:size)[as.logical(dag[, x])]
  
  if (!silent) {
    cat("x =", x, "\n")
    cat("xpar =", xpar, "\n")
    cat("xchi =", xchi, "\n")
    cat("type =", type, "\n")
    cat("xlast =", xlast, "\n \n")
  }
  
  if (any(type == c(0, 1)) & length(xpar) > 0) {
    for (par in xpar) {
      if (!is.null(xlast)) {
        if (any(xlast == par)) {
          next
        }
      }
      if (any(par == y)) {
        if (printpath) {
          cat("Path =", path, par,". \n")
        }
        return(TRUE)
      }
      if (!is.null(w)) {
        if (any(par == w)) {
          next
        }
      }
      result <- path_walker(dag, par, y, w, type = 0, xlast = c(x, xlast), 
                            silent = silent, des_mat = des_mat, printpath = printpath, path = path)
      if (result) {
        return(TRUE)
      }
    }
  }
  if (any(type == c(0, 2)) & length(xchi) > 0) {
    for (chi in xchi) {
      if (!is.null(xlast)) {
        if (any(chi == xlast)) {
          next
        }
      }
      if (any(chi == y)) {
        if (printpath) {
          cat("Path =", path, chi,". \n")
        }
        return(TRUE)
      }
      if (any(chi == w)) {
        result <- path_walker(dag, chi, y, w, type = 1, xlast = c(x, xlast), 
                              silent = silent, des_mat = des_mat, printpath = printpath, path = path)
        if (result) {
          return(TRUE)
        }
      } else {
        if (any(des_mat[w, chi] == 1)) {
          result <- path_walker(dag, chi, y, w, type = 0, xlast = c(x, xlast), 
                                silent = silent, des_mat = des_mat, printpath = printpath, path = path)
          if (result) {
            return(TRUE)
          }
        } else {
          result <- path_walker(dag, chi, y, w, type = 2, xlast = c(x, xlast), 
                                silent = silent, des_mat = des_mat, printpath = printpath, path = path)
          if (result) {
            return(TRUE)
          }
        }
      }
    }
  }
  return(FALSE)
}

condi_indep <- function(SCM, x, y, w = NULL) {
  mat_coef <- SCM$coef
  mat_dag  <- SCM$dag
  
  no_path <- T
  nx <-length(x)
  i <- 1
  while(no_path & i <= nx) {
    xx <- x[i]
    have_path <- path_walker(mat_dag, xx, y, w)
    no_path <- !have_path
    i <- i + 1
  }
  return(no_path)
}

