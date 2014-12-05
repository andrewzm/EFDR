###########################################################################
###########################################################################
###########################################################################

dwt.2d <- function(x, wf, J=4, boundary="periodic")
{
  m <- dim(x)[1]
  storage.mode(m) <- "integer"
  n <- dim(x)[2]
  storage.mode(n) <- "integer"

  dict <- wave.filter(wf)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  z <- matrix(0, m/2, n/2)
  storage.mode(z) <- "double"

  x.wt <- vector("list", 3*J+1)
  x.names <- NULL
  for(j in 1:J) {
    out <- .C("two_D_dwt", "Image"=as.double(x), "Rows"=m, "Cols"=n, 
                "filter.length"=L, "hpf"=h, "lpf"=g, "LL"=z, "LH"=z,
                "HL"=z, "HH"=z, PACKAGE="EFDR")[7:10]
    if(j < J) {
      index <- (3*j-2):(3*j)
      x.wt[index] <- out[-1]
      x.names <- c(x.names, sapply(names(out)[-1], paste, j, sep=""))
      x <- out[[1]]
      m <- dim(x)[1]
      storage.mode(m) <- "integer"
      n <- dim(x)[2]
      storage.mode(n) <- "integer"
      z <- matrix(0, m/2, n/2)
      storage.mode(z) <- "double"
    }
    else {
      index <- (3*j):(3*(j+1)) - 2
      x.wt[index] <- out[c(2:4,1)]
      x.names <- c(x.names, sapply(names(out)[c(2:4,1)], paste, j, sep=""))
    }
  }

  names(x.wt) <- x.names
  attr(x.wt, "J") <- J
  attr(x.wt, "wavelet") <- wf
  attr(x.wt, "boundary") <- boundary
  attr(x.wt, "class") <- "dwt.2d"
  x.wt
}

###########################################################################
###########################################################################
###########################################################################

idwt.2d <- function(y)
{
  J <- attributes(y)$J

  dict <- wave.filter(attributes(y)$wavelet)
  L <- dict$length
  storage.mode(L) <- "integer"
  h <- dict$hpf
  storage.mode(h) <- "double"
  g <- dict$lpf
  storage.mode(g) <- "double"

  LL <- paste("LL", J, sep="")
  y.in <- y[[LL]]

  for(j in J:1) {
    LH <- paste("LH", j, sep="")
    HL <- paste("HL", j, sep="")
    HH <- paste("HH", j, sep="")
    
    m <- dim(y.in)[1]
    storage.mode(m) <- "integer"
    n <- dim(y.in)[2]
    storage.mode(n) <- "integer"
    x <- matrix(0, 2*m, 2*n)
    storage.mode(x) <- "double"

    out <- .C("two_D_idwt", as.double(y.in), as.double(y[[LH]]),
              as.double(y[[HL]]), as.double(y[[HH]]), m, n, L, h, g,
              "Y"=x, PACKAGE="EFDR")
    y.in <- out$Y
  }
  zapsmall(y.in)
}


