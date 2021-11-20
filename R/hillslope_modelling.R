
#'  Hillslope steady state profiles
#'
#' Compute hillslope topography steady state profiles, in the case of linear and nonlinear sediment flux
#'
#' If the threshold hillslope gradient `Sc` is specified the profile is computed in the non-linear case, otherwise linear.
#'
#' @param x Horizontal position along the hillslope (x=0 at hilltop, m)
#' @param E Denudation rate (m/a)
#' @param K Diffusion coefficient (m2/a)
#' @param L hillslope length (m)
#' @param beta rock to teglith density ratio (dimensionless)
#' @param Sc Critical hillslope gradient (optional, m/m)
#'
#' @return Hillslope elevation (relative to base of hillslope) at the position specified in `x`.
#'
#' @export
#'
#' @examples
#' x = seq(0,100,by=1)
#' L =  max(x) # hillslope length
#' E = 100/1e6 # erosion rate (m/Ma -> m/a)
#' K = 0.01 # diffusion coefficient (m2/a)
#' beta = 2
#' Sc = 0.8 # (m/m)
#' znl = hillslope_profile(x,E,K,L,beta,Sc)
#' plot(x,znl,type="l",lwd=3,
#' xlab="Distance from hilltop (m)",ylab="Elevation (m)")
#' zl = hillslope_profile(x,E,K,L,beta)
#' lines(x,zl,col="blue")
hillslope_profile <- function(x,E,K,L,beta,Sc=NULL){
  if (is.null(Sc)){
    z = -1*E*beta/K/2*(x^2-L^2)
    z = z-min(z)
  } else {
    a = sqrt(1 + (2*beta*E*x/(K*Sc))^2)
    z = ((K*Sc^2)/(2*beta*E) * ( log(0.5*(a + 1)) - a + 1) )
    z = z-min(z)
  }
  return(z)
}


#' Ratio of the nonlinear to linear components of the sediment flux
#'
#' Compute ratio of the nonlinear to linear components of sediment flux at the base of the hillslope according to \insertCite{roering2001hillslope;textual}{geomorphology}
#'
#' @param U relative baselevel fall rate (m/a)
#' @param K hillslope diffusion coefficient (m2/a)
#' @param beta rock-to-regolith density ratio (-)
#' @param Sc critical hillslope gradient (m/m)
#' @param L hillslope length (m)
#'
#' @return Value of Psi, which can be a vector
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
hillslope_psi <- function(U,K,L,beta,Sc) {
  f1=(2*beta*U*L)/Sc
  f2=((-K+sqrt(K^2+f1^2))/f1)^2
  return(f2/(1-f2))
}




#' Hillslope adjustement time
#'
#'
#' Compute exponential equilibrium adjustment timescale for sediment flux or hillslope morphology according to \insertCite{roering2001hillslope;textual}{geomorphology}
#'
#' @param U relative baselevel fall rate (m/a)
#' @param K hillslope diffusion coefficient (m2/a)
#' @param beta rock-to-regolith density ratio (-)
#' @param Sc critical hillslope gradient (m/m), optional if not specified the linear adjustement timescale is returned
#' @param L hillslope length (m)
#' @param A constant from \insertCite{roering2001hillslope;textual}{geomorphology}, optional (default 0.405)
#' @param B constant from \insertCite{roering2001hillslope;textual}{geomorphology}, optional (default 1.74)
#'
#' @return
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
hillslope_tau <- function(U,K,L,beta,Sc=NULL,A=0.4052847,B=1.74) {
  tau_l = A*L^2/K
  if (is.null(Sc)){
    return(tau_l)
  } else {
    psi = hillslope_psi(U,K,L,beta,Sc)
    return(tau_l/(1+psi)^B)
  }
}



#' Compute hillslope non-dimensional Erosion rate E*
#'
#' Compute hillslope non-dimensional Erosion Rate E*  according to \insertCite{roering2007functional;textual}{geomorphology}
#'
#' @param U relative baselevel fall rate (m/a)
#' @param K hillslope diffusion coefficient (m2/a)
#' @param beta rock-to-regolith density ratio (-)
#' @param Sc critical hillslope gradient (m/m)
#' @param L hillslope length (m)
#'
#' @return Value of E*, which can be a vector
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
hillslope_Es <- function(U,K,L,beta,Sc) {
  return((U*2*L*beta)/(K*Sc))
  }









#' Nonlinear hillslope evolution solver
#'
#' Solve nonlinear hillslope evolution using \insertCite{perron2011numerical;textual}{geomorphology} implicit formulation over one time step
#'

#'
#'
#' @param z elevation vector (m)
#' @param K hillslope diffusion coefficient (m2/a)
#' @param beta rock-to-regolith density ratio (-)
#' @param Sc critical hillslope gradient (m/m)
#' @param U relative baselevel fall rate (m/a)
#' @param dx horizontal distance increment (m)
#' @param dt time step (a)
#'
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @return Elevation vector at the end of time step
perron_nlhs<-function(z,K,beta,Sc,U,dx,dt){
nx = length(z)
F=getF(z,K,beta,Sc,nx,dx)
f=getf(z,K,beta,Sc,nx,dx,U)
B=getB(z,f,F,nx,dt)
A=getA(z,f,F,nx,dt)
z=solve(A,B)
return(z)
}


#' Perron 2011
#' @param z elevation vector (m)
#' @param K hillslope diffusion coefficient (m2/a)
#' @param prps rock-to-regolith density ratio (-)
#' @param Sc critical hillslope gradient (m/m)
#' @param nx number of space increment
#' @param dx horizontal distance increment (m)
#'
#' @return
getF<-function(z,K,prps,Sc,nx,dx){
  a=K/prps
  b=1/Sc^2
  F=matrix(0,nx,3)
  # create arrays of the functions f and F in equation 23
  # column 1 is left-hand neighbor (j-1), column 2 is right-hand neighbor
  # (j+1), column 3 is j
  for(i in 2:(nx-1)) {
    zp1=z[i+1]
    zm1=z[i-1]
    z1 = 0.5*(zp1-zm1)/dx
    z2 = (zp1-2*z[i]+zm1)/dx^2
    c=1/(1-b*z1^2)
    d1=1+b*c*z1^2
    d2=1+2*b*c*z1^2
    #
    F[i,3]=-2*a*c*d2/dx^2 # i
    F[i,2]=a*c*d2/dx^2 + a*b*c^2*d2*z1*z2/dx + a*b*c^2*d1*z1*z2/dx # i+1 (right)
    F[i,1]=a*c*d2/dx^2 - a*b*c^2*d2*z1*z2/dx - a*b*c^2*d1*z1*z2/dx # i+1 (left)
  }
  return(F)
}

#' Perron 2011
#' @param z elevation vector (m)
#' @param K hillslope diffusion coefficient (m2/a)
#' @param prps rock-to-regolith density ratio (-)
#' @param Sc critical hillslope gradient (m/m)
#' @param nx number of space increment
#' @param dx horizontal distance increment (m)
#' @param U relative baselevel fall rate (m/a)
#'
#' @return
getf<-function(z,K,prps,Sc,nx,dx,U){
  a=K/prps
  b=1/Sc^2
  f=matrix(0,nx,1)
  # create arrays of the functions f and F in equation 23
  # column 1 is left-hand neighbor (j-1), column 2 is right-hand neighbor
  # (j+1), column 3 is j
  for(i in 2:(nx-1)) {
    zp1=z[i+1]
    zm1=z[i-1]
    z1 = 0.5*(zp1-zm1)/dx
    z2 = (zp1-2*z[i]+zm1)/dx^2
    f[i] = U + a*z2/(1-b*z1^2)*(1+2*b*z1^2/(1-b*z1^2))
  }
  return(f)
}

#' Perron 2011
#' @param z elevation vector (m)
#' @param f f
#' @param F F
#' @param nx number of space increment
#' @param dt time step (a)
#'
#' @return
getB<-function(z,f,F,nx,dt){
  B=matrix(0,nx,1)
  di = c(-1,1) # vector, and matrix column, offset for each neighbor [L R]
  for(i in 1:nx) { # for each elevation
    # first 2 terms of B
    B[i] = (1 - dt*F[i,3])*z[i] + dt*f[i]
    for(j in 1:2) { # loop through neighbors
      n = i + di[j] # vector and matrix column index of neighbor element
      if (n>0 & n<(nx+1)) {  # don't go off boundaries
        B[i] = B[i] - dt*F[i,j]*z[n] # neighbor terms of B
      }
    }
  }
  return(B)
}


#' Perron 2011
#' @param z elevation vector (m)
#' @param f f
#' @param F F
#' @param nx number of space increment
#' @param dt time step (a)
#'
#' @return
getA<-function(z,f,F,nx,dt){
  A=matrix(0,nx,nx)
  di = c(-1,1) # vector, and matrix column, offset for each neighbor [L R]
  for(i in 1:nx) { # for each elevation
    # diagonal elements of A
    A[i,i] = 1 - dt*F[i,3]
    for(j in 1:2) { # loop through neighbors
      n = i + di[j] # vector and matrix column index of neighbor element
      if (n>0 & n<(nx+1)) {  # don't go off boundaries
        A[i,n] = -dt*F[i,j] # off-diagonal elements of A
      }
    }
  }
  return(A)
}






