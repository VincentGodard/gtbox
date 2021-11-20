


#' Compute probability density of discharge according to an inverse gamma law
#'
#' Q and Qb must have the same units
#'
#' See equation 11 of Dibiase and Whipple (2011)
#'
#' @param k Variability parameter
#' @param Q Discharge
#' @param Qb Mean annual discharge
#'
#' @return Probability density values for discharge
#' @export
#'
#' @examples
Inv_Gamma_dist <- function(k,Q,Qb){
  pred_ig=k^(k+1)/gamma(k+1)*exp(-1*k/(Q/Qb))*(Q/Qb)^(-1*(2+k))/Qb #(/Qb !)
  return(pred_ig)
}



#' Compute long-term channel incision according to a stream power-based instantaneous incision law and an inverse gamma distribution for discharge
#'
#' This function is actually solving equation 10 from Dibiase and Whipple (2011)
#'
#' @param k Variability parameter of the inverse gamma distribution for discharge
#' @param ks Channel steepness index
#' @param Psi Incision threshold
#' @param K erodibility
#' @param n slope exponent of the stream power law
#' @param gamma discharge exponent
#' @return Long-term incision value (m/s)
#' @export
#'
#' @examples
I_long_term <- function(k,ks,Psi,K,n,gamma){
  # check
  l = c(length(k),length(ks),length(Psi),length(K),length(n),length(gamma))
  nr = max(l)
  l = l[l != 1]
  if ( (length(l)>1) & (var(l) != 0)) {stop("If more than one vector in the inputs all lengths must be the same")}
  # expand
  E = rep(NA,nr)
  if (length(k)==1){k = rep(k,nr)}
  if (length(ks)==1){ks = rep(ks,nr)}
  if (length(Psi)==1){Psi = rep(Psi,nr)}
  if (length(K)==1){K = rep(K,nr)}
  if (length(n)==1){n = rep(n,nr)}
  if (length(gamma)==1){gamma = rep(gamma,nr)}
  # run
  for (i in 1:nr){
    Qc=Qcritical(ks[i],Psi[i],K[i],n[i],gamma[i]) # normalized Qc
    E[i] = integrate(f_int_10,lower=Qc,upper=Inf,k[i],ks[i],Psi[i],K[i],n[i],gamma[i],subdivisions=1e6,abs.tol=1e-9)$value
  }

  # nout = 1e4
  # for (i in 1:length(ks)){
  #   Qc=Qcritical(ks[i],Psi[i],K[i],n[i],gamma[i])
  #   Qm=Qc*1e3
  #   Q1 = 10^seq(log10(Qc),log10(Qm),length.out = nout)
  #   Q2 =seq(Qc,Qm,length.out = nout)
  #   Qs= sort(unique(c(Q1,Q2)))
  #   I = I_instantaneous(Qs,ks[i],Psi[i],K[i],n[i],gamma[i])
  #   igd = Inv_Gamma_dist(k[i],Qs,1)
  #   #  igd[igd<1e-8] = 0 # !!!
  #   E[i] = pracma::trapz(Qs,I*igd)
  # }
  return(E)
}
#
f_int_10 <- function(Qs,k,ks,Psi,K,n,gamma){
  Ig = Inv_Gamma_dist(k,Qs,1)
  val=Ig*I_instantaneous(Qs,ks,Psi,K,n,gamma)
  return(val)
}



#' Compute instantaneous incision
#'
#' See equation 8 of Dibiase and Whipple (2011)
#'
#' @param Qs Normalized discharge (to mean annual discharge)
#' @param ks Channel steepness index
#' @param Psi Incision threshold
#' @param K erodibility
#' @param n slope exponent of the stream power law
#' @param gamma discharge exponent
#' @return Instantaneous incision value (m/s)
#' @export
#'
#' @examples
I_instantaneous <- function(Qs,ks,Psi,K,n,gamma){
  I = K * (Qs)^gamma * ks^n - Psi
  I[I<0] = 0
  return(I)
}

#' Compute critical discharge (normalized) for the onset of incision
#'
#' Determined by setting I=0 in equation 8 (instantaneous discharge) of Dibiase and Whipple (2011)
#'
#' @param ks Channel steepness index
#' @param Psi Incision threshold
#' @param K erodibility
#' @param n slope exponent of the stream power law
#' @param gamma discharge exponent
#' @return Normalized critical discharge
#' @export
#'
#' @examples
Qcritical <- function(ks,Psi,K,n,gamma){
  Qc = (Psi/K)^(1/gamma)*ks^(-n/gamma)
  return(Qc)
}


