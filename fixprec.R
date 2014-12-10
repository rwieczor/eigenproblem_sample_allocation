#
# This file contains R functions implementing the algorithms described  
# in the paper: Wesolowski J., Wieczorkowski R.:
# "An eigenproblem approach to optimal equal-precision sample allocation in subpopulations"
# (submitted to Communications in Statistics - Theory and Methods)
#
# actualization date: 10.12.2014
#


fixprec_HR_SRSWOR<-function(m,n,data,J="sub",H="h",N_jhi="N_jhi",M_jh="M_jh",
                     z_jhi="z_jhi",t_jhi="t_jhi",S2_jhi="S2_jhi",
                     omega_jhi="omega_jhi",T_j="T_j", gamma_jh="gamma_jh")
#
# Function for equal-precision optimal allocation in subpopulations in two-stage 
# sampling with Hartley-Rao pi-ps scheme at the first stage and the SRSWOR at the second,
# see Theorem 3.5 from the article
#
# Arguments:
#   m - sample size in terms of PSU
#   n - sample size in terms of SSU
#   data - data frame with input information for each stratum, subpopulation and PSU
#          (should contains identification variable for PSU)
#   J - name of the variable denoting subpopulations
#   H - name of the variable denoting strata
#   N_jhi - name of the variable with population sizes in terms of SSU
#   M_jh - name of the variable with population sizes in terms of PSU
#   z_jhi - name of the auxiliary variable \tilde{z} 
#   t_jhi - name of the variable with totals of surveyed variable
#   S2_jhi - name of the variable with population variances of surveyed variable
#   omega_jhi - name of the variable omega
#   T_j - name of the variable with totals of surveyed variable
#   gamma_jh - name of the variable gamma
#
# Values:
#   output data frame contains input data appended with new columns 'm_jhi' and 'n_jhi',
#   which give optimal allocation of numbers of PSU and SSU in subpulations, 
#   strata and PSU.
#
{ 
  stopifnot(all(data[[gamma_jh]]>0))
  
  d0<-unique(data[c(J,H,gamma_jh)])
  A_jh<-data.frame(d0[,c(J,H)],A_jh=d0[[gamma_jh]])
  
  alfa_jhi<-data.frame(data[,c(J,H)],alfa_jhi=data[[z_jhi]])
  B_jhi<-data.frame(data[,c(J,H)],B_jhi=(d[[N_jhi]]^2)*d[[S2_jhi]]/((d[[T_j]]^2)*d[[z_jhi]]))
  c_j<-aggregate(d[[omega_jhi]]*d[[z_jhi]]/(d[[T_j]]^2),list(d[[J]]),sum)
  names(c_j)=c(J,"c_j")
  
  M_jh<-unique(d[,c(J,H,M_jh)])
  N_jhi<-d[,c(J,H,PSU,N_jhi)]
  
  
  nJ<-length(table(c_j[[J]]))
  a_j<-matrix(tapply(sqrt(A_jh$A_jh/m),A_jh[[J]],sum),nJ,1)
  b_j<-matrix(tapply(sqrt(alfa_jhi$alfa_jhi*B_jhi$B_jhi/n),alfa_jhi[[J]],sum),nJ,1)
  
  D.matrix<-a_j%*%t(a_j) + b_j%*%t(b_j) - diag(c_j$c_j)
  
  m1<-eigen(D.matrix)
  # must be unique positive eigenvalue
  lambda<-m1$values
  lambda<-lambda[lambda>0]
  #print(lambda)
  #print(length(lambda))
  if (length(lambda)>1) stop("Positive eigenvalue is not unique - solution do not exists !")
  
  cv_optimal<-m1$values[1]  # maximum eigenvalue
  cat("CV optimal (in %) = ",100*sqrt(cv_optimal),"\n") 
  v<-(-m1$vectors[,1])  # corresponding eigenvector
  #sum(v)
  
  v_j<-data.frame(J=c_j[[J]],v=v)
  names(v_j)[1]<-J
  #v_j[[J]]<-as.factor(v_j[[J]])
  
  A_jh<-merge(A_jh,v_j,by=J)
  
  tmp<-with(A_jh,v*sqrt(A_jh))
  m_jh<-( m*((A_jh$v*sqrt(A_jh$A_jh))/sum(tmp)))
  #m_jh<-round(m_jh)
  sum(m_jh)
  
  B_jhi<-cbind(B_jhi,alfa_jhi=alfa_jhi$alfa_jhi)
  A_jh<-cbind(A_jh,m_jh=m_jh)
  B_jhi<-merge(B_jhi,A_jh)
  
  tmp<-with(B_jhi,v*sqrt(alfa_jhi*B_jhi))
  n_jhi<-((B_jhi$v*sqrt(B_jhi$B_jhi/B_jhi$alfa_jhi))/sum(tmp))*n/B_jhi$m_jh 
  n_jhi<-round(n_jhi)
  n_jhi<-pmax(n_jhi,1)
  sum(n_jhi)
  
  A_jh$m_jh<-round(A_jh$m_jh)
  out<-merge(N_jhi,M_jh)
  out<-cbind(out,n_jhi=n_jhi)
  out<-merge(out,A_jh[,c(J,H,"m_jh")])
  out<-out[order(out[[J]],out[[H]]),]
  
  return(out)
}






fixprec_SRSWOR<-function(n,data,J="sub",N_jh="N_jh",S2_jh="S2_jh",t_j="T_j")
#
# Function for equal-precision optimal allocation in single-stage sampling
# with subpopulations and strata, see Theorem 2.3 in the article
#
# Arguments:
#   n - sample size
#   data - data frame with input information for each stratum and subpopulation
#   J - name of the variable denoting subpopulations
#   N_jh - name of the variable with population sizes 
#   S2_jh - name of the variable with population variances of surveyed variable
#   t_j - name of the variable with totals of stratification variable in subpopulations
#
# Values:
#   output data frame contains input data appended with new column 'n_jh',
#   which gives optimal allocation in subpulations and strata.
#
{
  
  A_jh<-data[[N_jh]]*sqrt(data[[S2_jh]])/data[[t_j]]  
  B_jh<-(A_jh^2)/data[[N_jh]]
  
  A_j<-as.vector(tapply(A_jh,data[[J]],sum))
  a<-matrix(A_j,length(A_j),1)
  c<-as.vector(tapply(B_jh,data[[J]],sum))
  

  
  D.matrix<-(1/n)*(a%*%t(a))- diag(c)
  
  
  m1<-eigen(D.matrix)
  # must be unique positive eigenvalue
  lambda<-m1$values
  lambda<-lambda[lambda>0]
  #print(lambda)
  #print(length(lambda))
  if (length(lambda)>1) stop("Positive eigenvalue is not unique - solution do not exists !")
  
  cv_opt<-m1$values[1]  # maximum eigenvalue
  cat("cV optimal (in %)  = ",100*sqrt(cv_opt),"\n") 
  v<-(-m1$vectors[,1])  # corresponding eigenvalue
  
  v<-v*(n/as.numeric((t(a)%*%v)))
  out<-as.data.frame(table(data[[J]]))
  names(out)<-c(J,"freq")
  out<-cbind(out,n_jh=v)
  out<-merge(data,out,by=J)
  out$n_jh=round(A_jh*out$n_jh) # optimal allocation
  out$freq=NULL
  
  return(out)
}



