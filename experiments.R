# Code to do penalized spline with 
#   q-order difference penalties (q=0-3)
#    (if q is anything else, no penalty is applied)
#  Uses B-splines 
#
# scatterplot is (x[i],y[i]); k=number of knots (k-2 interior knots), 
# and pen=penalty
# returns GCV value for unconstrained & constrained splines
#   type:  1=monotone increasing 
#          2=monotone decreasing
#          3=convex
#          4=concave
#          5=convex increasing
#          6=convex decreasing
#          7=concave increasing
#          8=concave decreasing
#
#  returns: cfit = constrained fit
#           ucfit = unconstrained fit
#           cgcv = constrained GCV
#           ucgcv = unconstrained GCV
#           edfc = effective degrees of freedom for constrained fit
#           edfu = effective degrees of freedom for unconstrained fit
#           knots
#           xpl = grid of points for plotting smooth fits
#           cpl = constrained fit values at xpl
#           ucpl = unconstrained fit values at xpl
########################################################################
### need library(coneproj)
########################################################################
penspl=function(type,x,y,k,q,pen){
  n=length(y)
  sm=1e-8
  if(type<=2){
    ans=bqspline(x,k)
  }else{ans=bcspline(x,k)}
  m=length(ans$edges)/n
  delta=t(ans$edges)
  dmat=matrix(0,nrow=m-q,ncol=m)
  #  third-order
  if(q==3){
    for(i in 4:m){
      dmat[i-3,i-3]=1;dmat[i-3,i-2]=-3;dmat[i-3,i-1]=3;dmat[i-3,i]=-1
    }
  }
  # second order
  if(q==2){
    for(i in 3:m){
      dmat[i-2,i-2]=1;dmat[i-2,i-1]=-2;dmat[i-2,i]=1
    }
  }
  # first order
  if(q==1){
    for(i in 2:m){
      dmat[i-1,i-1]=1;dmat[i-1,i]=-1
    }
  }
  # zero order
  if(q==0){
    for(i in 1:m){
      dmat[i,i]=1
    }
  }
  # if q is anything else, no penalty
  qmat=delta%*%t(delta)+pen*t(dmat)%*%dmat
  umat=chol(qmat)
  if(type<=2){
    smat=ans$slopes
  }else{smat0=ans$d2}
  if(type==2){smat=-smat}
  if(type==3){smat=smat0}
  if(type==4){smat=-smat0}
  if(type==5){
    smat=matrix(0,ncol=m-1,nrow=m)
    smat[,1:(m-2)]=smat0
    smat[1,m-1]=-1;smat[2,m-1]=1
  }
  if(type==6){
    smat=matrix(0,ncol=m-1,nrow=m)
    smat[,1:(m-2)]=smat0
    smat[m-1,m-1]=1;smat[m,m-1]=-1
  }
  if(type==7){
    smat=matrix(0,ncol=m-1,nrow=m)
    smat[,1:(m-2)]=-smat0
    smat[m-1,m-1]=-1;smat[m,m-1]=1
  }
  if(type==8){
    smat=matrix(0,ncol=m-1,nrow=m)
    smat[,1:(m-2)]=-smat0
    smat[1,m-1]=1;smat[2,m-1]=-1
  }
  xpl=ans$xpl
  bpl=ans$bpl
  knots=ans$knots
  uinv=solve(umat)
  # make cone edges
  bmata=t(smat)%*%uinv
  bmat=matrix(0,ncol=m,nrow=m)
  if(type==3|type==4){
    perpmat=matrix(nrow=2,ncol=m)
    uvec=matrix(runif(2*m),nrow=m,ncol=2)
    perpmat=uvec-t(bmata)%*%solve(bmata%*%t(bmata))%*%bmata%*%uvec
    bmat[1:2,]=t(perpmat)
    bmat[3:m,]=bmata
    
  }else{
    uvec=runif(m)
    perpvec=uvec-t(bmata)%*%solve(bmata%*%t(bmata))%*%bmata%*%uvec
    bmat[1,]=perpvec
    bmat[2:m,]=bmata
  }
  edges=t(solve(bmat))
  ysend=t(uinv)%*%delta%*%y
  np=1;if(type==3|type==4){np=2}
  coef=coneproj::coneB(ysend,t(edges[(np+1):m,]),matrix(t(edges[1:np,]),ncol=np))
  beta=uinv%*%t(edges)%*%coef$coefs
  yhat=t(delta)%*%beta
  qinv=solve(qmat)
  pmat=t(delta)%*%qinv%*%delta
  ytil=pmat%*%y
  edfu=sum(diag(pmat))
  cv=sum((y-ytil)^2)/(1-edfu/n)^2
  # Compute cv for constrained fit
  index=coef$coefs>sm
  index[1]=TRUE
  gmat=edges[index,]
  if(length(gmat)/m==1){gmat=matrix(gmat,ncol=m)}
  pcmat=t(delta)%*%uinv%*%t(gmat)%*%solve(gmat%*%t(gmat))%*%gmat%*%t(uinv)%*%delta
  edfc=sum(diag(pcmat))
  cvc=sum((y-yhat)^2)/(1-edfc/n)^2
  ans=new.env()
  ans$cfit=yhat
  ans$ucfit=ytil
  ans$ucgcv=cv
  ans$cgcv=cvc
  ans$edfu=edfu
  ans$edfc=edfc
  ans$xpl=xpl
  ans$cpl=bpl%*%beta
  ans$ucpl=bpl%*%qinv%*%delta%*%y
  ans$knots=knots
  ans
}





############################################
#
# Here is the generic hinge code: 
# project y onto the cone defined by sigma 
# First p rows of sigma are unconstrained
#
# Next rows of sigma are edge vectors
# 
# returns coefficient vector
#
############################################
hinge=function(y,sigma,p){
  
  n=length(y)
  m=length(sigma)/n
  sm=0.00000001
  
  # Do Gram-Schmidt orthogonalization of first p rows
  
  if(p>0){
    q1=sigma[1,]/sqrt(sum(sigma[1,]^2))
  }
  q=matrix(q1,ncol=1,nrow=n)
  if(p>1){
    for(i in 2:p){
      q1=sigma[i,]-q%*%solve(t(q)%*%q)%*%t(q)%*%sigma[i,]
      q1=q1/sqrt(sum(q1^2))
      q=cbind(q,q1)
    }
  }
  l=p
  if(p>0){
    h=1:p
    ones=1:p*0+1
    r=diag(ones,nrow=p,ncol=p)
  }
  
  # add one edge to start
  
  check=0
  if(p==0){rhat=y}
  if(p>0){rhat=y-q%*%solve(t(q)%*%q)%*%t(q)%*%y}
  b2=sigma%*%rhat
  
  if(max(b2[(p+1):m])>sm){
    obs=(p+1):m
    i=min(obs[b2[(p+1):m]==max(b2[(p+1):m])])
    l=p+1
    if(p==0){
      q=sigma[i,]
      q=q/sqrt(sum(q^2))
      h=i
      r=matrix(1,nrow=1,ncol=1)
    }
    if(p>0){
      q1=sigma[i,]-q%*%solve(t(q)%*%q)%*%t(q)%*%sigma[i,]
      q1=q1/sqrt(sum(q1^2))
      q=cbind(q,q1)
      h[l]=i
      r=t(q)%*%t(sigma[h,])
    }
  }
  if(max(b2[(p+1):m])<(sm)){
    check=1
  }
  if(check==1){a=t(q)%*%y}
  # LOOP starts here:
  nrep=0
  while(check==0 & nrep<1000){
    nrep=nrep+1
    # Fit data to current EDGES
    
    a=t(q)%*%y
    
    # check if convex: 
    #first find the b vector:
    b=1:l*0
    if(l>1){b[l]=a[l]/r[l,l]
    for( j in (l-1):1){
      b[j]=a[j]
      for(i in (j+1):l){
        b[j]=b[j]-r[j,i]*b[i]
      }
      b[j]=b[j]/r[j,j]
    }
    }else{b[l]=a[l]/r}
    
    #check to see if b positive
    
    if(l>p){
      obs=(p+1):l
      i=obs[b[(p+1):l]==min(b[(p+1):l])]
      if(b[i]<(-sm)){
        # if not, remove hinge, make new q and r
        c1=0
        if(i>1){h=c(h[1:(i-1)],h[(i+1):l])}else{h=h[2:l]}
        l=l-1
        if(i>1){
          q=q[,1:(i-1)]
          for( j in i:l){
            qnew=sigma[h[j],]-q%*%t(q)%*%sigma[h[j],]
            qnew=qnew/sqrt(sum(qnew^2))
            q=cbind(q,qnew)
          }
          r=t(q[,1:l])%*%t(sigma[h,])
        }
        if(i==1&l>1){
          q[,1]=sigma[h[1],]/sqrt(sum(sigma[h[1],]^2))
          for(j in 2:l){
            qnew=sigma[h[j],]-q[,1:(j-1)]%*%t(q[,1:(j-1)])%*%sigma[h[j],]
            qnew=qnew/sqrt(sum(qnew^2))
            q=cbind(q[,1:(j-1)],qnew)
          }
          r=t(q[,1:l])%*%t(sigma[h,])
        }	
        if(i==1&l==1){
          q[,1]=sigma[h[1],]/sqrt(sum(sigma[h[1],]^2))
          r=matrix(1,nrow=1,ncol=1)
        }	
        q=q[,1:l]
      }
    }
    if(b[i]>(-sm)) {
      c1=1
      #
      # now see if we need to add another hinge
      #
      theta=q%*%t(q)%*%y
      rhat=y-theta
      b2=sigma%*%rhat
      
      # check to see if b2 negative
      
      obs=(p+1):m
      i=min(obs[b2[(p+1):m]==max(b2[(p+1):m])])
      if(l<m&l>0){if(b2[i]>sm){
        l=l+1
        qnew=sigma[i,]-q%*%t(q)%*%sigma[i,]
        qnew=qnew/sqrt(sum(qnew^2))
        q=cbind(q,qnew)
        h=c(h,i)
        r=t(q)%*%t(sigma[h,])
        c2=0
      }}
      if(b2[i]<sm){c2=1}
      check=c1*c2
      h
    }
  }
  # find coefficient vector
  b=1:l*0
  b[l]=a[l]/r[l,l]
  if(l>1){
    for( j in (l-1):1){
      b[j]=a[j]
      for(i in (j+1):l){
        b[j]=b[j]-r[j,i]*b[i]
      }
      b[j]=b[j]/r[j,j]
    }
  }
  coef=1:m*0
  coef[h]=b
  coef
}
########################################
#       MAKE THE EDGE VECTORS          #
########################################

##############################################################
# B-spline quadratic basis
# returns basis functions and slopes of basis functions
# at the knots
bqspline=function(x,m){
  tk=0:(m-1)/(m-1)*(max(x)-min(x))+min(x)
  k=3
  t=1:(m+2*k-2)*0
  t[1:(k-1)]=min(x);t[(m+k):(m+2*k-2)]=max(x)
  t[k:(m+k-1)]=tk
  n=length(x)
  sm=1e-8
  h=t[4]-t[3]
  
  bmat=matrix(1:(n*(m+k-2))*0,nrow=n)
  index=x>=t[3]&x<=t[4]
  bmat[index,1]=(t[4]-x[index])^2
  bmat[index,2]=2*(x[index]-t[2])*(t[4]-x[index])+(t[5]-x[index])*(x[index]-t[3])
  index=x>=t[4]&x<=t[5]
  bmat[index,2]=(t[5]-x[index])^2
  for( j in 3:(m-1) ){
    index=x>=t[j]&x<=t[j+1]
    bmat[index,j]=(x[index]-t[j])^2
    index=x>=t[j+1]&x<=t[j+2]
    bmat[index,j]=(x[index]-t[j])*(t[j+2]-x[index])+(x[index]-t[j+1])*(t[j+3]-x[index])
    index=x>=t[j+2]&x<=t[j+3]
    bmat[index,j]=(t[j+3]-x[index])^2
  }
  index=x>=t[m]&x<=t[m+1]
  bmat[index,m]=(x[index]-t[m])^2
  index=x>=t[m+1]&x<=t[m+2]
  bmat[index,m]=(x[index]-t[m])*(t[m+2]-x[index])+2*(x[index]-t[m+1])*(t[m+3]-x[index])
  index=x>=t[m+1]&x<=t[m+2]
  bmat[index,m+1]=(x[index]-t[m+1])^2
  
  #################################################
  # plotting splines
  
  xpl=0:1000/1000*(max(x)-min(x))+min(x)
  bpl=matrix(1:(1001*(m+k-2))*0,nrow=1001)
  index=xpl>=t[3]&xpl<=t[4]
  bpl[index,1]=(t[4]-xpl[index])^2
  bpl[index,2]=2*(xpl[index]-t[2])*(t[4]-xpl[index])+(t[5]-xpl[index])*(xpl[index]-t[3])
  index=xpl>=t[4]&xpl<=t[5]
  bpl[index,2]=(t[5]-xpl[index])^2
  for( j in 3:(m-1) ){
    index=xpl>=t[j]&xpl<=t[j+1]
    bpl[index,j]=(xpl[index]-t[j])^2
    index=xpl>=t[j+1]&xpl<=t[j+2]
    bpl[index,j]=(xpl[index]-t[j])*(t[j+2]-xpl[index])+(xpl[index]-t[j+1])*(t[j+3]-xpl[index])
    index=xpl>=t[j+2]&xpl<=t[j+3]
    bpl[index,j]=(t[j+3]-xpl[index])^2
  }
  index=xpl>=t[m]&xpl<=t[m+1]
  bpl[index,m]=(xpl[index]-t[m])^2
  index=xpl>=t[m+1]&xpl<=t[m+2]
  bpl[index,m]=(xpl[index]-t[m])*(t[m+2]-xpl[index])+2*(xpl[index]-t[m+1])*(t[m+3]-xpl[index])
  index=xpl>=t[m+1]&xpl<=t[m+2]
  bpl[index,m+1]=(xpl[index]-t[m+1])^2
  
  
  #################################################
  
  slopes=matrix(0,ncol=m,nrow=m+k-2)
  slopes[1,1]=-2*h
  slopes[m+k-2,m]=2*h
  slopes[2,1]=4*h
  slopes[2,2]=-2*h
  if(m==4){slopes[3,2]=2*h;slopes[3,3]=-2*h}
  slopes[m+k-3,m]=-4*h
  slopes[m+k-3,m-1]=2*h
  if(m>4){
    for(j in 3:(m+k-4)){
      slopes[j,j-1]=2*h
      slopes[j,j]=-2*h
    }
  }
  
  bmat[,1]=bmat[,1]*2
  bmat[,m+1]=bmat[,m+1]*2
  slopes[1,]=slopes[1,]*2
  slopes[m+1,]=slopes[m+1,]*2
  bpl[,1]=bpl[,1]*2
  bpl[,m+1]=bpl[,m+1]*2
  mb=max(bpl)
  slopes=slopes/mb
  bpl=bpl/mb
  bmat=bmat/mb
  
  
  ans=new.env()
  ans$edges=bmat
  ans$slopes=slopes
  ans$knots=tk
  ans$xpl=xpl
  ans$bpl=bpl
  ans
}

##############################################################
# cubic b-splines
# returns basis functions and 2nd derivative of basis functions
# at the knots
##############################################################
bcspline=function(x,m){
  tk=0:(m-1)/(m-1)*(max(x)-min(x))+min(x)
  k=4
  t=1:(m+2*k-2)*0
  t[1:(k-1)]=min(x);t[(m+k):(m+2*k-2)]=max(x)
  t[k:(m+k-1)]=tk
  n=length(x)
  sm=1e-8
  
  m0=matrix(1:(n*(m+k+1))*0,nrow=n)
  for( j in 4:(m+k-1) ){
    if(t[j]<t[j+1]-sm){
      index=x<=t[j+1]&x>t[j]
      m0[index,j]=1
    }
  }
  m0[1,k]=1
  
  m1=matrix(1:(n*(m+k))*0,nrow=n)
  for( j in 3:(m+k-1) ){
    index=x>t[j]&x<=t[j+2]
    if(t[j+1]>t[j]+sm){
      p1=(x[index]-t[j])/(t[j+1]-t[j])*m0[index,j]
    }else{p1=0}
    if(t[j+2]>t[j+1]+sm){
      p2= (t[j+2]-x[index])/(t[j+2]-t[j+1])*m0[index,j+1]
    }else{p2=0}
    m1[index,j]=p1+p2
  }
  imin=x==min(x)
  m1[imin,k-1]=1
  
  
  m2=matrix(1:(n*(m+k-1))*0,nrow=n)
  for( j in 2:(m+k-1) ){
    index=x>t[j]&x<=t[j+3]
    if(t[j+2]>t[j]+sm){
      p1=(x[index]-t[j])/(t[j+2]-t[j])*m1[index,j]
    }else{p1=0}
    if(t[j+3]>t[j+1]+sm){
      p2=(t[j+3]-x[index])/(t[j+3]-t[j+1])*m1[index,j+1]
    }else{p2=0}
    m2[index,j]=p1+p2
  }
  m2[imin,k-2]=1
  
  
  m3=matrix(1:(n*(m+k-2))*0,nrow=n)
  for( j in 1:(m+k-2) ){
    index=x>=t[j]&x<=t[j+4]
    if(t[j+3]>t[j]+sm){
      p1=(x[index]-t[j])/(t[j+3]-t[j])*m2[index,j]
    }else{p1=0}
    if(t[j+4]>t[j+1]+sm){
      p2=(t[j+4]-x[index])/(t[j+4]-t[j+1])*m2[index,j+1]
    }else{p2=0}
    m3[index,j]=p1+p2
  }
  
  #  plotting splines
  
  np=1000
  xpl=0:(np-1)/(np-1)*(max(x)-min(x))+min(x)
  m0pl=matrix(1:(np*(m+k+1))*0,nrow=np)
  for( j in 4:(m+k-1) ){
    if(t[j]<t[j+1]-sm){
      index=xpl<=t[j+1]&xpl>t[j]
      m0pl[index,j]=1
    }
  }
  m0pl[1,k]=1
  m1pl=matrix(1:(np*(m+k))*0,nrow=np)
  for( j in 3:(m+k-1) ){
    index=xpl>t[j]&xpl<=t[j+2]
    if(t[j+1]>t[j]+sm){
      p1=(xpl[index]-t[j])/(t[j+1]-t[j])*m0pl[index,j]
    }else{p1=0}
    if(t[j+2]>t[j+1]+sm){
      p2= (t[j+2]-xpl[index])/(t[j+2]-t[j+1])*m0pl[index,j+1]
    }else{p2=0}
    m1pl[index,j]=p1+p2
  }
  m1pl[1,k-1]=1
  m2pl=matrix(1:(np*(m+k-1))*0,nrow=np)
  for( j in 2:(m+k-1) ){
    index=xpl>t[j]&xpl<=t[j+3]
    if(t[j+2]>t[j]+sm){
      p1=(xpl[index]-t[j])/(t[j+2]-t[j])*m1pl[index,j]
    }else{p1=0}
    if(t[j+3]>t[j+1]+sm){
      p2=(t[j+3]-xpl[index])/(t[j+3]-t[j+1])*m1pl[index,j+1]
    }else{p2=0}
    m2pl[index,j]=p1+p2
  }
  m2pl[1,k-2]=1
  m3pl=matrix(1:(np*(m+k-2))*0,nrow=np)
  for( j in 1:(m+k-2) ){
    index=xpl>=t[j]&xpl<=t[j+4]
    if(t[j+3]>t[j]+sm){
      p1=(xpl[index]-t[j])/(t[j+3]-t[j])*m2pl[index,j]
    }else{p1=0}
    if(t[j+4]>t[j+1]+sm){
      p2=(t[j+4]-xpl[index])/(t[j+4]-t[j+1])*m2pl[index,j+1]
    }else{p2=0}
    m3pl[index,j]=p1+p2
  }
  
  # matrix of second derivatives
  
  secder=matrix(0,ncol=m,nrow=m+k-2)
  secder[1,1]=6
  secder[2,1]=-9
  secder[2,2]=3/2
  secder[3,1]=3
  secder[3,2]=-5/2
  secder[3,3]=1
  if(m>4){
    for(j in 4:(m-1)){
      secder[j,j-2]=1;secder[j,j-1]=-2;secder[j,j]=1
    }
  }
  secder[m,m-2]=1
  secder[m,m-1]=-5/2
  secder[m,m]=3
  secder[m+1,m-1]=3/2
  secder[m+1,m]=-9
  secder[m+2,m]=6
  
  ans=new.env()
  ans$bpl=m3pl
  ans$xpl=xpl
  ans$edges=m3
  ans$d2=secder
  ans$knots=tk
  ans
}


############################################################

# Section 0                                                #

###########################################################


# This skript requieres the user to use the code provided by 
# Mery C. MEyer on her Webpage: 
#   https://www.stat.colostate.edu/~meyer/penspl.htm



#############################################################

#               Section 1                                  #

############################################################

# Here we define a few functions we use to compute the errors 
# in the simulation

divide_ans <- function(ans,runs){
  
  #           cfit = constrained fit
  #           ucfit = unconstrained fit
  #           cgcv = constrained GCV
  #           ucgcv = unconstrained GCV
  #           edfc = effective degrees of freedom for constrained fit
  #           edfu = effective degrees of freedom for unconstrained fit
  #           knots
  #           xpl = grid of points for plotting smooth fits
  #           cpl = constrained fit values at xpl
  #           ucpl = unconstrained fit values at xpl
  
  ans$cfit = ans$cfit / runs
  ans$ucfit = ans$ucfit / runs
  ans$cgcv = ans$cgcv / runs
  ans$ucgv = ans$ucgv / runs
  ans$cpl = ans$cpl / runs
  ans$ucpl = ans$ucpl / runs
  return(ans)
}

add_ans <- function(ans_1,ans_2){
  #           cfit = constrained fit
  #           ucfit = unconstrained fit
  #           cgcv = constrained GCV
  #           ucgcv = unconstrained GCV
  #           edfc = effective degrees of freedom for constrained fit
  #           edfu = effective degrees of freedom for unconstrained fit
  #           knots
  #           xpl = grid of points for plotting smooth fits
  #           cpl = constrained fit values at xpl
  #           ucpl = unconstrained fit values at xpl
  ans_1$cfit = ans_1$cfit+ans_2$cfit
  ans_1$ucfit = ans_1$ucfit + ans_2$ucfit
  ans_1$cgcv = ans_1$cgcv + ans_1$cgcv
  ans_1$ucgv = ans_1$ucgv + ans_2$ucgv
  ans_1$cpl = ans_1$cpl + ans_2$cpl
  ans_1$ucpl = ans_1$ucpl + ans_2$ucpl
  return(ans_1)
   
}




# THe function that runs the simulations


run_simulation <-function(runs,n_candidates, objective_function, objective_function_name){
  p = 3 # degree of splines used
  confidence_level = 0.99
  unconstrained_errors <- character(runs) # here we safe the error for each run
  constrained_errors <- character(runs)
  n_tracker = 0 # here we keep track which n we are using a.t.m.
  
  
                                  
  # Create the data frames.
  sim_results <- data.frame(
    runs="",
    datapoints="",
    knots="",
    objective_function="",
    constrained_errors = "",
    unconstrained_errors = "",
    error_difference_variance  = "",
    confidence_interval = "",
    p_value = "",
    stringsAsFactors = FALSE
  )


  
  for (n in n_candidates){
    n_tracker = n_tracker + 1
    
    for (i in 1:runs) {

      # create data
      x = runif(n)
      y = objective_function(x) + rnorm(n)
      # compute the estimator
      ans = penspl(1,x,y,round(3 * (n^(1/(2*p + 3))),digits = 0),3,0.1)
      # safe the computed estimator
      if (i == 1){
        ans_mean = ans
      }
      else{
        ans_mean = add_ans(ans_mean,ans)
      }
    
      y_clean = objective_function(ans$xpl)
      constrained_errors[i] = mean(abs(objective_function(ans$xpl)-ans$cpl))
      unconstrained_errors[i] = mean(abs(objective_function(ans$xpl)-ans$ucpl))
    }
    # compute the average of the estimators
    ans_mean = divide_ans(ans_mean,runs)
    unconstrained_errors = as.numeric(unconstrained_errors)
    constrained_errors = as.numeric(constrained_errors)
    
    
    
    if (n_tracker == 1){
      ans_means = list(ans_mean)
    }
    else{
      ans_means = c(ans_means, ans_mean)
    }

 
    test_result = t.test((constrained_errors-unconstrained_errors),c(),"two.sided",0, FALSE, confidence_level)
    
    sim_results = rbind(sim_results, data.frame(
      runs = runs,
      datapoints = n,
      knots = round(3*(n^(1/(2*p + 3)))),
      objective_function = c(objective_function_name),
      constrained_errors = c(mean(constrained_errors)),
      unconstrained_errors = c(mean(unconstrained_errors)),
      error_difference_variance = c(var(constrained_errors - unconstrained_errors)),
      confidence_interval = test_result$conf.int,
      p_value = test_result$p.value,
                               
      stringsAsFactors = FALSE)
    )
    
  
    
  }

  return(list("stats" = sim_results, "data" = ans_means))
} # simulation



# For the results to be plottet correctly one has to call run_simulation and plot_results in succsession for each n
# otherwise the results for n_1 are overwritten with the results of n_2 before they are plottet
plot_results <- function(runs, n, objective_function_name,  objective_function, data){
  
  true_y = objective_function(data$xpl)
  
  plot(data$xpl,true_y,type= "l",col = "black", ylab = "y", xlab = "x")
  
  lines(data$xpl,data$cpl, col = "blue")
  
  lines(data$xpl,data$ucpl,col="red")
  
  #lines(data$xpl,true_y,col = "black")
  
  legend("topleft", legend=c("constrained","unconstrained","objective"), col = c("blue","red","black"),lty=1:2,cex=0.8,title = paste(objective_function_name,", n =",as.character(n)))
  
} # plot routine
  


  
  
# The final six functions we used in the simulation
  f_1 <- function(x){
    return(10*(x-0.5)^3+0.5)
  }
  f_1_name = "x -> 10 * (x - 0.5)^3 +0.5"
  f_2 <- function(x){
    return(0.1*log(50 * x))
  }
  f_2_name = "x -> 0.1 * ln(50 * x)"
  

  
  f_5 <- function(x){
    return(x^3 - 0.2 *x + 2)
  }
  f_5_name = "x -> x^3 - 0.2 * x + 2"
  f_6_name = "x -> -1.75 * x^2 + 2.65 * x"
  f_6 <- function(x){
    return(-1.75 * x^2 + 2.65 * x)
  }
  
  function_8 <- function(x){
    f = 0
    if (x<0.35){
      f = 25 * (x-0.35)^3
    }
    if (x>0.65){
      f = 25* (x-0.65)^3
    }
    return(f)
  }
  function_8 = Vectorize(function_8, vectorize.args = "x")
  function_8_name = "x -> f(x)"
  
  
  function_9 <- function(x){
    f = 0
    if (x< 0.5){
      f = -(1 / exp(1/(0.5-x)))
    }
    if (x>0.5){
      f =  (1 / exp(1/(x - 0.5)))
    }
    return(f)
  }
  function_9 = Vectorize(function_9, vectorize.args = "x")
  
  f_3 = Vectorize(function_9,vectorize.args = "x")
  f_3_name = "x -> f_3(x)"
  
  f_4 = Vectorize(function_8,vectorize.args = "x")
  f_4_name = "x -> f_4(x)"
  
  
  # A list containing the functions their names and space to save the data from the simulation
  experiment_list = list(
    list(f = function(x) f_1(x), name = f_1_name, data_1 = NULL, data_2 = NULL, data_3 = NULL, data_4 = NULL, data_5 = NULL, data_6= NULL, data_7 = NULL, data_8 = NULL),
    list(f = function(x) f_2(x), name = f_2_name, data_1 = NULL, data_2 = NULL, data_3 = NULL, data_4 = NULL, data_5 = NULL, data_6= NULL, data_7 = NULL, data_8 = NULL),
    list(f = function(x) f_3(x), name = f_3_name, data_1 = NULL, data_2 = NULL, data_3 = NULL, data_4 = NULL, data_5 = NULL, data_6= NULL, data_7 = NULL, data_8 = NULL),
    list(f = function(x) f_4(x), name = f_4_name, data_1 = NULL, data_2 = NULL, data_3 = NULL, data_4 = NULL, data_5 = NULL, data_6= NULL, data_7 = NULL, data_8 = NULL),
    list(f = function(x) f_5(x), name = f_5_name, data_1 = NULL, data_2 = NULL, data_3 = NULL, data_4 = NULL, data_5 = NULL, data_6= NULL, data_7 = NULL, data_8 = NULL),
    list(f = function(x) f_6(x), name = f_6_name, data_1 = NULL, data_2 = NULL, data_3 = NULL, data_4 = NULL, data_5 = NULL, data_6= NULL, data_7 = NULL, data_8 = NULL)
    )
  
  





###############################################################

#              Section 2                                      #
  
###############################################################
  
  
  
# THis seciton is dedicated to creating figure 5
  
  
  figure_5_seed = 3441
  set.seed(figure_5_seed)
 

  
  
  n = 50
  x= runif(n)
  p = 3
  q = 3
  pen = 1
  type = 1
  
  knots = round(3 * (n^(1/(2*p + 3))))
  

  
  g <- function(x){
    return(100 * (0.6*(x-0.5))^5-0.05* 0.6*x - exp(0.6*(x- 1)))
  }
  
  g <- function(x){
    a = 10
    b = 0.5
    c = 0.5
    return(a * (x - b)^3 + c)
  }
  g <- function(x){
    a = 10
    b = 0.25
    c = 0.75
    d = 9/16
    if(x < 0.5){
      return(a * (x - b)^3 + b)
    }
    else{
      return( a * (x - c)^3 + d)
    }
    
  }
  g = Vectorize(g,vectorize.args = "x")
  
  name = "x -> 60 * (x - 1/2)^0.5 - 0.5 * x - exp(x - 1)"
  name = "x -> g(x)"
  
  y=g(x)+ rnorm(n)
  
  ans=penspl(type,x,y,knots,q,pen)
  
  y_clean = g(ans$xpl)
  
  plot(x,y)
  
  cfit = ans$cfit[order(x, decreasing=FALSE)]
  ucfit = ans$ucfit[order(x, decreasing = FALSE)]
  x = sort(x,decreasing = FALSE)
  error_c = mean(abs(g(ans$xpl)- ans$cpl))
  error_uc = mean(abs(g(ans$xpl)- ans$ucpl))
  
  lines(ans$xpl,ans$cpl, col = "blue")
  
  lines(ans$xpl,ans$ucpl,col= "red")
  
  lines(ans$xpl,y_clean,col = "black")
  
legend("bottomleft", 
         legend=c("constrained","unconstrained","objective"), 
         col = c("blue","red","black"),
         lty=1:2,
         cex=0.8,
         title = paste(name,", n =",as.character(n), ", knots = ", as.character(knots),", error =", as.character(error_c - error_uc))
)
  
##################################################################
  
#             Section 3                                         #
  
###################################################################
 
  
  
# The following paragra?h runs the simulation

sim_results = run_simulation(runs, c(50), function_1, "function 1")$stats
runs = 500
n_candidates = c(25,50,100,250,500,1000,2000)
i = 1
for (bundle in experiment_list){
    j = 1
    print(bundle[[2]])
    for (n in n_candidates){
      print(n)
      sim_results_new = run_simulation(runs,c(n),bundle[[1]],bundle[[2]])
      sim_results = rbind(sim_results,sim_results_new$stats)
      experiment_list[[i]][[2 + j]] = sim_results_new$data
      j = j +1
    }
    
    i = i + 1
}

  

# In the following are all other plots but figure 10 from the termpaper created
  
for (i in c(1,2,3,4,5,6)){
  x = experiment_list[[1]]$data_1[[1]]$xpl
  y = experiment_list[[i]]$f(x)
  name = experiment_list[[i]]$name
  plot(x, y, xlim = c(-0.2,1.2), type = "l", col = "black", ylab = "y", xlab = "x") 
  lines(x = experiment_list[[i]]$data_1[[1]]$xpl, y = experiment_list[[i]]$data_1[[1]]$cpl, xlim = c(0,1), type = "l", col = "green")
  lines(x = experiment_list[[i]]$data_2[[1]]$xpl, y = experiment_list[[i]]$data_2[[1]]$cpl, xlim = c(0,1), type = "l", col = "yellow")
  lines(x = experiment_list[[i]]$data_3[[1]]$xpl, y = experiment_list[[i]]$data_3[[1]]$cpl, xlim = c(0,1), type = "l", col = "blue")
  lines(x = experiment_list[[i]]$data_4[[1]]$xpl, y = experiment_list[[i]]$data_4[[1]]$cpl, xlim = c(0,1), type = "l", col = "red")
  lines(x = experiment_list[[i]]$data_5[[1]]$xpl, y = experiment_list[[i]]$data_5[[1]]$cpl, xlim = c(0,1), type = "l", col = "orange")
  lines(x = experiment_list[[i]]$data_6[[1]]$xpl, y = experiment_list[[i]]$data_6[[1]]$cpl, xlim = c(0,1), type = "l", col = "brown")
  lines(x = experiment_list[[i]]$data_7[[1]]$xpl, y = experiment_list[[i]]$data_7[[1]]$cpl, xlim = c(0,1), type = "l", col = "pink")
  legend("topleft", legend=c("n = 25","n = 50","n = 100", "n = 250", "n = 500", "n = 1000", "n = 2000", "objective"), col = c("green","yellow","blue","red","orange","brown","pink","black"),lty=1:2,cex=0.8,title = name)

}
  
for (i in c(1,2,3,4,5,6)){
    x = experiment_list[[1]]$data_1[[1]]$xpl
    y = experiment_list[[i]]$f(x)
    name = experiment_list[[i]]$name
    plot(x, y, xlim = c(-0.2,1.2), type = "l", col = "black", ylab = "y", xlab = "x") 
    lines(x = experiment_list[[i]]$data_1[[1]]$xpl, y = experiment_list[[i]]$data_1[[1]]$ucpl, xlim = c(0,1), type = "l", col = "green")
    lines(x = experiment_list[[i]]$data_2[[1]]$xpl, y = experiment_list[[i]]$data_2[[1]]$ucpl, xlim = c(0,1), type = "l", col = "yellow")
    lines(x = experiment_list[[i]]$data_3[[1]]$xpl, y = experiment_list[[i]]$data_3[[1]]$ucpl, xlim = c(0,1), type = "l", col = "blue")
    lines(x = experiment_list[[i]]$data_4[[1]]$xpl, y = experiment_list[[i]]$data_4[[1]]$ucpl, xlim = c(0,1), type = "l", col = "red")
    lines(x = experiment_list[[i]]$data_5[[1]]$xpl, y = experiment_list[[i]]$data_5[[1]]$ucpl, xlim = c(0,1), type = "l", col = "orange")
    lines(x = experiment_list[[i]]$data_6[[1]]$xpl, y = experiment_list[[i]]$data_6[[1]]$ucpl, xlim = c(0,1), type = "l", col = "brown")
    lines(x = experiment_list[[i]]$data_7[[1]]$xpl, y = experiment_list[[i]]$data_7[[1]]$ucpl, xlim = c(0,1), type = "l", col = "pink")
    legend("topleft", legend=c("n = 25","n = 50","n = 100", "n = 250", "n = 500", "n = 1000", "n = 2000", "objective"), col = c("green","yellow","blue","red","orange","brown","pink","black"),lty=1:2,cex=0.8,title = name)
    
}
  
  
####################################################################
  
#             Section 4                                           #
  
###################################################################
  
  
  
# In this section figure 10 is created


if(TRUE == FALSE){
  c_int_1 = c(0.045,0.025,0.015,0.01,0.01,0.005)
  c_int_2 = c(0.065,0.045,0.035,0.025,0.02,0.01)
  c_int_3 = c(0.075,0.045,0.035,0.03,0.02,0.015)
  c_int_4 = c(0.045,0.025,0.02,0.015,0.01,0.01)
  c_int_5 = c(0.06,0.035,0.025,0.02,0.01,0.01)
  c_int_6 = c(0.044,0.024,0.019,0.014,0.009,0.009)
  
  n_points = c(25,100,250,500,1000,2000)
  
  plot(n_points,c_int_1,type = "l", col = "black", ylab = "- error", xlab = n)
  lines(n_points,c_int_2, type = "l", col = "green")
  lines(n_points,c_int_3, type = "l", col = "red")
  lines(n_points,c_int_4, type = "l", col = "blue")
  lines(n_points,c_int_5, type = "l", col = "orange")
  lines(n_points,c_int_6, type = "l", col = "yellow")
  legend("topright", legend=c("f_1","f_2", "f_3", "f_4", "f_5", "f_6"), col = c("black","green","red","blue","orange","yellow"),lty=1:2,cex=0.8)
  
}


print(simulation)
#test_list = readRDS("simulation_results.RData")

#print(test_list[[2]]$data_4[[1]]$cpl)
print(ans$cgcv)
print(ans$xpl)
#  returns: cfit = constrained fit
#           ucfit = unconstrained fit
#           cgcv = constrained GCV
#           ucgcv = unconstrained GCV
#           edfc = effective degrees of freedom for constrained fit
#           edfu = effective degrees of freedom for unconstrained fit
#           knots
#           xpl = grid of points for plotting smooth fits
#           cpl = constrained fit values at xpl
#           ucpl = unconstrained fit values at xpl