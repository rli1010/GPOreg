#transformed linear model
trans.reg=function (formula, data, trans = TRUE) 
{#function for fitting the transformed linear model with the GPO
  call <- match.call()
  formula <- as.formula(formula)
  mf <- model.frame(formula, data)
  Ys <- model.response(mf)
  Zmat <- model.matrix(formula, data)
  n <- nrow(data)
  K.x <- ncol(Zmat) - 1
  pihat <- rowMeans(apply(Ys, 2, rank))/n
  pihat[pihat == 1] <- 0.9999
  pihat[pihat == 0] <- 1e-04
  if (trans) {
    etaihat <- qtruncnorm(pihat, -3, 3)
    g.deriv <- 1/dtruncnorm(etaihat, -3, 3)
  }
  else {
    etaihat <- pihat
    g.deriv <- 1
  }
  fit <- lm(as.formula(paste0("etaihat ~", as.character(formula)[[3]])), 
            data = data)
  coefs <- coef(fit)
  var.naive <- vcov(fit)
  s1i <- Zmat * as.vector((etaihat - Zmat %*% coefs))
  Y.ij1 <- apply(Ys, 2, function(x) matrix(rep(x, n), nrow = n, 
                                           byrow = T) <= matrix(rep(x, n), nrow = n, byrow = F))
  Y.ij2 <- apply(Ys, 2, function(x) (matrix(rep(x, n), nrow = n, 
                                            byrow = T) < matrix(rep(x, n), nrow = n, byrow = F)) + 
                   1/n)
  Y.ij <- (Y.ij1 + Y.ij2)/2
  ksi.ij <- matrix(rowMeans(Y.ij), nrow = n) - pihat
  s2i.ij <- array(NA, c(n, n, K.x + 1))
  for (p in 1:(K.x + 1)) s2i.ij[, , p] <- Zmat[, p] * g.deriv * 
    ksi.ij
  s2i <- apply(s2i.ij, 3, colMeans)
  negA <- crossprod(Zmat)/n
  negA.inv <- solve(negA)
  Var.B <- crossprod(s1i + s2i)/n^2
  var.ABA <- negA.inv %*% Var.B %*% t(negA.inv)
  z.value <- coefs/sqrt(diag(var.ABA))
  p.value <- 2 * (1 - pnorm(abs(z.value)))
  result=data.frame(cbind(coefs,sqrt(diag(var.ABA)),z.value,p.value))
names(result)=c("Coef","SE","Z","p")
return(result)
}

bbeta.norm=function (bbeta) { 
  if (sum(bbeta^2) > 1) {b1.2 <- -100}
  if (sum(bbeta^2) <= 1) {b1.2 <- sqrt(1 - sum(bbeta^2))}
  return(c(b1.2, bbeta))}

##monotonic index model
gen.init=function(x.formula,da,phat,n.init,ptrg=1:4,bb.init=NULL){
  da1=da;  da1$phat=phat
  lm.fit=summary(lm(as.formula(paste("phat~",x.formula,sep="")),data=da1))
  bb0=lm.fit$coefficients[-1,1:2]
  
  grid <- matrix(c(t(outer(rep(ptrg, each = n.init), bb0[,2]))),dim(bb0)[1])
  grid.pt=grid*matrix(runif(length(grid),-1,1),dim(grid)[1],dim(grid)[2])
  
  if(is.null(bb.init)){bb.init1=bb0[,1]}
  if(!is.null(bb.init)){bb.init1=sqrt(sum(bb0[,1]^2))*bb.init}
  
  b.init=cbind(as.vector(bb.init1),as.vector(bb.init1)+grid.pt)
  gen.init=apply(b.init,2,function(x) x/sqrt(sum(x^2))) }

mono_obj=function(bbeta,Zmat,phat,wi,type,Cn=NULL){
  # objective, kernel version, can set Cn extremely small to get indicator version
  Zb=as.vector(Zmat%*%matrix(bbeta.norm(bbeta)))
  mono_obj=-kernelC2(Zb,phat,wi,Cn)
  return(mono_obj)}

optim.m=function(par0,Zmat,phat,wi,type=1,Cn=NULL,control,method="Nelder-Mead"){
  fit=optim(par=par0,fn=mono_obj,Zmat=Zmat,phat=phat,wi=wi,type=type,Cn=Cn,control=control,method=method)
  optim.m=c(fit$convergence,fit$value,fit$par)}

find_sol=function(init.mat,Zmat,da,phat,wi,C0,naive.sigma,n.init,ptrg,method="Nelder-Mead"){
  n=dim(Zmat)[1];  
  Cn=naive.sigma*C0*n^(-1/3)
  
  out0=t(apply(init.mat,2,optim.m,Zmat=Zmat,phat=phat,wi=wi,type=2,Cn=Cn,control=list(maxit=1000),method=method)) #kernel
  out=out0[out0[,1]==0,]
  bb.fit.k=out[which.min(out[,2]),-(1:2)]
  bb.fit.k=bbeta.norm(bb.fit.k)
  return(list(bb.fit=bb.fit.k))}

mono.reg.sub=function(Ys,da,Zmat,wi=NULL,C0,ptrg=1:4,n.init,bb.init=NULL,varlist){
  n=dim(Ys)[1];   
  if(is.null(wi)){wi=rep(1,n)};
  
  phat0=rowMeans(apply(Ys, 2, myrank, wi=wi))/n; 
  phat=truncnorm:::qtruncnorm( pmin(pmax(phat0,0),1) ,-3,3)
  
  init.mat=gen.init(paste(varlist,collapse="+"),da,phat,n.init=n.init,ptrg,bb.init)[-1,]
  naive.b=bbeta.norm(init.mat[,1])
  bb.fit.k=naive.b; out2=data.frame(naive.b); row.names(out2)=varlist;
  naive.sigma=sd(Zmat%*%naive.b)
  
  fit05=find_sol(init.mat,Zmat,da,phat,wi,C0,naive.sigma=naive.sigma,n.init=5,ptrg=ptrg)
  out1=data.frame(varlist=varlist,coef=fit05$bb.fit)
  return(list(out1=out1))  }

mono.reg=function(Ylist,varlist,da,n.init0=20,resamp,C0=0.5,seed=1234){
  set.seed(seed)
  
  Ys=as.matrix(da[,Ylist])
  Zmat=as.matrix(da[,varlist])
  nn=dim(Zmat)[1]
  mono.est=mono.reg.sub(Ys,da,Zmat,wi=rep(1,nn),C0,ptrg=1:4,n.init=n.init0,varlist=varlist)
  print(mono.est$out1)
  print("Resampling in progress")
  
  bb.rmat=c()
  for(r.i in 1:resamp){
    mono.r=mono.reg.sub(Ys,da,Zmat,wi=rexp(nn),C0,ptrg=1:4,n.init=4,bb.init=mono.est$out1$coef,varlist)
    bb.rmat=rbind(bb.rmat,mono.r$out1$coef)
    print(r.i)}
  
  mono.est$out1$SE=apply(bb.rmat,2,sd)
  mono.est$out1$p=round(pnorm(-abs(mono.est$out1$coef/mono.est$out1$SE))*2,4)
  set.seed(NULL)
  return(list(est=mono.est,stgst=mono.est$stgst,phat=mono.est$phat))}
