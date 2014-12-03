#'@export
getVolume<-function(xyz,id,alpha,cas) {

if (nrow(xyz)!=length(id)) {stop("The xyz and id do not have the same length")}
if (cas!=TRUE & cas!=FALSE) {stop("The cas parameter is invalid. Please, use TRUE or FALSE")}

  xrange1<-range(xyz[,1])
  yrange2<-range(xyz[,2])
  ntree<-as.factor(id)
  N<-nlevels(ntree)
  
  xyzid<-cbind(xyz,id)
  
  repmat = function(X,m,n){

    mx = dim(X)[1]
    nx = dim(X)[2]
    matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
  }
  
  allc<-numeric(3)
  len=numeric()
  suppressWarnings(
  for (i in 1:N) {
    len[i]=dim(subset(xyzid,id==i))
    allc<-rbind(allc,subset(xyzid, id==i))
    
  })
  allc<-allc[-1,1:3]
  allc<-na.omit(rbind(allc,allc,allc))
  

  real_range=diff(apply(allc,2,range))
  all_scaled1=allc-repmat(t(as.matrix(apply(allc,2,min))),nrow(allc),1)
  all_scaled=all_scaled1/repmat(t(as.matrix(apply(all_scaled1,2,max))),nrow(all_scaled1),1)

  alpha=0.25
  alpha <- rep(alpha,N) 
  cuts<-c(0,cumsum(len)) 
  v<-numeric()
  volume_final<-numeric()

  for (i in 1:N){
    obs<-all_scaled[(cuts[i]+1):cuts[i+1],] 
    b=ashape3d(as.matrix(obs),alpha=alpha[i],pert = TRUE, eps = 1e-09) 
    borig<-b 
    borig$x<-borig$x*repmat(t(as.matrix(apply(all_scaled1,2,max))),len[i],1)
    borig$x<-borig$x+repmat(t(as.matrix(apply(allc,2,min))),len[i],1)
    if (cas==TRUE) {
      bg3d("white")
      plot(borig,indexAlpha = "all",clear=FALSE,col=rep(i,N),transparency = 1)}
    v[i]=volume_ashape3d(b)
    volume_final[i]=prod(v[i],real_range)
  }
  volume3d.trees<-NULL
  
  for( i in 1:N) {
    volume3d.trees[i]<-(volume_final[i]) 
  }
  volume3d.trees<-cbind(id=as.numeric(levels(factor(id))),V=round(volume3d.trees, digits=2))
  colnames(volume3d.trees)<-c("id","V")
  
 return(volume3d.trees)
 }

