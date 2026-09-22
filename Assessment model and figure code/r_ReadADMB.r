#################################################### read report file
read.admb <- function(ifile)
  {  
    ret=read.fit(ifile)
    
    fn=paste(ifile,'.rep', sep='')
    A=read.rep(fn)
    A$fit=ret
    
    pfn=paste(ifile,'.psv',sep='')
    if(file.exists(pfn))
      A$post.samp=read.psv(pfn)
    
    return(A)
  }

read.fit <- function(ifile)
  {
    # __Example:             
    #	file <-("~/admb/simple")
    #	A <- reptoRlist(file)
    #	Note there is no extension on the file name.
    
    ## The following is a contribution from:
    ## Anders Nielsen that reads the par & cor files.
    ret<-list() 
    parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
                             what='', n=16, quiet=TRUE)[c(6,11,16)]) 
    ret$nopar<-as.integer(parfile[1]) 
    ret$nlogl<-parfile[2] 
    ret$maxgrad<-parfile[3] 
    file<-paste(ifile,'.cor', sep='') 
    lin<-readLines(file) 
    ret$npar<-length(lin)-2 
    ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
    sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
    ret$names<-unlist(lapply(sublin,function(x)x[2])) 
    ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
    ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
    ret$cor<-matrix(NA, ret$npar, ret$npar) 
    corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
    ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
    ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
    ret$cov<-ret$cor*(ret$std%o%ret$std)
    return(ret)
  }

read.rep <- function(fn)
  {
    # The following reads a report file
    # Then the 'A' object contains a list structure
    # with all the elemements in the report file.
    # In the REPORT_SECTION of the AMDB template use 
    # the following format to output objects:
    #  	report<<"object \n"<<object<<endl;
    #
    # The part in quotations becomes the list name.
    # Created By Steven Martell
    options(warn=-1)  #Suppress the NA message in the coercion to double
       
    ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
    idx=sapply(as.double(ifile),is.na)
    vnam=ifile[idx] #list names
    nv=length(vnam) #number of objects
    A=list()
    ir=0
    for(i in 1:nv)
    {
      ir=match(vnam[i],ifile)
      if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
      dum=NA
      if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
      if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))
      
      if(is.numeric(dum))#Logical test to ensure dealing with numbers
      {
        A[[vnam[i]]]=dum
      }
    }
    options(warn=0)
    
    return(A)
  }

read.psv <- function(fn, nsamples=10000)
  {
    #This function reads the binary output from ADMB
    #-mcsave command line option.
    #fn = paste(ifile,'.psv',sep='')
    filen <- file(fn, "rb")
    nopar <- readBin(filen, what = integer(), n = 1)
    mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
    mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
    close(filen)
    return(mcmc)
  }
  
############################################## read std file
###Read from a .std file produced by admb
###Read from a .std file produced by admb or a file of identical structure.
###@param name  path of file to read from
###@example: par_read("sea16.std")
###@value A list with components est and std for estimates and standard deviations.
###@author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}

read.std <- function(name)
{
  file_name = paste(name,".std",sep="") 
  tmp = read.table(file_name,header=T)
  est = tmp[,'value']
  std = tmp[,'std.dev']
  
  name.temp=as.character(tmp[,'name'])
  
  for(i in 1:nrow(tmp)){
  	if(grepl('[',tmp[i,'name'],fixed=T)){
  		j=regexpr('[',tmp[i,'name'],fixed=T)[1]
  		name.temp[i]=substring(tmp[i,'name'],1,j-1)
  	}
  }
  
  names(est) = name.temp
  names(std) = name.temp
   
  L1 = list()
  L2 = list()
   
  for(i in unique(names(std)))
  {
    L1[[i]] = est[names(est)==i]
    L2[[i]] = std[names(std)==i]
  }
  
  list(est=L1,std=L2)
}

############################################################# read par file
#' Read from a .par file produced by admb
#'
#' Read from a .par file produced by admb or a file of identical structure.
#' This reads the values of the objective function and the named parameters.
#' If some of the parameters are matrices, then their number of columns
#' should be supplied in a named list. (See example)
#' 
#' @param name  path of file to read from
#' @param ncols named list with the number of columns for matrix parameters
#' @example: par_read("sea16.par")
#' @example: par_read("sea16.rep",ncols=list(N=3,ogives=2,age_dist=50)) # N,ogives,age_dist are matrices
#' @value A list with components named as in the .par file, along with three components at the end: $n_par, $loglik, and $gradient which are extracted from the first line of the .par file.
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}

# par_read <- function(name,ncols=list()) # matrices must be specified by name and ncol
read.par <- function(name)
{
  file_name = paste(name,".par",sep="")
  
  tmp = scan(file_name,what="",quiet=T)
  tmp2 = split(tmp,cumsum(tmp=="#"))
  x = tmp2  
  x = x[-1]
  
  for(i in 1:length(x))
  {
    y = x[[i]]
    n = nchar(y[2])
    x[[i]] = as.numeric(y[-(1:2)])
    
    j1=regexpr('[',y[2],fixed=T)[1] ##has to do this to get rid of the quotes in list element name;otherwise, will cause problem when call list with variable name
    j2=regexpr(':',y[2],fixed=T)[1]
    
    part.1=substring(y[2],1,j1-1)
    part.2=substring(y[2],j1+1,j2-2)
            
    names(x)[i] = paste(part.1,'_',part.2,sep='')
   }
  
  # # Convert to matrix for those arguments relevant
  # if(length(ncols)>0)
    # for(i in 1:length(ncols))
    # {
      # NN = names(ncols)[i]
      # x[[NN]] <- matrix(x[[NN]],ncol=ncols[[i]],byrow=T)
    # }
  
    x$n_par = -as.numeric(tmp2[[1]][6])
    x$loglik = -as.numeric(tmp2[[1]][11])
    x$gradient = -as.numeric(tmp2[[1]][16])  
  x
} 

######################################################### read other files
#' make a correlation matrix from a covariance matrix
#'
#' If m is a covariance matrix, this returns the corresponding correlation matrix.
#' @param m a covariance matrix
#' @example: par_read("sea16.std")
#' @value A list with components est and std for estimates and standard deviations.
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
#' @note There is a similar function in base R or a similar name, but Hans still calls this in places.
cov2corr <- function(m) diag(1/sqrt(diag(m))) %*% m %*% diag(1/sqrt(diag(m)))


#' test each element of x to see if it is in y
#'
#' A function that delivers the same functionality as x %in% y
#'
#' @param x a vector
#' @param y another vector
#' @example: member(1:5, 2:10)
#' @value A logical vector of the same length as x
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
#' @note There is similar functionality in base R, but Hans calls this in places so we include it here.
member <- function(x,y) !is.na(match(x,y))  


#' return logical vector picking out the lower triangle of an n x n square matrix
#'
#' A function that delivers the same functionality as lower.tri(matrix(1,n,n))
#'
#' @param n dimension of a square matrix
#' @param strictly  TRUE means do not include the diagonal and FALSE means do.
#' @example: below(10)
#' @value A logical vector (matrix) with TRUE on the lower diagonal (and on the diagonal with strictly=F)
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no}
#' @note There is similar functionality in base R, but Hans calls this in places so we include it here.
below <- function(n,strictly=F)
{
  M <- matrix(T,n,n)
  M[rep(1:n,n)<rep(1:n,rep(n,n))] <- F
  if(strictly)
    diag(M) = F
  M
}

#' read Hessian of dimension n from admb .cor file
#'
#' 
#' @param file path to file holding the Hessian
#' @param n number of rows and columns of the Hessian
#' @param cor logical indicating whether to return...
#' @example: readH("boing.cor", 10, T)
#' @value A matrix
#' @author Hans J. Skaug \email{Hans.Skaug@@math.uib.no} 
readH <- function(file,n,cor=F)
{
  N = n*(n+1)/2+4*n
  tmp = scan(file,what="",skip=2,quiet=T)
  if(length(tmp)<N) stop("n is too large")
  tmp = tmp[1:N]
  
  stdtab = numeric(n)
  H = diag(n)
  
  for(i in 1:n)
  {
    stdtab[i] = as.numeric(tmp[4])
    tmp = tmp[-(1:4)]
    H[i,1:i] = as.numeric(tmp[1:i])
    tmp = tmp[-(1:i)]
  }
  
  if(length(tmp)!=0)
  {
    print(length(tmp))
    stop("Det er noe galt")
  }
  
  H = H+t(H)	# Fill in upper diagonal
  diag(H) = 1
  
  
  if(!cor)
    H = diag(stdtab) %*% H %*% diag(stdtab)
  
  H
}

