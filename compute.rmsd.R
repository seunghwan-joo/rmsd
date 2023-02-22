# Dichotomous vs Polytomous
dich <- T

# Number of Items
J <- 40

# Sample Size
N <- 500

# Number of Groups
G <- 80


# 2PL function
twopl <- function(a,b,theta){
  prob=exp(1.7*(a*(theta-b)))/(1+exp(1.7*(a*(theta-b))))
  return(prob)
}

# GPCM function 
gpcm <- function(a,b,d,theta,cat){
  num <- rep(0,cat)
  for(i in 1:cat){
    if(i==1){
      z <- 1.7*a*(theta-b)
      num[i] <- exp(z)
    }
    if(i>1){ 
      z <- z + 1.7*a*(theta-b+d[i-1])
      num[i] <- exp(z)
    }
  }
  deno <- sum(num)
  p <- num/deno
  return(p)
}

# 2PL information function
twoplinfo <- function(a,b,theta){
  D=1.7
  p=twopl(a,b,theta)
  info=D^2*a^2*p*(1-p)
  return(info)
}

# GPCM information function
gpcminfo <- function(a,b,d,theta,cat){
  p <- gpcm(a,b,d,theta,cat)
  tmp1 <- rep(0,cat)
  tmp2 <- rep(0,cat)
  for(i in 1:cat){
    tmp1[i] <- (i-1)^2*p[i] 
    tmp2[i] <- (i-1)*p[i]
  }
  info <- 1.7^2*a^2*(sum(tmp1)-(sum(tmp2))^2)
  return(info)
}

# Read in data file
datafile <- "SimData.dat"
data <- read.fwf(datafile,c(3,rep(1,J)))

# Read in lst file
lstfile <- "mdltm.lst"
lst <- read.mdltm(lstfile)

# Item Parameters
item <- lst$pars
if(dich==T) itempar <- item[,c("Group","Item","Slope","Difficulty","Cats")]
if(dich!=T) itempar <- item[,c("Group","Item","Slope","Difficulty","IRT_Step1","IRT_Step2","Cats")]
itempar$Cats <- itempar$Cats - 1

# Group Means and SDs
group <- lst$mv
grpmn <- group[,"Mean"]
grpsd <- group[,"SD"]

# Compute group weights
n.grp <- G
n.item <- J
qt <- seq(-4,4,by=.2)
n.qt <- length(qt)
grpwt <- matrix(NA,nr=n.qt,nc=n.grp)
for(g in 1:n.grp){
  for(q in 1:n.qt){
    grpwt[q,g] <- dnorm(qt[q],grpmn[g],grpsd[g])/sum(dnorm(qt,grpmn[g],grpsd[g]))
  }
}

# Compute item information (normalized)
info <- array(NA,dim=c(n.qt,n.item,n.grp))
for(g in 1:n.grp){
  ipar <- itempar[itempar[,"Group"]==100+g,]
  for(i in 1:n.item){
    for(q in 1:n.qt){
      if(ipar$Cats[i]==2){
        info_sum <- matrix(0,nr=n.item,nc=n.grp)
        for(qq in 1:n.qt){
          info_sum[i,g] <- info_sum[i,g] + twoplinfo(ipar$Slope[i],ipar$Difficulty[i],qt[qq])
        }
        info[q,i,g] <- twoplinfo(ipar$Slope[i],ipar$Difficulty[i],qt[q])/info_sum[i,g]
        #info[q,i,g] <- info_tmp[q,i,g]/sum(twoplinfo(ipar$Slope[i],ipar$Difficulty[i],qt))
      }
      if(ipar$Cats[i]>2){
        info_sum <- matrix(0,nr=n.item,nc=n.grp)
        for(qq in 1:n.qt){
          info_sum[i,g] <- info_sum[i,g] + gpcminfo(ipar$Slope[i],ipar$Difficulty[i],c(ipar$IRT_Step1[i],ipar$IRT_Step2[i]),qt[qq],ipar$Cats[i])
        }
        info[q,i,g] <- gpcminfo(ipar$Slope[i],ipar$Difficulty[i],c(ipar$IRT_Step1[i],ipar$IRT_Step2[i]),qt[q],ipar$Cats[i])/info_sum[i,g]
      }
    }
  }
}

# Compute posterior and observed icc 
obs <- array(NA,dim=c(n.qt,n.item,n.grp))
for(g in 1:n.grp){
  ipar <- itempar[itempar[,"Group"]==100+g,]
  gdata <- data[data[,1]==100+g,2:(n.item+1)]

  itemtrace <- matrix(NA,nr=n.qt,nc=n.item)
  itemtrace0 <- matrix(NA,nr=n.qt,nc=n.item)
  itemtrace1 <- matrix(NA,nr=n.qt,nc=n.item)
  itemtrace2 <- matrix(NA,nr=n.qt,nc=n.item)
  pseudo <- matrix(0,nr=n.qt,nc=n.item)
  deno <- matrix(0,nr=n.qt,nc=n.item)
  for(i in 1:n.item){
    for(q in 1:n.qt){
      if(ipar$Cats[i]==2){
        itemtrace[q,i] <- twopl(ipar$Slope[i],ipar$Difficulty[i],qt[q])
      }
      if(ipar$Cats[i]>2){
        tmp.prob <- gpcm(ipar$Slope[i],ipar$Difficulty[i],c(ipar$IRT_Step1[i],ipar$IRT_Step2[i]),qt[q],ipar$Cats[i])
        itemtrace0[q,i] <- tmp.prob[1]
        itemtrace1[q,i] <- tmp.prob[2]
        itemtrace2[q,i] <- tmp.prob[3]
      }
    }
  }
  n.ex <- nrow(gdata)
  for(p in 1:n.ex){
    likelihood <- rep(1,n.qt)
    for(i in 1:n.item){
      if(ipar$Cats[i]==2){
        if(gdata[p,i]==1) likelihood <- likelihood*itemtrace[,i]
        else likelihood <- likelihood*(1-itemtrace[,i])
      }
      if(ipar$Cats[i]>2){
        if(gdata[p,i]==2) likelihood <- likelihood*itemtrace2[,i]
        if(gdata[p,i]==1) likelihood <- likelihood*itemtrace1[,i]
        if(gdata[p,i]==0) likelihood <- likelihood*itemtrace0[,i]
      }
    }
    posterior <- likelihood*grpwt[,g]
    
    # normalize posterior 
    expd <- sum(posterior)
    posterior <- (posterior/expd)
    
    # put this response pattern (pseudo counts)
    for(i in 1:n.item){
      deno[,i] <- deno[,i] + posterior
      if(ipar$Cats[i]==2){
        if(gdata[p,i]==1) pseudo[,i] <- pseudo[,i] + posterior
      }
      if(ipar$Cats[i]>2){
        if(gdata[p,i]==2) pseudo[,i] <- pseudo[,i] + posterior
        if(gdata[p,i]==1) pseudo[,i] <- pseudo[,i] + 0.5*posterior
      }
    }
  }
  obs[,,g] <- pseudo/deno
}

# Compute expected icc
exp <- array(NA,dim=c(n.qt,n.item,n.grp))
for(g in 1:n.grp){
  ipar <- itempar[itempar[,"Group"]==100+g,]
  for(i in 1:n.item){
    if(ipar$Cats[i]==2){
      for(q in 1:n.qt){
        exp[q,i,g] <- twopl(ipar$Slope[i],ipar$Difficulty[i],qt[q])
      }
    }
    if(ipar$Cats[i]>2){
      for(q in 1:n.qt){
        tmp.p <- gpcm(ipar$Slope[i],ipar$Difficulty[i],c(ipar$IRT_Step1[i],ipar$IRT_Step2[i]),qt[q],ipar$Cats[i])
        exp[q,i,g] <- 0.5*tmp.p[2] + tmp.p[3]
      }
    }
  }
}

# Compute MD and RMSDs
md <- matrix(NA,nr=n.item,nc=n.grp)
rmsd <- matrix(NA,nr=n.item,nc=n.grp)
rmsd_eqw <- matrix(NA,nr=n.item,nc=n.grp)
rmsd_infow <- matrix(NA,nr=n.item,nc=n.grp)
rmsd_bnorm <- matrix(NA,nr=n.item,nc=n.grp)
for(g in 1:n.grp){
  ipar <- itempar[itempar[,"Group"]==100+g,]
  grp.range.ind <- which(round(grpwt[,g],2) %in% 0.01)
  grp.range.length <- length(grp.range.ind)
  minind <- grp.range.ind[1]
  maxind <- grp.range.ind[grp.range.length]
  for(i in 1:n.item){
    md[i,g] <- sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])*grpwt[minind:maxind,g])
    rmsd[i,g] <- sqrt(sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])^2*grpwt[minind:maxind,g]))
    rmsd_eqw[i,g] <- sqrt(sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])^2/length(minind:maxind)))
    rmsd_infow[i,g] <- sqrt(sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])^2*info[minind:maxind,i,g]))
    rmsd_bnorm[i,g] <- sqrt(sum((obs[minind:maxind,i,g]-exp[minind:maxind,i,g])^2*dnorm(qt[minind:maxind],ipar$Difficulty[i],1)/sum(dnorm(qt[minind:maxind],ipar$Difficulty[i],1))))
  }
}

# MD
round(md,4)
# RMSD
round(rmsd,4)
# RMSD w/ equal weight
round(rmsd_eqw,4)
# RMSD w/ item information weight
round(rmsd_infow,4)
# RMSD w/ b-norm weight
round(rmsd_bnorm,4)

p.rmsd <- 0
p.rmsd_eqw <- 0
p.rmsd_infow <- 0
p.rmsd_bnorm <- 0
for(i in 1:n.item){
  for(j in 1:n.grp){
    if(rmsd[i,j] > .12) p.rmsd <- p.rmsd + 1
    if(rmsd_eqw[i,j] >.12) p.rmsd_eqw <- p.rmsd_eqw + 1
    if(rmsd_infow[i,j] >.12) p.rmsd_infow <- p.rmsd_infow + 1
    if(rmsd_bnorm[i,j] >.12) p.rmsd_bnorm <- p.rmsd_bnorm + 1
  }
}
p.rmsd
p.rmsd_eqw
p.rmsd_infow
p.rmsd_bnorm
