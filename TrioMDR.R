
require(MASS)
require(compiler)
require(glmnet)

## read trio data
readData <- function(fileName){
  
  dat <- read.table(fileName,sep=" ",header = F)
  nc <- ncol(dat)
  nSnp <- (nc-6)/2 
  
  hap.dat <- dat[,7:nc]
  a <- seq(from=1,to=2*nSnp-1,by=2)
  b <- 2*c(1:nSnp)
  geno.dat <- hap.dat[,a]+hap.dat[,b]
  geno.dat <- geno.dat-2
  return(geno.dat)
}

## partition genotype combinations into high-risk group or low-risk group based on geno-PDT statistics
class_HL <- function(genodata,classm="PDT"){
  
  genodata<- as.matrix(genodata)
  nTrio <-  nrow(genodata)/3
  k = ncol(genodata)
  n = 2*nTrio
  
  dat <- matrix(data=NA,nrow= 2*nTrio,ncol=k)
  
  # transmitted alleles
  dat[1:nTrio,] <- genodata[3*c(1:nTrio),]
  
  # Nontransmitted alleles
  dat[(nTrio+1):(2*nTrio),] <- genodata[(3*c(1:nTrio)-2),]+genodata[(3*c(1:nTrio)-1),] - genodata[3*c(1:nTrio),]

  ids = which(complete.cases(dat))
  tlist = vector('list', k)
  i <- 1
  for(i in 1:k) tlist[[i]] = dat[, i]
  cells = split(data.frame(cbind(ids, dat)), tlist)
  
  # delete NULL cells
  obs.cell = sapply(cells, function(x) nrow(x))
  cell.null = which(obs.cell == 0)
  if(length(cell.null) > 0 ) cells = cells[-cell.null]
  
  ncells = length(cells)
  PDTstatistics <- rep(0,ncells)
  clabels <- rep(0,nTrio)
  plabels <- rep(0,nTrio)
  threshold <- 0
  
  for(i in 1:ncells){
    ctrl.idxs <- cells[[i]]$ids[which(cells[[i]]$ids>nTrio)]-nTrio
    case.idxs <- cells[[i]]$ids[which(cells[[i]]$ids<=nTrio)]
    
    a <- length(setdiff(case.idxs,ctrl.idxs))
    b <- length(setdiff(ctrl.idxs,case.idxs))
    if( a==0 && b==0 ){
      PDTstatistics[i] <- 0
    }else {
      PDTstatistics[i] <- (a-b)/sqrt(a+b)
      if (PDTstatistics[i]>threshold){
	    # labels indicate high/low-risk, 1 for high-risk and 0 for low-risk
        clabels[case.idxs] <- 1
        plabels[ctrl.idxs] <- 1
      }
    }
  } 
  labels <- c(clabels,plabels)
  return(list("labels"=labels,"cdat"=dat))
}
class_HL <- cmpfun(class_HL)

## generate permutated data based on Mendelian inheritance law
permuData <- function(genodata){
  
  nloci <- ncol(genodata)
  nTrio <- nrow(genodata)/3 
  perm.genodata <- genodata
  i <- 1
  for(i in 1:nloci){
    parent1 <- genodata[(3*c(1:nTrio)-2),i]
    parent2 <- genodata[(3*c(1:nTrio)-1),i]
    offspring <- rep(0,nTrio)
    dat <- cbind(parent1,parent2)
    
    ids = which(complete.cases(dat))
    tlist = vector('list', 2)
    j <- 1
    for(j in 1:2) tlist[[j]] = dat[, j]
    cells = split(data.frame(cbind(ids, dat)), tlist)
    
    # delete NULL cells
    obs.cell = sapply(cells, function(x) nrow(x))
    cell.null = which(obs.cell == 0)
    if(length(cell.null) > 0 ) cells = cells[-cell.null]
    
    ncells = length(cells)
    k <- 1
    for(k in 1:ncells){
      genotypes <- strsplit(names(cells[k]),split='[.]')
      genotypes <- as.numeric(genotypes[[1]])
      tmp <- sum(genotypes)
      index <- cells[[k]]$ids
      count <- length(index)
      probs <- runif(count,0,1)
	  
	  # parents have same homozygous genotypes
      if(tmp==0 || tmp==4){
        offspring[index] <- genotypes[1]
        next
      }
	  
	  #parents have a homozygous genotype and a heterozygous genotype
      if(tmp==1){
        offspring[index[which(probs>=0.5)]] <- 1
        offspring[index[which(probs< 0.5)]] <- 0
        next
      }
	  
	  # parents have same heterozygous genotypes
      if(tmp==2){
        if(genotypes[1]==1){
          offspring[index[which(probs< 0.25)]] <- 0
          offspring[index[which(probs>=0.25 & probs < 0.75)]] <- 1
          offspring[index[which(probs>=0.75)]] <- 2
        }else{
          offspring[index] <- 1
        }
        next
      }
	  
	  #parents have a homozygous genotype and a heterozygous genotype
      if(tmp==3){
        offspring[index[which(probs>=0.5)]] <- 2
        offspring[index[which(probs < 0.5)]] <- 1
      }
    }
    perm.genodata[3*(1:nTrio),i] <- offspring
  } 
  return(perm.genodata)
}
permuData <- cmpfun(permuData)


## estimat non-center parameter of the null distribution for the Wald type statistics
est_nullcenter <- function(phe, snp_pair, classm = 'PDT', adj.main='FALSE', lbd = 0, nperm = 10){
  if (exists(".Random.seed", .GlobalEnv)){
    oldseed <- .GlobalEnv$.Random.seed
  }else {
    oldseed <- NULL
  }
  
  set.seed(12345)
  ncase = nrow(snp_pair)/3
  centers = rep(0, nperm)
  coefs = rep(0, nperm)
  
  for(i in 1:nperm){
    
    permuted.pair <- permuData(snp_pair)
    res = class_HL(permuted.pair, classm)
    labels <- res$labels
    cdat <- res$cdat
    if(adj.main){
      xx <- as.matrix(cbind(cdat,labels))
    }else{
      xx =labels
    }
    if(lbd == 0){
      sumstat = coef(summary(glm(phe ~ xx, family ='binomial')))
      centers[i] = sumstat[nrow(sumstat), 3]
      coefs[i] = sumstat[nrow(sumstat), 1]
    }
  }
  if (!is.null(oldseed))
    .GlobalEnv$.Random.seed <- oldseed
  else
    rm(".Random.seed", envir = .GlobalEnv)
  
  return(list('w_hat' = centers^2, 'beta' = coefs))
}
est_nullcenter = cmpfun(est_nullcenter)

##Detect SNP-SNP interaction by TrioMDR
#' @param geno.dat, snp matrix, n by p.
#' @param K, K-locus interaction, default 2.
#' @param classm classification rule to define H/L, PDT indicates geno-PDT statistic
#' @param adj.main logical value, adjust marginal effect or not, default FALSE.
#' @param nperm small number of permutation time to estimate non-central parameters, default 5
#' @return a list including the elements:
#' \item{beta}{coefficients.}
#' \item{pvs}{p-values corrected by semiparamtric correction procedure.}
#' \item{raw.pvs}{raw p-values.}
#' \item{best_comb}{best k-snps, which gives smallest corrected p-value.}

TrioMDR <- function(geno.dat, K=2, classm='PDT', adj.main='FALSE', nperm=5){

  nloci <- ncol(geno.dat)
  snp.combs <- combn(nloci, K)
  ns <-  ncol(snp.combs)
  ncase <- nrow(geno.dat)/3
  phe <- rep(0,2*ncase)
  phe[1:ncase] <- 1
  phe <- as.matrix(phe)
  
  # save all coefs and pvs for S
  rpvs = pvs = coefs = rep(0.5, ns)
  
  # select best model(i.e. snp combination)
  # use glm
  est_w = stats =  stats0 = NULL
  
  for(j in 1:ns){
    res <- class_HL(geno.dat[,snp.combs[,j]],classm)
    labels <- res$labels
    cdat <- res$cdat
    if(adj.main){
      xx<- as.matrix(cdat)
      labels <- cbind(cdat,labels)
    }
    mod = glm(phe ~ labels, family ='binomial')
    res = coef(summary(mod))
    pp = nrow(res)
    coefs[j] = res[pp, 1]  # the coefficent vector
    rpvs[j] = res[pp, ncol(res)]
    vv = res[pp, 2]^2
    stat = res[pp, 1]^2/vv
    
    tmp2 = est_nullcenter(phe, geno.dat[,snp.combs[,j]], classm, adj.main, lbd = 0, nperm)
    
    lambda = mean(tmp2$w_hat) - 1
    if(lambda >= 0) pvs[j] = pchisq(stat, df = 1, ncp = lambda, lower.tail = FALSE )
    
    stat0 = abs(coefs[j] - mean(tmp2$beta))/sqrt(vv)
    if(lambda < 0) pvs[j] = 2 * pnorm(-stat0)
    est_w = c(est_w, tmp2$w_hat)
    stats = c(stats, stat)
    stats0 = c(stats0, stat0)
  }

  best_comb = snp.combs[, which.min(pvs)]
  return(list('beta' = coefs, 'pvs' = pvs, 'raw.pvs' = rpvs, 'best_comb' = best_comb))
}
TrioMDR <- cmpfun(TrioMDR)