
# count mismatch between 2 samples, in percet ##
count_n_loci <- function(i,j,data_genotypes_N_combined){
  
  if(i!=j & i<j){ # only needs to compare i!=j & i<j 
    s1 = data_genotypes_N_combined[i,]
    s2 = data_genotypes_N_combined[j,]
    
    locus_position = which(is.na(s1)==FALSE & is.na(s2)==FALSE)
    n_12 = length(locus_position)
    # m_12 = sum(s1!=s2,na.rm=T)
    
    # mismatch_perct = round(m_12/n_12,2)
  }
  
  else{
    # mismatch_perct = NA
    n_12 = NA
    
  }
}

# count mismatch between 2 samples, in percet ##
count_mismatch <- function(i,j,data_genotypes_N_combined){
  
  if(i!=j & i<j){ # only needs to compare i!=j & i<j 
    s1 = data_genotypes_N_combined[i,]
    s2 = data_genotypes_N_combined[j,]
    
    locus_position = which(is.na(s1)==FALSE & is.na(s2)==FALSE)
    # n_12 = length(locus_position)
    m_12 = sum(s1!=s2,na.rm=T)
    
    # mismatch_perct = round(m_12/n_12,2)
  }
  
  else{
    # mismatch_perct = NA
    m_12 = NA
    
  }
}

## function to generate P_G (relatively arbitrarily), the probs of each element in genotypes_pool ##
generate_p_G <- function(ind_major_pairs, genotypes_pool, which_major_pairs, p_G_major_pairs, p_other_pairs){
  
  p_G = rep(NA, length(genotypes_pool))
  p_G[which_major_pairs[[ind_major_pairs]]] = p_G_major_pairs
  
  y = runif(length(genotypes_pool)-3,1,100)
  p_G[-which_major_pairs[[ind_major_pairs]]] = (y/sum(y))*(p_other_pairs) # assign probs for remaining possible (non-major) allele pairs
  
  p_G
  
}  

## function to return a sample incorporating genotyping error ##
generate_geno_error <- function(sample,rho, genotypes_pool){
  
  n = length(sample)
  
  change_or_not = rbinom(n,size=1,prob=rho)
  
  
  sample_with_geno_error = mapply(function(sample_i,change_or_not_i,genotypes_pool){
    
    ifelse(change_or_not_i==0,sample_i,sample(genotypes_pool[-sample_i],1))
    
  }, sample, change_or_not, MoreArgs = list(genotypes_pool=genotypes_pool))
  
  sample_with_geno_error
  
}


# ## function to return Pr(gA,gB) on a locus for same-individual samples ##
compute_pr_gAgB_II <- function(q_and_g,rho,probs_pool_locus){ # qA and qB (in q) are the genotype freqs of gA and gB in the cohort population
  
  
  qA=q_and_g[1]
  qB=q_and_g[2]
  
  gA=q_and_g[3]
  gB=q_and_g[4]
  
  
  n_probs = length(probs_pool_locus)
  probs_pool_locus_sub = probs_pool_locus[-c(gA,gB)]
  
  
  (2*( (1/(n_probs-1))*(1-rho)*rho*(qB+qA) + (1/(n_probs-1))^2*(rho^2)*sum(probs_pool_locus_sub) ))
  
}

## function to return Pr(gA,gB) for gA anf gB on a locus from different individuals##
compute_pr_gAgB_IJ <- function(q,rho,n_probs){ # qA and qB (in q) are the genotype freqs of gA and gB in the cohort population
  
  qA=q[1]
  qB=q[2]
  
  (2*( qA*qB*(1-rho)^2 + (1/(n_probs-1))*((1-qA)*qB + qA*(1-qB))*(1-rho)*rho + (1/(n_probs-1))^2*(1-qA)*(1-qB)*(rho^2))) # as of 4 Jan 2019
  
  
  # qA*qB*(1-rho)^2 + 
  #   (1/(n_probs-1))*((1-qA)*qB + qA*(1-qB))*(1-rho)*rho + 
  #   (1/(n_probs-1))^2*(1-qA)*(1-qB)*(rho^2) +
  #   qA^2*rho*(1-rho)*(1/(n_probs-1)) + 
  #   qB^2*rho*(1-rho)*(1/(n_probs-1))+
  #   (1-qA-qB)^2*(rho^2)*(1/(n_probs-1))^2 
  
  # (2*( qA*qB*(1-rho)^2 + (1/(n_probs-1))*((1-qA)*qB + qA*(1-qB))*(1-rho)*rho ))
  
  # 2*(qA*qB*(1-rho)^2 + qA^2*rho*(1-rho)*(1/(n_probs-1)) + qB^2*rho*(1-rho)*(1/(n_probs-1)) )
  
}



## function to return Pr(gA,gB) for gA anf gB on a locus from different individuals (Germline vs Germline, or Tunor vs Tumor)##
# compute_pr_gAgB_IJ_same_type <- function(q,rho,n_probs){ # qA and qB (in q) are the genotype freqs of gA and gB in the cohort population
#   
#   qA=q[1]
#   qB=q[2]
#   
#   (2*( qA*qB*(1-rho)^2 + (1/(n_probs-1))*((1-qA)*qB + qA*(1-qB))*(1-rho)*rho + (1/(n_probs-1))^2*(1-qA)*(1-qB)*(rho^2)))
#   
# }
# 
# 
# ## function to return Pr(gA,gB) for gA anf gB on a locus from different individuals (Germline vs Tumor, germline places at 1st)##
# compute_pr_gAgB_IJ_mixed_type <- function(q,rho_GL, rho_TM,n_probs){ # qA and qB (in q) are the genotype freqs of gA and gB in the cohort population
#   
#   qA=q[1] # on Germline sample
#   qB=q[2] # on Tumor sample
#   
#   # (2*( qA*qB*(1-rho)^2 + (1/(n_probs-1))*((1-qA)*qB + qA*(1-qB))*(1-rho)*rho + (1/(n_probs-1))^2*(1-qA)*(1-qB)*(rho^2)))
# 
#   Pr_AB = qA*qB*(1-rho_GL)*(1-rho_TM) + qA*(1-qB)*(1/(n_probs-1))*(1-rho_GL)*rho_TM +  (1-qA)*qB*(1/(n_probs-1))*rho_GL*(1-rho_TM) + (1-qA)*(1-qB)*(1/(n_probs-1))^2*(rho_GL*rho_TM)
#   Pr_BA = qB*qA*(1-rho_GL)*(1-rho_TM) + qB*(1-qA)*(1/(n_probs-1))*(1-rho_GL)*rho_TM +  (1-qB)*qA*(1/(n_probs-1))*rho_GL*(1-rho_TM) + (1-qB)*(1-qA)*(1/(n_probs-1))^2*(rho_GL*rho_TM)
# 
#   Pr = Pr_AB + Pr_BA
# 
# }
# 

## model for mismatch between samples from same individual, a function of rho ##
model_II <- function(rho,genotypes_pool,probs_pool_matrix,n_II,m_II){
  
  # pG_II =sum(apply(rbind(combn(probs_pool,2),combn(genotypes_pool,2)),MARGIN=2,compute_pr_gAgB_II,rho=rho,probs_pool=probs_pool))
  # dbinom(m_II,n_II,pG_II)
  
  n_bp= ncol(probs_pool_matrix)
  
  pG_II = sapply(1:n_bp, FUN=function(i){sum(apply(rbind(combn(probs_pool_matrix[,i],2),combn(genotypes_pool,2)),MARGIN=2,compute_pr_gAgB_II,rho=rho,probs_pool_locus=probs_pool_matrix[,i]))}) # a vector 
  
  max(0,dpoisbinom(m_II,pp=pG_II, log_d = FALSE))
  
}

## model for mismatch between samples from different individuals, a function of rho ##
model_IJ <- function(rho,genotypes_pool,probs_pool_matrix,n_IJ,m_IJ){
  
  # pG_IJ =sum(apply(combn(probs_pool,2),MARGIN=2,compute_pr_gAgB_IJ,rho=rho,n_probs=length(probs_pool))) 
  # dbinom(m_IJ,n_IJ,pG_IJ)
  
  n_bp= ncol(probs_pool_matrix)
  
  pG_IJ = sapply(1:n_bp, FUN=function(i){sum(apply(combn(probs_pool_matrix[,i],2),MARGIN=2,compute_pr_gAgB_IJ,rho=rho,n_probs=nrow(probs_pool_matrix)))})
  
  # max(0,dpoisbinom(m_IJ,pp=pG_IJ, log_d = FALSE))
  
  dpoibin (m_IJ,pp=pG_IJ)
  
}


##vectorized model_II ##
model_II_vectorized <- function(rho,genotypes_pool,probs_pool_matrix,n_II,m_II){
  
  # pG_II =sum(apply(rbind(combn(probs_pool,2),combn(genotypes_pool,2)),MARGIN=2,compute_pr_gAgB_II,rho=rho,probs_pool=probs_pool))
  # dbinom(m_II,n_II,pG_II)
  
  n_bp= ncol(probs_pool_matrix)
  
  # pG_II = sapply(1:n_bp, FUN=function(i){rowSums(apply(rbind(combn(probs_pool_matrix[,i],2),combn(genotypes_pool,2)),MARGIN=2,compute_pr_gAgB_II,rho=rho,probs_pool_locus=probs_pool_matrix[,i]))}) # a vector 
  # sapply(1:nrow(pG_II),function(i){dpoibin(m_II,pp=pG_II[i,])} ) # as of Jan 4, 2019
  
  ### sapply(1:nrow(pG_II),function(i){max(0,dpoisbinom(m_II,pp=pG_II[i,], log_d = FALSE))} )
  
  
  sapply(1:length(rho),function(i){dbinom(m_II,n_II, prob=(1-(1-rho[i])^2) )} )
  
}



## vectorized modle IJ ##
model_IJ_vectorized <- function(rho,genotypes_pool,probs_pool_matrix,n_IJ,m_IJ){
  
  n_bp= ncol(probs_pool_matrix)
  
  pG_IJ = sapply(1:n_bp, FUN=function(i){rowSums(apply(combn(probs_pool_matrix[,i],2),MARGIN=2,compute_pr_gAgB_IJ,rho=rho,n_probs=nrow(probs_pool_matrix)))})
  
  
  # pG_IJ = sapply(1:n_bp, FUN=function(i){ # when considering only heterogeneous genotypes (3 possible combinations)
  #   tmp = combn(unique(probs_pool_matrix[,i][which(probs_pool_matrix[,i]!=0)]),2)
  #   rowSums(apply(tmp,MARGIN=2,compute_pr_gAgB_IJ,rho=rho,n_probs=3))
  # })
  
  
  # sapply(1:nrow(pG_IJ),function(i){max(0,dpoisbinom(m_IJ,pp=pG_IJ[i,]))} )
  
  sapply(1:nrow(pG_IJ),function(i){dpoibin(m_IJ,pp=pG_IJ[i,])} )
  
}


# model_II_vectorized  <- Vectorize(model_II,vectorize.args = "rho")
# model_IJ_vectorized <- Vectorize(model_IJ,vectorize.args = "rho")

## compute BF (H0 vs H1), for sample i and j ##
generate_BF <- function(i,j,data_genotypes_N_combined,data_reads_N_combined,genotypes_pool,probs_pool_matrix,reads_cut_off,rho_levels, rho_mid_pts, n_rhos,shape_1,shape_2, max_rho){ 
  
  if(i!=j & i<j){ # only needs to compare i!=j & i<j 
    gA = data_genotypes_N_combined[i,]
    gB = data_genotypes_N_combined[j,]
    
    locus_position = which(is.na(gA)==FALSE & is.na(gB)==FALSE)
    
    if(length(locus_position)>=20){ 
      
      gA= gA[locus_position]
      gB = gB[locus_position]
      
      if(is.na(reads_cut_off)==FALSE){ # 
        rA = data_reads_N_combined[i,]
        rB = data_reads_N_combined[j,]
        n_AB = length(which(rA>reads_cut_off&rB>reads_cut_off))
        m_AB = sum(gA[which(rA>reads_cut_off&rB>reads_cut_off)]!=gB[which(rA>reads_cut_off&rB>reads_cut_off)])
      }
      
      
      else{
        n_AB = length(gA)
        m_AB = sum(gA!=gB)
        
      }
      
      
      # l_AB_H0 = integrate(model_II_vectorized,lower=0,upper=1,genotypes_pool=genotypes_pool,probs_pool_matrix=probs_pool_matrix,n_II=n_AB,m_II=m_AB)[["value"]]
      # l_AB_H1 = integrate(model_IJ_vectorized,lower=0,upper=1,genotypes_pool=genotypes_pool,probs_pool_matrix=probs_pool_matrix,n_IJ=n_AB,m_IJ=m_AB)[["value"]] # still use n_12 and m_12 as the data have not changed
      # BF_AB = l_AB_H0/l_AB_H1
      
      # l_AB_H0 = sum(diff(rho_levels)*sapply(rho_mid_pts,model_II,genotypes_pool=genotypes_pool,probs_pool_matrix=probs_pool_matrix,n_II=n_AB,m_II=m_AB))
      # l_AB_H1 = sum(diff(rho_levels)*sapply(rho_mid_pts,model_IJ,genotypes_pool=genotypes_pool,probs_pool_matrix=probs_pool_matrix,n_IJ=n_AB,m_IJ=m_AB)) # still use n_12 and m_12 as the data have not changed
      # BF_AB = l_AB_H0/l_AB_H1
      
      l_AB_H0 = sum(diff(rho_levels)*model_II_vectorized(rho_mid_pts,genotypes_pool=genotypes_pool,probs_pool_matrix=probs_pool_matrix[,locus_position],n_II=n_AB,m_II=m_AB))
      l_AB_H1 = sum(diff(rho_levels)*model_IJ_vectorized(rho_mid_pts,genotypes_pool=genotypes_pool,probs_pool_matrix=probs_pool_matrix[,locus_position],n_IJ=n_AB,m_IJ=m_AB)) # still use n_12 and m_12 as the data have not changed
      BF_AB = l_AB_H0/l_AB_H1
      
      # l_AB_H0 = (1/n_rhos)*sum(model_II_vectorized(rbeta(rtrunc(n_rhos,"beta",a=0,b=max_rho,shape1=shape_1,shape2=shape_2),shape_1,shape_2),genotypes_pool=genotypes_pool,probs_pool_matrix=probs_pool_matrix[,locus_position],n_II=n_AB,m_II=m_AB))
      # l_AB_H1 = (1/n_rhos)*sum(model_IJ_vectorized(rbeta(rtrunc(n_rhos,"beta",a=0,b=max_rho,shape1=shape_1,shape2=shape_2),shape_1,shape_2),genotypes_pool=genotypes_pool,probs_pool_matrix=probs_pool_matrix[,locus_position],n_IJ=n_AB,m_IJ=m_AB)) # still use n_12 and m_12 as the data have not changed
      # BF_AB = l_AB_H0/l_AB_H1
      
    } else{
      BF_AB=NA
    }

    
  } else{
    
    BF_AB = NA
  }
  
  # ind = ifelse(ind_individual[i]==ind_individual[j],1,0)
  # c(BF_AB,ind)
  
  BF_AB
  
}

### unused codes ###


##### explore Likelihood ratio test #####

# pG =  sum(apply(combn(probs_pool,2),MARGIN=2,compute_pr_gAgB,rho=rho)) # pr of observing a mismatch when the pair does not come from same individual
# 
# s1 = data[[1]]$genotypes[1,]
# s2 = data[[1]]$genotypes[2,]
# s3 = data[[5]]$genotypes[4,]
# 
# r1 = data[[1]]$reads[1,]
# r2 = data[[1]]$reads[2,]
# r3 = data[[5]]$reads[4,]
# 
# n_12 = length(which(r1>reads_cut_off&r2>reads_cut_off))
# m_12 = sum(s1[which(r1>reads_cut_off&r2>reads_cut_off)]!=s2[which(r1>reads_cut_off&r2>reads_cut_off)])
# 
# n_13 = length(which(r1>reads_cut_off&r3>reads_cut_off))
# m_13 = sum(s1[which(r1>reads_cut_off&r3>reads_cut_off)]!=s3[which(r1>reads_cut_off&r3>reads_cut_off)])
# 
# 
# dist_LRT_12 = log(dbinom(0:n_12,n_12,2*rho)) - log(dbinom(0:n_12,n_12,pG))
# LRT_12 =  log(dbinom(m_12,n_12,2*rho)) - log(dbinom(m_12,n_12,pG))
# 
# 
# dist_LRT_13 = log(dbinom(0:n_13,n_13,2*rho)) - log(dbinom(0:n_13,n_13,pG))
# LRT_13 =  log(dbinom(m_13,n_13,2*rho)) - log(dbinom(m_13,n_13,pG))


#### explore nimble ####
# 
# require(nimble)
# 
# mismatchCode <- nimbleCode({
#   
#   for (i in 1:length(m)){ 
#   
#       
#     if (z[i]==0) m[i] ~ dbinom(size=n_comp[i],prob=2*rho) # when corresponds to two samples from same individual (n_comp is the number of comparable snps)
#     
#     if(z[i]==1){
#       
#         pG <-  sum(apply(combn(probs_pool,2),MARGIN=2,compute_pr_gAgB,rho=rho)) # pr of observing a mismatch when the pair does not come from same individual
#         
#         m[i] ~ dbinom(size=n_comp[i],prob=pG)
#     }
#     
#   }
# 
#   rho ~ dunif(0,1) # prior?
# })
# 
# 
# n_comp=rep(50,10)
# m= c(rep(10,5),rep(30,5))
# z=rbinom(length(m),size=1,prob=0.8)
# 
# mismatchConstant <- list(n_comp=n_comp, probs_pool=probs_pool)
# mismatchData <- list(m=m)
# mismatchInits <- list(rho=rho,z=z)
# 
# ##
# 
# mismatch <- nimbleModel(code=mismatchCode,name="mismatch",constants=mismatchConstant,data=mismatchData,inits=mismatchInits)
# 
# 

