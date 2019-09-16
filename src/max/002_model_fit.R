### perform probability computations: most important output is "posterior_H0_matrix" which contains pairwise probabilities any two samples belong to a same individual ##
#########################################################################################################################################################################

start_time <- Sys.time()

##

set.seed(1)


source(file=paste(dir,"/code/R/src_functions.R",sep=""))

ind_subset_data = FALSE # only test on a subset of the full data, but may still use the full data to estimate the population genotype freqs
if (ind_subset_data==TRUE) n_subset_data = 200

ind_subset_locus = FALSE # only test on a subset of locus
threshold_major_allele_freq = 0.5 # threshold to exclude some loci

mc_cores = 32 # number of cores for parallelization (if used)
mc_cores = detectCores(all.tests = FALSE, logical = TRUE)

max_rho= 0.02 # max level of noise

n_rhos = 10 # number of rhos to sample in estimating the marginal likelihood in BF calculation
shape_1 = NA # shape1 parameter for the Beta dist, priror for rho
shape_2 = NA

# rho_levels = c(seq(0,0.1,0.001),seq(0.1,1,0.02)) # rho levels for integration for BF calculations; take more frequent subdivision at smaller rho
rho_levels = c(seq(0,max_rho,0.002)) # rho levels for integration for BF calculations; take more frequent subdivision at smaller rho
rho_mid_pts = rho_levels[-length(rho_levels)] + diff(rho_levels)/2

genotypes_pool = c("A/C","A/G","A/T","C/G","C/T","G/T", "A/A","C/C","G/G","T/T")

observed_genotypes = as.matrix(test_data[,sampleTosubject$sample_ID])
for (i in seq(ncol(observed_genotypes))){
  message("[", i, "/", ncol(observed_genotypes), "] Scanning observed genotypes")
  # observed_genotypes[,i] = as.numeric(sapply(observed_genotypes[,i] , function(x){min(which(genotypes_pool==x),which(genotypes_pool==stringi::stri_reverse(x)))})) # use numeric indexes to represent genotypes
  observed_genotypes[,i] = as.numeric(sapply(observed_genotypes[,i] , function(x){
    ifelse(is.na(x)==FALSE, min(which(genotypes_pool==x),which(genotypes_pool==stringi::stri_reverse(x))),NA)
  })) # use numeric indexes to represent genotypes (this version will contain NA)
}

sample_ID = sampleTosubject$sample_ID
subject_ID = sampleTosubject$subject_ID

probs_pool_matrix = as.matrix(apply(observed_genotypes,1,function(x){prop.table(table(factor(x,levels=1:10)))})) # matrix, each column contains the probs of elements in genotypes_pool for a position (number of alaways rows = 10, the number of possible genotypes)
colnames(probs_pool_matrix) = NULL

if(ind_subset_locus==TRUE) {
  # n_subset_locus = 1000
  # ind_subset_locus = sample(1:ncol(probs_pool_matrix),n_subset_locus)
  
  # ind_subset_locus = which(apply(probs_pool_matrix,2,function(x){sum(x>0)})==3)
  
  # tmp = which( (apply(probs_pool_matrix,2,function(x){min(x[x!=0])})>(threshold_major_allele_freq^2))  & (apply(probs_pool_matrix,2,function(x){sum(x>0)})==3) )
  tmp = which((apply(probs_pool_matrix,2,function(x){sum(x>0)})==3) )

  # tmp = which( (apply(probs_pool_matrix,2,function(x){max(x[x!=1])})>(threshold_major_allele_freq^2))  & (apply(probs_pool_matrix,2,function(x){sum(x>0)})==3) )
  # tmp = which( (apply(probs_pool_matrix,2,function(x){max(x)})>(threshold_major_allele_freq^2)) )
  
  observed_genotypes = observed_genotypes[tmp,]
  probs_pool_matrix = probs_pool_matrix[,tmp]
}

if(ind_subset_data==TRUE){
  # tmp = sample(1:ncol(observed_genotypes),n_subset_data)
  tmp = 1:n_subset_data
  observed_genotypes = observed_genotypes[,tmp]
  # probs_pool_matrix = as.matrix(apply(observed_genotypes,1,function(x){prop.table(table(factor(x,levels=1:10)))})) # matrix, each column contains the probs of elements in genotypes_pool for a position (number of alaways rows = 10, the number of possible genotypes)
  colnames(probs_pool_matrix) = NULL
  sample_ID= sampleTosubject$sample_ID[tmp]
  subject_ID = sampleTosubject$subject_ID[tmp]
}


n_size= table(factor(subject_ID, levels=unique(subject_ID)))
prior_H0 = sum(choose(n_size,2)) /choose(sum(n_size),2)
prior_H1 = 1 - prior_H0


###################################################################################################################


### Bayes factor comparing two hypotheses, H0= two samples come from a same individual, H1= two samples come from different individual ###

## decision rule based on 2*ln(BF), Kass and Raftery ##
# 0-2 not worth mentioning
# 2-6 positive
# 6-10 strong
# >10 very strong
# in coding, 1=negative, 2=not worth mentioning ,3=positive,4=strong, 5=very strong 

ID_1 = sample_ID[1]
ID_2 = sample_ID[3]

s1 = observed_genotypes[,ID_1]
s2 = observed_genotypes[,ID_2]

locus_position = which(is.na(s1)==FALSE & is.na(s2)==FALSE)
n_12 = length(locus_position)
m_12 = sum(s1!=s2,na.rm=T)

# l_12_H0 = integrate(Vectorize(model_II,vectorize.args = "rho"),lower=0,upper=1,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix,n_II=n_12,m_II=m_12)[["value"]]
# l_12_H1 = integrate(Vectorize(model_IJ,vectorize.args = "rho"),lower=0,upper=1,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix,n_IJ=n_12,m_IJ=m_12)[["value"]] # still use n_12 and m_12 as the data have not changed
# BF_12 = l_12_H0/l_12_H1

# l_12_H0 = sum(diff(rho_levels)*sapply(rho_mid_pts,model_II,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix,n_II=n_12,m_II=m_12))
# l_12_H1 = sum(diff(rho_levels)*sapply(rho_mid_pts,model_IJ,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix,n_IJ=n_12,m_IJ=m_12))
# BF_12 = l_12_H0/l_12_H1

l_12_H0 = sum(diff(rho_levels)*model_II_vectorized(rho_mid_pts,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix[,locus_position],n_II=n_12,m_II=m_12))
l_12_H1 = sum(diff(rho_levels)*model_IJ_vectorized(rho_mid_pts,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix[,locus_position],n_IJ=n_12,m_IJ=m_12))
BF_12 = l_12_H0/l_12_H1

# l_12_H0 =(1/n_rhos)*sum(model_II_vectorized(rtrunc(n_rhos,"beta",a=0,b=max_rho,shape1=shape_1,shape2=shape_2),genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix[,locus_position],n_II=n_12,m_II=m_12))
# l_12_H1 = (1/n_rhos)*sum(model_IJ_vectorized(rtrunc(n_rhos,"beta",a=0,b=max_rho,shape1=shape_1,shape2=shape_2),genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix[,locus_position],n_IJ=n_12,m_IJ=m_12))
# BF_12 = l_12_H0/l_12_H1
# 


posterior_H0 = 1/(1+(1/BF_12)*(prior_H1/prior_H0))
posterior_H1 = 1 - posterior_H0

m_12
n_12
posterior_H0

##

layout(matrix(1:2))
# tmp = sort(rtrunc(n_rhos,"beta",a=0,b=max_rho,shape1=shape_1,shape2=shape_2))
# plot(tmp,model_II_vectorized(tmp,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix[,locus_position],n_II=n_12,m_II=m_12))
# plot(tmp,model_IJ_vectorized(tmp,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix[,locus_position],n_IJ=n_12,m_IJ=m_12))
# # 
plot(rho_mid_pts,model_II_vectorized(rho_mid_pts,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix[,locus_position],n_II=n_12,m_II=m_12))
plot(rho_mid_pts,model_IJ_vectorized(rho_mid_pts,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix[,locus_position],n_IJ=n_12,m_IJ=m_12))
# 

###

I = rep(1:ncol(observed_genotypes),each=ncol(observed_genotypes))
J = rep(1:ncol(observed_genotypes),ncol(observed_genotypes))


# I = rep(1,each=ncol(observed_genotypes))
# J = rep(1:ncol(observed_genotypes),1)



n_loci_matrix = matrix(mapply(count_n_loci, I, J, MoreArgs=list(data_genotypes_N_combined=t(observed_genotypes))), nrow=ncol(observed_genotypes),byrow=T) # observed mismatch 
observed_mismatch_matrix = matrix(mapply(count_mismatch, I, J, MoreArgs=list(data_genotypes_N_combined=t(observed_genotypes))), nrow=ncol(observed_genotypes),byrow=T) # observed mismatch 

ind_II_matrix = matrix(mapply(
  function(i,j,subject_ID) {
    if(i!=j & i<j){ # only needs to compare i!=j & i<j 
      same_individual = ifelse(subject_ID[i]==subject_ID[j],"Same Individual","Different Individual")
    }
    else{
      same_individual = NA
    }
  }, I, J, MoreArgs=list(subject_ID=subject_ID)), nrow=ncol(observed_genotypes),byrow=T) # obserced mismatch in perct

# Note: This step takes a long time
BF_matrix = matrix( mcmapply(generate_BF,I,J, MoreArgs=list(data_genotypes_N_combined=t(observed_genotypes),data_reads_N_combined=NA,genotypes_pool=1:10,probs_pool_matrix=probs_pool_matrix,reads_cut_off=NA, rho_levels, rho_mid_pts,n_rhos,shape_1,shape_2, max_rho),mc.cores=mc_cores), nrow=ncol(observed_genotypes),byrow=T)

# n_sample = ncol(observed_genotypes)
# n_pairs = (factorial(n_sample)/(factorial(2)*factorial(n_sample-2)))
# prior_H0 = 1/n_pairs # when no paitent ID, assume 
# prior_H1 = 1 - prior_H0

# prior_H0 = 0.5 # when no paitent ID, assume
# prior_H1 = 1 - prior_H0

# prior_H0 = sum(choose(n_size,2)) /choose(sum(n_size),2)
# prior_H1 = 1 - prior_H0

posterior_H0_matrix = 1/(1+(1/BF_matrix)*(prior_H1/prior_H0))
posterior_H1_matrix = 1 - posterior_H0_matrix

# posterior_H0_matrix = matrix(posterior_H0_matrix,byrow=T,nrow=ncol(observed_genotypes))
rownames(posterior_H0_matrix) = colnames(posterior_H0_matrix) = sample_ID

posterior_H0_matrix

data_posterior = data.frame(observed_mismatch_perct= unlist(as.numeric(observed_mismatch_matrix))/unlist(as.numeric(n_loci_matrix)),posterior_H0=unlist(as.numeric(posterior_H0_matrix)),posterior_H1=unlist(as.numeric(1-posterior_H0_matrix)),posterior_odds=unlist(as.numeric(posterior_H0_matrix)/as.numeric(1-posterior_H0_matrix)), ln_BF=2*log(unlist(as.numeric(BF_matrix))),ind_II=unlist(as.character(ind_II_matrix)),I=as.numeric(matrix(I,ncol=ncol(observed_genotypes),byrow=TRUE)),J=as.numeric(matrix(J,ncol=ncol(observed_genotypes),byrow=TRUE))) %>% drop_na(ln_BF)

data_posterior = data_posterior[order(data_posterior$ind_II),]


### some figures ####

# layout(matrix(1:1))
# # ggplot(data = melt(posterior_H0_matrix), aes(x=Var1, y=Var2, fill=value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90,hjust = 1))
# plot(observed_mismatch_matrix,posterior_H0_matrix,xlab="Observed Proportion of mismatch", ylab="Probability belongs to a same individual")
# 


p1 = ggplot(data=data_posterior) + geom_jitter(aes(x=observed_mismatch_perct,y=posterior_H0),alpha=0.6) +
  theme_bw() + scale_color_discrete(name="Observed labelling") +   xlab("Observed Proportion of mismatch") + ylab("Probability belongs to \n a same individual") + theme(legend.position="bottom")



p2 = ggplot(data=data_posterior) + geom_jitter(aes(x=observed_mismatch_perct,y=posterior_H0, color=factor(ind_II)),alpha=0.6) +
  facet_grid(~ind_II) + theme_bw() + scale_color_discrete(name="Observed labelling") +   xlab("Observed Proportion of mismatch") + ylab("Probability belongs to \n a same individual") + theme(legend.position="bottom", strip.text = element_blank())

# p3 = ggplot(data=data_posterior) + geom_histogram(aes(x=posterior_H0,y=..density.., fill=factor(ind_II)),color=gray(0.9),binwidth = 0.01) +  xlab("Posterior probability of H0: two samples are from a same patient") + ylab("Density")+ scale_fill_discrete(name = "Observed labelling") + theme_bw() + theme(legend.position="bottom")

p3 = ggplot(data=data_posterior) + geom_histogram(aes(x=posterior_H0, y=..density..,fill=factor(ind_II)),color=gray(0.9),binwidth = 0.01) +  xlab("Posterior probability of H0: two samples are from a same patient") + ylab("Density")+ scale_fill_discrete(name = "Observed labelling") + theme_bw() + theme(legend.position="bottom", strip.text = element_blank()) + facet_wrap(~ind_II, scales = "free_y")


pushViewport(viewport(layout = grid.layout(3, 1)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1, vp = vplayout(1, 1))
print(p2, vp = vplayout(2, 1))
print(p3, vp = vplayout(3, 1))
grid.text("(a)",x=unit(0.01,"npc"),y=unit(0.98,"npc"))
grid.text("(b)",x=unit(0.01,"npc"),y=unit(0.7,"npc"))
grid.text("(c)",x=unit(0.01,"npc"),y=unit(0.35,"npc"))



## some tables ###

data_full_compare = data_posterior[,c("I","J","posterior_H0")]
# colnames(data_full_compare) = c("sample_1","sample_2","prob_from_same_subject")
data_full_compare$sample_ID_1 = sample_ID[data_full_compare$I]
data_full_compare$sample_ID_2 = sample_ID[data_full_compare$J]

data_full_compare$subject_ID_1 = subject_ID[data_full_compare$I]
data_full_compare$subject_ID_2 = subject_ID[data_full_compare$J]

data_full_compare =  data_full_compare[,c("sample_ID_1","sample_ID_2","subject_ID_1","subject_ID_2","posterior_H0")]
data_full_compare = data_full_compare[order(data_full_compare$sample_ID_2,data_full_compare$posterior_H0, decreasing=T),]
colnames(data_full_compare) = c("sample_ID_1","sample_ID_2","subject_ID_1","subject_ID_2","prob_from_same_subject")

write.csv(data_full_compare,file=paste(results_dir,"/data_full_compare.csv",sep=""),row.names = FALSE)

##
unique_J = sort(unique(J))
best_matched_J = {}
matching_prob = {}
observed_mismatch_perct = {}

for (j in unique_J) {
  
  tmp_1 = subset(data_posterior,data_posterior$I==j | data_posterior$J==j)
  
  max_posterior_H0 = max(tmp_1$posterior_H0)
  tmp_2 = subset(tmp_1,tmp_1$posterior_H0==max_posterior_H0)
  
  tmp_3 = c(tmp_2$I,tmp_2$J)
  tmp_4 = tmp_3[tmp_3!=j] 
  
  if (length(tmp_4)>1) cat("ooops, more than 1 best matched sample")
  
  best_matched_J = c(best_matched_J, tmp_4[1])
  observed_mismatch_perct = c(observed_mismatch_perct, tmp_2$observed_mismatch_perct[1])
  matching_prob = c(matching_prob, max_posterior_H0)
}


data_best_matched = data.frame(sample=sample_ID[unique_J], best_matched=sample_ID[best_matched_J],subject_ID_for_sample = subject_ID[unique_J],subject_ID_for_best_matched=subject_ID[best_matched_J], prob_from_same_subject= matching_prob)

data_mismatched = subset(data_best_matched,as.character(data_best_matched$subject_ID_for_sample)!=as.character(data_best_matched$subject_ID_for_best_matched))
data_mismatched = data_mismatched[order(data_mismatched$sample),]

write.csv(data_best_matched,file=paste(results_dir,"/data_best_matched.csv",sep=""),row.names = FALSE)
write.csv(data_mismatched,file=paste(results_dir,"/data_mismatched.csv",sep=""),row.names = FALSE)

##

end_time <- Sys.time()
end_time - start_time


