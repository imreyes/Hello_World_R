# Confounding, Batch Effects.
f3# Sources are from the below link on 10/30/2016:
# https://courses.edx.org/courses/course-v1:HarvardX+PH525.4x+2T2016/courseware/

# Firstly let's look at a demo - Small dataset of school admission, categorized with major & gender.
setwd('G:/R/Project1/')
load('admissions (1).rda')
print(admissions)   # Explore the data.
index<-which(admissions$Gender==1)    # Pick men's data.
accepted<-sum(admissions$Number[index]*admissions$Percent[index]/100)
applied<-sum(admissions$Number[index])
accepted/applied          # Rate of men acceptance.

# Look at women now.
index<-which(admissions$Gender==0)    # Pick women's data.
accepted<-sum(admissions$Number[index]*admissions$Percent[index]/100)
applied<-sum(admissions$Number[index])
accepted/applied          # Rate of women acceptance.

idx<-admissions$Gender==1
acp<-round(
  matrix(
    c(
      # Male accepted:
      sum(admissions$Number[idx]*admissions$Percent[idx]),
      # Male rejected:
      sum(admissions$Number[!idx]*admissions$Percent[!idx]),
      # Female accepted:
      sum(admissions$Number[idx]*(100-admissions$Percent[idx])),
      # Female rejected:
      sum(admissions$Number[!idx]*(100-admissions$Percent[!idx]))
      ),
    2,2)/100,
  0)
chisq.test(acp)$p.val
# Note the p-value indicates the significance of gender difference in admission,
# which was proven to be confounding.

# In order to reveal the reason behind the issue, we stratify the entire data by major.
# Here we define 'easy' and 'hard' majors by the acceptance rate of each major.
H<-sapply(levels(admissions$Major),function(i){
  idx<-admissions$Major==i
  accepted<-sum(admissions$Number[idx]*admissions$Percent[idx]/100)
  accepted/sum(admissions$Number[idx])
})

# Look at the correlation between # of application and H.
cor(admissions$Number[admissions$Gender==1],H)  # For men.
cor(admissions$Number[admissions$Gender==0],H)  # For women.


#==================================================================================================#
# Now it's real gene data.
library(Biobase)
load('GSE5859.rda')
# We can extract the gene expression data and sample info table
# using the Bioconductor fucntion exprs() and pData():
geneExpression<-exprs(e)
sampleInfo<-pData(e)

# See how many and which years of data contain more than one ethnic group.
year<-format(sampleInfo$date,'%y')
unique(year)
eths<-sapply(unique(year),function(i){
  length(unique(sampleInfo$ethnicity[year==i]))
})
sum(eths>=2)

# Now consider months as well.
month.year<-format(sampleInfo$date,'%m%y')
eths2<-sapply(unique(month.year),function(i){
  length(unique(sampleInfo$ethnicity[month.year==i]))
})
mean(eths2>=2)

# Still look into years: subset the CEU ethnic group out, and look at
# differences between data obtained in 2002 vs 2003.
library(genefilter)
idx<-which(sampleInfo$ethnicity=='CEU'&year%in%c('02','03'))
fac<-factor(format(sampleInfo$date[idx],'%y'))
pvals<-rowttests(geneExpression[,idx],fac)$p.val
library(qvalue)
qvals<-qvalue(pvals)$qval
sum(qvals<0.05)
qvalue(pvals)$pi0
# Note there are ~50% data called significant, which is unexpected.

# Now compare 2003 vs 2004:
idx<-which(sampleInfo$ethnicity=='CEU'&year%in%c('03','04'))
fac<-droplevels(factor(year[idx]))          # A more efficient way to get 'fac' factor; see above block for comparison.
pvals<-rowttests(geneExpression[,idx],fac)$p.val
qvals<-qvalue(pvals)$qval
sum(qvals<0.05)

# One additional: how does it look like within a year?
idx<-which(sampleInfo$ethnicity=='CEU'&year=='02')
fac<-factor(sample(c(0,1),length(idx),replace=T))
pvals<-rowttests(geneExpression[,idx],fac)$p.val
qvals<-qvalue(pvals)$qval
sum(qvals<0.05)
# There's exactly 0 genes called significant - different years do give biased results.

# This time, let's look at ethnicity.
# We'll compare between ASN and CEU.
idx<-which(sampleInfo$ethnicity%in%c('CEU','ASN'))
fac<-droplevels(factor(sampleInfo$ethnicity[idx]))    # Here droplevels() is redundant; see below block.
pvals<-rowttests(geneExpression[,idx],fac)$p.val
qvals<-qvalue(pvals)$qval
sum(qvals<0.05)
# Note over 80% of genes expressed differently when comparing between ethnic groups.

# Now we need to find a solution, by stratification.
# For illustration and test purposes, we compare between CEU vs ASN,
# but only the data recorded in the year of 2005.
idx<-which(sampleInfo$ethnicity%in%c('CEU','ASN')&year=='05')
fac<-factor(sampleInfo$ethnicity[idx])
pvals<-rowttests(geneExpression[,idx],fac)$p.val
qvals<-qvalue(pvals)$qval
sum(qvals<0.05)

# Finally, let's compare between CEU at 2002 vs ASN at 2005. Pick 3 random samples for CEU group.
set.seed(3)
idx1<-sample(which(sampleInfo$ethnicity=='CEU'&year=='02'),3)
idx2<-which(sampleInfo$ethnicity=='ASN'&year=='05')
idx<-c(idx1,idx2)
fac<-factor(sampleInfo$ethnicity[idx])
pvals<-rowttests(geneExpression[,idx],fac)$p.val
qvals<-qvalue(pvals)$qval
sum(qvals<0.05)



#==============================================================================================================#
# Now that have we seen the influence of batch effect, we hope to find a solution to parse it out of our models.
# Let's work with the subset again.
setwd('G:/R/Project1/')
load('GSE5859Subset.rda')
sex<-sampleInfo$group
month<-factor(format(sampleInfo$date,'%m'))
table(sex,month)

# Calculate the sex-based difference in genes.
library(qvalue)
pvals<-rowttests(geneExpression,factor(sex))$p.val
qvals<-qvalue(pvals)$qval
sum(qvals<0.1)

# It's believed that gender-correlated differences are mostly on chrY and chrX:
mean(geneAnnotation$CHR[qvals<0.1]%in%c('chrX','chrY'))

# And now look at the date influences
# on the autosomes called significant by sex(not chrY or chrX):
chridx<-which(qvals<0.1&!geneAnnotation$CHR%in%c('chrX','chrY'))
# Note:which(!geneAnnotation$CHR[qvals<0.1]%in%c(chr'X','chrY'))
# reads indeces of new vector of geneAnnotation$CHR[qvals<0.1], not geneAnnotation$CHR.
pvals<-rowttests(geneExpression[chridx,],factor(month))$p.val
# Note here qvalue() gives error.
mean(pvals<0.05)

# Now we try to make a regression model
# taking both sex and month into account.
X<-model.matrix(~sex+month)   # matrix of factor sex & month.
pvals<-sapply(1:dim((geneExpression))[1],function(i){
  y<-geneExpression[i,]
  fit<-lm(y~X-1)
  summary(fit)$coef[,4]
})
pvals <- t(pvals)
qvals.Sex <- qvalue(pvals[, 2])$qval                # Grouped by sex.
sum(qvals.Sex<0.1)
mean(geneAnnotation$CHR[qvals.Sex < 0.1] %in% c('chrX', 'chrY'))

qvals.Month <- qvalue(pvals[, 3])$qval              # Now by month - much bigger!.
sum(qvals.Month < 0.1)
mean(geneAnnotation$CHR[qvals.Month < 0.1] %in% c('chrX', 'chrY'))      # And independent of sex.
