#####################"
### Statistics for Genomic Data Science
### Mooc John Hopkin
### Week 1
### Module 1
### Exploratoy analysis
### http://jtleek.com/genstats/inst/doc/01_10_exploratory-analysis.html

install.packages(c("devtools","gplots"))
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("Biobase","org.Hs.eg.db","AnnotationDbi"))
biocLite("alyssafrazee/RSkittleBrewer")


library(gplots)
library(devtools)
library(Biobase)
library(org.Hs.eg.db)
library(AnnotationDbi)

################################################################"
#'RSkittleBrewer' is not available (for R version 3.3.3)#######
#############################################################""

# Colors more info on : http://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html
library(RSkittleBrewer)
# Make the colors pretty
# Load the library and 
#set the color palette with the palette function. Now when I type col = 1 it will look for the first color in the trop colors. We also set the character to be a filled dot with  par(pch=19).
trop = RSkittleBrewer("tropical")
palette(trop)
par(pch=19)             ## POur modifier les ronds sur les graphs : ronds pleins plutôt de vides


con = url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)  ##############  Je charge un fichier d'internet
bm = bodymap.eset
pdata=pData(bm)   #extraction de phenotype table
edata=exprs(bm)   # extract out the expression data  genomic
fdata = fData(bm) # extract out the feature data
ls()              # Voir ce que R a en mémoire

#pdata
names(pdata)
summary(pdata)
table(pdata$gender)   
table(pdata$age)
table(pdata$age,useNA="ifany")  # Attention il peut y avoir des NA
sum(pdata$age=="", na.rm=TRUE)              # POur savoir si il y a des vides en ne tenant pas compte des NA

table(pdata$gender, pdata$race)   # Tableau croisé

is.na(edata)[1,]        # Pour savoir si il y a des NA
sum(is.na(edata))
genes_na = rowSums(is.na(edata))  # si il y a des NA pour les isolers par row (par gene)
table(genes_na)

sample_na = colSums(is.na(edata))  
table(sample_na)



#######################
# Verifier si les dimensions matche
##################
dim(fdata)   # dimension of feature data
# The number of rows of the feature data should match the number of rows of the expression data, since feature data describes the genes, and phenotype data describes the samples.
dim(pdata)   # dimension of the phenotyp data
dim(edata)   # dimension of the expression data

###############################
####### Exploratry analysis part II
####### Plot

boxplot(edata[,1]) # first column
boxplot(log2(edata[,1]+1))
boxplot(log2(edata+1),col=2,range=0)  #one bow plot per column of the date set
help(log2)

par(mfrow=c(1,2))  #  The other thing that you can do is you can make histograms. So I'm actually going to show two histograms. And if I want to show plot side by side, I can use the par(mfrow=c(1,2)). This is basically saying set up the screen so that it has one row and two columns of plots. So it's side by side plots. So I can make a histogram of the values from the first sample. So again I'm going to do this transform here because I think it makes it a little easier to see. And then I'm going to set the color equal to two, so this is going to make this nice blue histogram. And so you can see here right away that almost all the values are equal to zero. And then you see some values out to the right there. And so, then, I can make that same plot. But for the second sample, and so, when I do that, 
hist(log2(edata[,1]+1),col=2)
hist(log2(edata[,2]+1),col=2)

par(mfrow=c(1,1))  # one by one plot
plot(density(log2(edata[,1]+1)),col=2)
plot(density(log2(edata[,2]+1)),col=2)

lines(density(log2(edata[,2]+1)),col=3) # overlay superposer les graphs
qqplot(log2(edata[,1]+1),log2(edata[,2]+1),col=3)
abline(c(0,1))

## QQ PLOT 
#I can see a bunch of data points here and so one thing to keep in mind is that when you make it this plot now, this q-q plot is making one dot for every gene. And so, or, sorry, for every quantile of these two distributions. And so it depends on the number of quantiles that it has to calculate how long that'll take. 
#And so for example this dot right here this says this is say the 5th percentile so this is the 5th percentile for the second sample is on the Y axis. And the fifth percentile for first samples on the x axis. So we can see that it's above the 45 degree line here and so you can see that the second sample has a higher fifth percentile than the first sample so it has higher values for low values. If you want to be able to see that a little more clearly you can add a 45 degree line. And one way you can do that is with this abline command, abline. And you can tell it, you use this intercept and this slope and so you can see here this has an intercept of zero and a slope of one. That's the 45 degree line. And so, if I add that on top so there's the 45 degree line, you can see here the quantiles are a little bit larger for the second sample down low and they're a little bit lower up high. And so you can kind of see how the two distributions compare to each other. The other thing you can do when comparing samples, is you can make what's called an MA or a Bland-Altman plot. So here I'm going to take the difference between the two samples and I'm going to make that the Y axis so I'm going to take the difference between sample one and sample two and then I'm going to add the two samples up. 

mm = log2(edata[,1]+1)-log2(edata[,2]+1)
aa = log2(edata[,1]+1)+log2(edata[,2]+1)
plot(aa,mm,col=2)
# MA or a Bland-Altman plot. 
# Pour notamment des réplicats techniques
#So on the X axis we have the sum of the two samples. So basically moving from left to right you get lower expression or lower counts to higher counts. And on that Y axis, you're to take the difference. So if it's at zero, that means there's no difference between the two samples. So you can see for example trends here. You can see as you get higher and higher accounts there's sort of a trend that appears to be that the samples get closer and closer together so each dot here is one gene and so I'm just taking the difference between the two samples for that gene. So this MA kind of plot people make very often you want to see that it's, especially for technical replicates, there's some kind of replicates that should be similar. You want to see it centered on the zero line, and you'd like to see it with low variability and no trends that are dependent on the total number of counts for those samples. And so, that's a way that you can sort of make those plots to compare them. Now, the next thing that 



# Dataset -> data frame
# pour faire des filtres pour enlever 
edata = as.data.frame(edata) 
#edata = as.matrix(edata) pour revenir à une matrice
filt_edata = filter(edata,rowMeans(edata) > 1)
dim(filt_edata)
# Probleme la fonction filter n'a rien enlever alors qu'elle aurait dû
boxplot(as.matrix(log2(filt_edata+1)),col=2)


########################################
############ Exploratory analysis part 3
############ CHECK OF CONSISTENCY