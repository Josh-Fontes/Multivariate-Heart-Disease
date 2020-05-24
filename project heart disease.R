



## analysis ideas:

# responses - resting blood pressure, cholestoral


install.packages("dplyr")
library(dplyr)

install.packages("MASS")
library(MASS)





setwd("C:\\Users\\Julio\\Documents\\JOSH\\SJSU\\math 257") 


heartdisease <- read.table("heart.csv", header = TRUE, sep = ",")
View(heartdisease)


objects(heartdisease)



class(heartdisease$sex)
heartdisease$sex <- as.factor(heartdisease$sex)

class(heartdisease$ca)
heartdisease$ca <- as.factor(heartdisease$ca)


# grouping ca4 with ca3 because zero females are in the ca4 group
levels(heartdisease$ca) <- c("0", "1", "2", "2", "2")
heartdisease$ca

# number observations per cell (for each combination of levels) - Unbalanced
k <- matrix(0,ncol=2,nrow=3)
for (b in 0:1){
  for (a in 0:2){
    k[a+1,b+1] <- nrow(heartdisease%>%filter(ca==a & sex==b))
  }
}


k

nl <- as.vector(matrix(k, ncol = 6, nrow = 1))
nl



## creating the TWO-WAY MANOVA model

# Specify the linear relationship in MANOVA
fit.lm.complete <- lm(cbind(chol, trestbps, thalach) ~  sex + ca + sex*ca, data = heartdisease) 
summary(fit.lm.complete)
manova(fit.lm.complete)
fit.manova.complete <- manova(fit.lm.complete)
summary(fit.manova.complete, test = "Wilks") # interaction not significant



# without interaction now
fit.lm <- lm(cbind(chol, trestbps, thalach) ~  sex + ca, data = heartdisease) 
summary(fit.lm)



# Run the Manova
manova(fit.lm)
fit.manova <- manova(fit.lm)
summary(fit.manova) # default test statistic used is Pillai


# See results for each of the tests
summary(fit.manova, test = "Pillai")
summary(fit.manova, test = "Wilks")
summary(fit.manova, test = "Hotelling-Lawley")
summary(fit.manova, test = "Roy")


# univariate examinations
summary.aov(fit.lm)



# summary examinations of manova model
summary(fit.manova)$SS  # sum of squares
summary(fit.manova)$Eigenvalues







############## ###################### ##################################
# just checking to see that we get the Wilks' Lambda statistic from our SS matrices


summary(fit.manova, test = "Wilks")
summary(fit.manova)$SS


SSf1 <- as.matrix(as.data.frame(summary(fit.manova)$SS[1]))
SSres <- as.matrix(as.data.frame(summary(fit.manova)$SS[3]))


det(SSres)/(det(SSres + SSf1)) # and we do

























##### ASSUMPTIONS check:

#### multivariate normality


## we need to check residuals of our full model correct? - yes - also do mahalanobias distances
par(mfrow = c(1,3))
# normal quantile plots
students.residuals <- rstudent(fit.manova)
qqnorm(students.residuals[,1], main = "Cholesterol")
qqline(students.residuals[,1])

qqnorm(students.residuals[,2], main = "RBS")
qqline(students.residuals[,2])

qqnorm(students.residuals[,3], main = "Maximum Heart Rate")
qqline(students.residuals[,3])





# compute the statistical distiances


# First convert the data frame to a matrix

hd.ex <- heartdisease[,c(5,4,8,2,12)]
head(hd.ex)

hd.M <- data.matrix(hd.ex[,1:3])

head(hd.M)




(Si <- solve(var(hd.M)))

# Two ways find the sample mean vector of multivariate data
(xbar <- colMeans(hd.ex[1:3]))




# Finding all the statistical distances 

# First, get dimensions of our data
n <- nrow(hd.M)
p <- ncol(hd.M)


(t(hd.M[1,] - xbar)%*% Si %*%(hd.M[1,] - xbar))




# Then, calculate all Mahalanobis distances
d <- sapply(1:n, function(k) (t(hd.M[k,] - xbar)%*%Si%*%(hd.M[k,] - xbar)))



# Compute quantiles of a chi-square distribution
q1 <- qchisq((1:n-0.5)/n,p)



par(mfrow = c(1,1))
qqplot(q1,d,xlab = "Chi-square quantiles", ylab = "Sample statistical distances", main = "Chi-square Probability Plot")
lines(q1,q1)


# maybe without the outlier?
max(d)
d2 <- d[-86]


qqplot(q1,d2,xlab = "Chi-square quantiles", ylab = "Sample statistical distances", main = "Chi-square Probability Plot")
lines(q1,q1)


hd.ex[86,]

# one female just has suuuuper high cholesterol
max(hd.ex[,1])
boxplot(hd.ex[,1])

# not really any really high outliers for the other responses
boxplot(hd.ex[,2])
boxplot(hd.ex[,3])












#### need to compare the 10 3x3 covariance matrices of our responses paired with each combination of levels from our factors - S1, S2, ..., S8


head(hd.ex)


class(hd.ex)

# 3x3 covariance matrices for each combination of levels from factor A and factor B - need to check that all of these are equal to each other
S1 <- var(hd.ex[hd.ex$sex == 0 & hd.ex$ca == 0, 1:3])
S2 <- var(hd.ex[hd.ex$sex == 0 & hd.ex$ca == 1, 1:3])
S3 <- var(hd.ex[hd.ex$sex == 0 & hd.ex$ca == 2, 1:3])
S4 <- var(hd.ex[hd.ex$sex == 0 & hd.ex$ca == 3, 1:3])

S5 <- var(hd.ex[hd.ex$sex == 1 & hd.ex$ca == 0, 1:3])  # looking at sample covariance matrix for the 3 response variables for males with 0 blood vessels 
S6 <- var(hd.ex[hd.ex$sex == 1 & hd.ex$ca == 1, 1:3])
S7 <- var(hd.ex[hd.ex$sex == 1 & hd.ex$ca == 2, 1:3])
S8 <- var(hd.ex[hd.ex$sex == 1 & hd.ex$ca == 3, 1:3])

# maybe treat this test as a one-way MANOVA problem with 8 levels of one factor??
g <- 8
p <- 3



C1 <- 1 - (sum(1/(nl-1)) - 1/sum(nl-1))*((2*p^2 + 3*p -1)/(6*(p + 1)*(g-1)))




Spool <- 1/sum(nl-1)*((nl[1]-1)*S1 + (nl[2]-1)*S2 + (nl[3]-1)*S3 + (nl[4]-1)*S4 + (nl[5]-1)*S5 + (nl[6]-1)*S6 + (nl[7]-1)*S7 + (nl[8]-1)*S8)



M <- (sum(nl-1))*log(det(Spool)) - ((nl[1]-1)*log(det(S1)) + (nl[2]-1)*log(det(S2)) + (nl[3]-1)*log(det(S3)) + (nl[4]-1)*log(det(S4)) + (nl[5]-1)*log(det(S5)) + (nl[6]-1)*log(det(S6)) + (nl[7]-1)*log(det(S7)) + (nl[8]-1)*log(det(S8)))

  
teststat <- M*C1



crit <- qchisq(0.95, (p*(p+1)*(g-1))/2)


teststat > crit


# reject Ho -> conclude that the covariance matrics are not equal




  
  
  
  

# checks out
sum(1/(nl-1))
1/(64-1) + 1/(15-1) + 1/(13-1) + 1/(4-1) + 1/(111-1) + 1/(50-1) + 1/(25-1) + 1/(21-1)



######### Univariate normality for each combination of levels   ###########  do not need to do this

heartdisease.factors <- heartdisease[,c(2,12)]
head(heartdisease.factors)


residuals.M <- cbind(students.residuals, heartdisease.factors)

summary(residuals.M)
head(residuals.M)


par(mfrow=c(2,3))
for (a in 0:1){
  for (b in 0:3){
    h1 <- residuals.M%>%filter(ca==b & sex==a)
    for (i in 1:3){
      qqnorm(h1[,i],main=paste(c('Normal QQ Plot of ', colnames(h1)[i])))
      qqline(h1[,i])
    }
  }
}

h1 <- residuals.M%>%filter(ca==4 & sex==1)
for (i in 1:3){
  qqnorm(h1[,i],main=paste(c('Normal QQ Plot of ', colnames(h1)[i])))
  qqline(h1[,i])
}









######## testing for multivariate normaility
install.packages("mvnormtest")
library(mvnormtest)

attach(heartdisease)
Y <- as.matrix(cbind(chol, trestbps, thalach))
head(Y)

mshapiro.test(t(Y)) # no multivariate normality



Y <- as.matrix(cbind(thalach, trestbps))

mshapiro.test(t(Y))  # much better p-value without cholesterol

detach(heartdisease)




################################################################################












############ Simultaneous confidence intervals  ###############

p <- 3   # number of response variables
g <- 2   # number of levels of factor 1 (sex)
b <- 4   # number of levels of factor 2 (ca)


# number observations per cell (for each combination of levels) - Unbalanced
k <- matrix(0,ncol=2,nrow=4)
for (b in 0:1){
  for (a in 0:3){
    k[a+1,b+1] <- nrow(heartdisease%>%filter(ca==a & sex==b))
  }
}

k



# all different sample sizes CI's (sample size for each combination of levels)
nfemaleca0 <- k[1,1]  # number of females with 0 blood vessels colored...
nfemaleca1 <- k[2,1]  # number of females with 1 blood vessel collored...
nfemaleca2 <- k[3,1]
nfemaleca3 <- k[4,1]
nmaleca0 <- k[1,2]  # number of males with 0 blood vessels colored...
nmaleca1 <- k[2,2]  # number of males with 1 blood vessel colored...
nmaleca2 <- k[3,2]
nmaleca3 <- k[4,2]







hd.ex <- heartdisease[,c(5,4,8,2,12)]
head(hd.ex)

grandmean <- colMeans(hd.ex[,1:3])
meanfemaleca0 <- colMeans(hd.ex[hd.ex$sex == 0 & hd.ex$ca == 0, 1:3])
meanfemaleca1 <- colMeans(hd.ex[hd.ex$sex == 0 & hd.ex$ca == 1, 1:3])
meanfemaleca2 <- colMeans(hd.ex[hd.ex$sex == 0 & hd.ex$ca == 2, 1:3])
meanfemaleca3 <- colMeans(hd.ex[hd.ex$sex == 0 & hd.ex$ca == 3, 1:3])

meanmaleca0 <- colMeans(hd.ex[hd.ex$sex == 1 & hd.ex$ca == 0, 1:3])
meanmaleca1 <- colMeans(hd.ex[hd.ex$sex == 1 & hd.ex$ca == 1, 1:3])
meanmaleca2 <- colMeans(hd.ex[hd.ex$sex == 1 & hd.ex$ca == 2, 1:3])
meanmaleca3 <- colMeans(hd.ex[hd.ex$sex == 1 & hd.ex$ca == 3, 1:3])



## mean differences in response between males and females at different levels of Factor 2

meandiffca0 <- meanfemaleca0 - meanmaleca0  # mean difference in response between males and females, who had zero major blood vessels colored by flouroscopy...
meandiffca1 <- meanfemaleca1 - meanmaleca1
meandiffca2 <- meanfemaleca2 - meanmaleca2
meandiffca3 <- meanfemaleca3 - meanmaleca3




## error terms
summary(fit.manova)$SS[3]
Errors <- as.data.frame(summary(fit.manova)$SS[3])
Errors  # need diagonal elements for confidence intervals
errorChol <- Errors[1,1]
errorRbps <- Errors[2,2]
errorThalach <- Errors[3,3]


# critical value example (for comparing mean difference in response between males and females who had 0 blood vessels colored by flouroscopy)

vca0 <- g*b*(nfemaleca0 + nmaleca0 - 2) # degrees of freedom

tmultca0 <- qt(0.05/(p*g*(g-1)), vca0, lower.tail=FALSE)
tmultca0




## Example CI for difference in mean cholestoral between females and males who had 0 major blood vessels colored by flouroscopy
LB <- meandiffca0[1] - tmultca0*sqrt(errorChol/vca0*(1/(b*nfemaleca0) + 1/(b*nmaleca0)))
UB <- meandiffca0[1] + tmultca0*sqrt(errorChol/vca0*(1/(b*nfemaleca0) + 1/(b*nmaleca0)))
c(LB,UB)













################ PC ANALYSIS ####################################################

head(heartdisease)

Y <- heartdisease[,c(4,5,8,10)]
head(Y)


Y.M <- as.matrix(Y)
head(Y.M)






# i think it makes more sense to use the correlation matrix
(hd.pcr <- prcomp(Y.M, scale = T))


(cumvar <- cumsum(hd.pcr$sdev^2)/sum(hd.pcr$sdev^2))
(propvar <- hd.pcr$sdev^2/sum(hd.pcr$sdev^2))





## checking stuff


## not what i would expect, but the following averages are consistent with PC1
# target = 1 means the patient

# average age
mean(heartdisease[which(heartdisease$target == 1),1])
mean(heartdisease[which(heartdisease$target == 0),1])

# average cholestoral
mean(heartdisease[which(heartdisease$target == 1),5])
mean(heartdisease[which(heartdisease$target == 0),5])

# average resting blood pressure
mean(heartdisease[which(heartdisease$target == 1),4])
mean(heartdisease[which(heartdisease$target == 0),4])

# average maximum heart rate achieved
mean(heartdisease[which(heartdisease$target == 1),8])
mean(heartdisease[which(heartdisease$target == 0),8])




## isnt this the opposite of what we would expect?
## isnt it known that higher cholesterol is at least correlated with heart disease?
# - apparently it is not entirely clear.  After looking into it, they may be associated with each other but causation has not been clear
# so confused






# Plot component scores of the PCs:

par(pch=5,fin=c(5,5), mfrow = c(1,1))
choose<-c(1,2,3,4)

pairs(hd.pcr$x[ ,choose],labels=c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3]), expression(lambda[4])))




# Plot proportion of variance explained by each component & cumulative proportion of variance
plot(propvar, ylim=c(0,1), xaxt = "n", main = "Proportion of variance explained by each PC", xlab = "principal component", ylab = "Proportion of variance explained", pch = 16, bty = "n")
axis(1, at = c(1:4), labels = c(expression(lambda[1]), expression(lambda[2]), expression(lambda[3]), expression(lambda[4])))
lines(propvar, lty = 19)
lines(cumvar, lty = 20)
points(cumvar, pch = 21, bg = "white")
legend(3, 0.5, legend=c("Proportion of variance", "Cumulative variance"), lty = c(19, 20), pch = c(16, 21), pt.bg = c(NA,"white"), cex = 0.8, bty = "n")





###### INTERPRET  ----  PC's are less about direction and more about the magnitude of maximum variability in the data - more about grouping variables together - can multiply PC1 by (-1) if it helps to interpret




# The first 3 principal components capture about 84.4% of variability in the correlation matrix.  It may be possible to adequately summarize the data with only 3 variables, rather than 4.  However, the last principal component captures about 15.6% of the variability in the correlation matrix, which is arguably important.


# The first principal component is an aggregate "heart disease" measure.  Cholesterol, resting blood pressure, and old peak are contrasted against maximum heart rate achieved.  This makes practical sense because higher cholesterol, resting blood pressure, and ST depression are generally associated with an unhealthy human body.  Higher values of maximum heart rate are good because this generally indicates a strong and healthy heart.


# Principal component 2 is a contrast between old peak and all other response variables, although old peak has a rather small coefficient.  PC2 is mainly dominated by the grouping of resting blood pressure, cholesterol, and maximum heart rate.


# Principal component 3 is contrast between cholesterol and all other response variables.  Although, PC3 is dominated by cholesterol and resting blood pressure.


# Principal component 4 is a contrast between resting blood pressure and all other response variables.  Although, PC4 is dominated by maximum heart rate achived and old peak.








## just trying to understand

## comparing scatter plots
plot(heartdisease$chol ~ heartdisease[,1]) # cholesterol by age
plot(heartdisease$trestbps ~ heartdisease[,1])
plot(heartdisease$oldpeak ~ heartdisease[,1])
plot(heartdisease$thalach ~ heartdisease[,1])
plot(heartdisease$chol ~ heartdisease$thalach) # cholesterol by maximum heart rate achieved
plot(heartdisease$trestbps~ heartdisease$chol)




## comparing SLR for each response based on age to see how they relate

# positive correlations
summary(lm(chol ~ heartdisease[,1], data = heartdisease))
summary(lm(trestbps ~ chol, data = heartdisease))
summary(lm(oldpeak ~ heartdisease[,1], data = heartdisease))


# negative correlation
summary(lm(thalach ~ heartdisease[,1], data = heartdisease)) 




