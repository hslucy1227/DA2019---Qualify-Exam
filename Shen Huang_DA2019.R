rm(list=ls())
library(nlme)
library(lattice) 
library(lme4)    
library(pastecs)
library(MuMIn)
library(ggplot2)
library(pROC)
library(devtools)               
library(ggbiplot)
options(scipen = 200)

# load data
data = read.table("D:/UCI/19spring/DA2019/DA2019.csv", sep=",", header = T)
dim(data)
head(data)
str(data)
# transfer gender as binary variable
data$Gender = as.numeric(data$Gender) - 1
summary(data[, c(1:5,dim(data)[2])])
stat.desc(data[, c(1:5,dim(data)[2])])


##### EDA
# EDA-age
ggplot(data[,c(1,3)], aes(x=data$Age)) +
  geom_histogram(color = "darkblue", fill = "lightblue", binwidth = 1) +
  labs(x = "Age", title = "                Distribution of Age") +
  theme(plot.title = element_text(colour = "black", face = "bold", size = 19, vjust = 1))
table(data$Age)

# EDA-RACE
ggplot(data[,c(1,5)], aes(x=data$RACE)) +
  geom_histogram(color = "darkblue", fill = "lightblue",binwidth=1) +
  labs(x = "RACE score", title = "       Distribution of RACE score") +
  theme(plot.title = element_text(colour = "black", face = "bold", size = 19, vjust = 1))
table(data$RACE)

# EDA-gender
table(data$Gender)

# EDA-age/gender
means.age = tapply(data$Stroke, list(data$Age, data$Gender), mean)
means.age[is.na(means.age)] = 0
barplot(t(means.age), beside=TRUE, col=c("Red","Blue"),
        ylim = c(0,1), legend = c("Female","Male"),
        xlab = "Age",ylab = "Probability of stroke",
        main = "Probability of Stroke by Age and Gender",cex.main = 1.4,las = 1,
        args.legend = list(x = "topleft", bty = "n"))

# EDA-LKW
ggplot(data[,c(1,4,106)], aes(x = data$LKW)) +
  geom_histogram(color="darkblue", fill="lightblue",binwidth=1) +
  labs(x = "Last Known Well time(hour)", title = "        Distribution of Last Known Well time") +
  theme(plot.title = element_text(colour = "black", face = "bold", size = 15, vjust = 1))
table(data$LKW)


##### Missing Data Imputation
# find number of missing value for each column
missing.num = vector()
for(i in 1:105){
  missing.num = append(missing.num, sum(is.na(data[,i])))
}
missing.num

# impute missing data by mean of each column
for(i in 1:105){
  data[,i][is.na(data[,i]) == 1] = mean(na.omit(data[,i]))
}
summary(data)


##### Correlation of EEG
corr = cor(na.omit(data[,6:105]))
summary(corr)
summary(unique(corr[abs(corr) >= 0.5]))


##### Modelling
set.seed(3)  
samp = sample(nrow(data), nrow(data) * 0.9)
da.tr = data[samp,]
da.va = data[-samp,]

# association between RACE & Stroke
mod1 = glm(Stroke ~ 1 + Gender + Age + LKW + RACE, family = binomial(link = "logit"), data = da.tr)
summary(mod1)

mod12 = glm(Stroke ~ 1 + Gender + Age + LKW, family = binomial(link = "logit"),data = da.tr)
summary(mod12)
anova(mod12, mod1)
pchisq(9.4438, 1, lower.tail=F)

# goodness of fit
pchisq(mod1$deviance, mod1$df.residual, lower.tail = F) 

# compared to the null model
pchisq(121.14, 89, lower.tail = F)

# prediction
mod1.pred = predict(mod1, da.va)
mod1p = exp(mod1.pred)/(exp(mod1.pred) + 1)
auc = roc(da.va$Stroke, mod1p, plot = TRUE) 
legend('bottomright', legend = paste('AUC =', round(auc$auc, 3)), bty = 'n')

ggroc(auc, color = "blue") + ggtitle("ROC curve of Model 1") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "darkgrey", linetype = "dashed") 

# add interaction term
mod2 = glm(Stroke ~ 1 + Gender + Age + LKW + RACE + RACE*Age, family = binomial(), data = da.tr)
summary(mod2)
anova(mod1, mod2)
pchisq(2.0167, 1, lower.tail = F)

# goodness of fit
pchisq(mod2$deviance, mod2$df.residual, lower.tail = F) 

# compared to the null model
pchisq(121.14,89,lower.tail = F)
mod2.pred = predict(mod2, da.va)
mod2p = exp(mod2.pred)/(exp(mod2.pred) + 1)
auc = roc(da.va$Stroke, mod2p, plot = TRUE) 
legend('bottomright', legend = paste('AUC =', round(auc$auc, 3)), bty = 'n')


##### PCA
# conduct PCA on training dataset
pca = prcomp(da.tr[,6:105], retx = TRUE, center = TRUE, scale = TRUE)
ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = factor(da.tr$Stroke),ellipse = TRUE, circle = TRUE) + 
  theme(legend.direction = 'horizontal', legend.position = 'top')
screeplot(pca,type = "line",main = "Scree plot of PCA", npcs = 30)
# percent explained variance
expl.var = pca$sdev^2/sum(pca$sdev^2) * 100 
#summary(pca)

num = c(1:50)
sum.var = cumsum(expl.var[1:50])
sumvar = data.frame(num, sum.var)
ggplot(data = sumvar, mapping = aes(x = num, y = sum.var), color = 'blue')+
  geom_smooth(method="loess") +
  theme(axis.text.x = element_text(size = 11,color = "black"), 
        axis.text.y = element_text(size = 11,color = "black")) + 
  labs(x = "Number of Components", y = "Cumulative Explained Variance (%)") +
  theme(axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12) )

# training component
pca.tr = predict(pca, newdata = da.tr[,6:105])
comp.tr = data.frame(comp1 = rep(0,90),comp2 = rep(0,90),comp3 = rep(0,90),comp4 = rep(0,90),
                     comp5 = rep(0,90),comp6 = rep(0,90),comp7 = rep(0,90),comp8 = rep(0,90),
                     comp9 = rep(0,90),comp10 = rep(0,90),comp11 = rep(0,90),comp12 = rep(0,90),
                     comp13 = rep(0,90),comp14 = rep(0,90),comp15 = rep(0,90),comp16 = rep(0,90),
                     comp17 = rep(0,90),comp18 = rep(0,90),comp19 = rep(0,90),comp20 = rep(0,90))
for(i in 1:20){
  comp.tr[,i] = pca.tr[, i]
}
da.tr = cbind(da.tr, comp.tr)

# predict
pca.pred = predict(pca, newdata = da.va[, 6:105])
comp = data.frame(comp1 = rep(0,10),comp2 = rep(0,10),comp3 = rep(0,10),comp4 = rep(0,10),
                  comp5 = rep(0,10),comp6 = rep(0,10),comp7 = rep(0,10),comp8 = rep(0,10),
                  comp9 = rep(0,10),comp10 = rep(0,10),comp11 = rep(0,10),comp12 = rep(0,10),
                  comp13 = rep(0,10),comp14 = rep(0,10),comp15 = rep(0,10),comp16 = rep(0,10),
                  comp17 = rep(0,10),comp18 = rep(0,10),comp19 = rep(0,10),comp20 = rep(0,10))
for(i in 1:20){
  comp[, i] = pca.pred[, i]
}
da.va = cbind(da.va, comp)

# glm with EEG
mod3 = glm(Stroke ~ 1 + Gender + Age + LKW + RACE + comp1 + comp2+ comp3 + comp4 + comp5 +
             comp6 + comp7 + comp8+ comp9 + comp10 + comp11 + comp12 + comp13 + comp14 + 
             comp15 + comp16 + comp17 + comp18 + comp19 + comp20, family = binomial(), data = da.tr)
summary(mod3)

# goodness of fit
pchisq(mod3$deviance, mod3$df.residual, lower.tail = F) 

# conpare to the null model
pchisq(121.14,89, lower.tail = F)
mod3.pred = predict(mod3, da.va)
mod3p = exp(mod3.pred)/(exp(mod3.pred) + 1)
auc2 = roc(da.va$Stroke, mod3p, plot = TRUE) 
legend('bottomright', legend = paste('AUC =', round(auc$auc, 3)), bty = 'n')

ggroc(auc2, color = "blue") +  ggtitle("ROC curve of Model 3") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "darkgrey", linetype = "dashed") 

