
library(lme4)
library(MCMCpack)
library(lattice)
library(plotrix)
								



scd.data <- read.csv('../data/scd_RAL_data_final.csv', na.strings="NA", stringsAsFactors = TRUE)
scd.data <- na.omit(scd.data)

#Filter out the WT's from the dataset
scd.data <- scd.data[!((1:length(scd.data[,1])) %in% grep("WT", as.character(scd.data$genotype))),]
#And get rid of the factor levels that are no longer used
scd.data$genotype <- droplevels(scd.data$genotype)

#####################################################
#1. Penetrance of ts2 comb differs among lines
#####################################################
penetrance.data <- data.frame(	genotype = scd.data$genotype, 
								ectopic_present = (scd.data$t2_teeth > 0)
							  )
							  
penetrance.data.tabulated <- data.frame(aggregate(penetrance.data[,2],
    by=list(genotype=penetrance.data$genotype), FUN=sum), 
        aggregate(penetrance.data[,2], 
        by=list(genotype=penetrance.data$genotype), FUN=length)[2])
        
penetrance.model <- glm(ectopic_present ~ genotype, 
    family = binomial(link="logit"), 
    data = penetrance.data)

chi2 <- penetrance.model$null.deviance - penetrance.model$deviance

chi2
pchisq(chi2, df=1, lower.tail=F)
# df = number of parameters - 1, right?
pchisq(chi2, df=length(penetrance.model$coefficients)-1, lower.tail=F)

#Will want to re-do this using mle2 or a Bayesian model, and then plot confidence intervals for the parameter estimates


background.model <- glmer(ectopic_present ~ 1 + (1|genotype), family = "binomial", data=penetrance.data)


summary(background.model)
ranef(background.model)
boot::inv.logit(2.048)

null.model <- glm(ectopic_present ~ 1, family="binomial", data=penetrance.data)


###########
#Fit confidence intervals individually

penetrance.data.tabulated$estimate <- rep(0, 18)
penetrance.data.tabulated$ci.min <- rep(0,18)
penetrance.data.tabulated$ci.max <- rep(0,18)
penetrance.data.tabulated$p.est <- rep(0,18)
penetrance.data.tabulated$p.min <- rep(0,18)
penetrance.data.tabulated$p.max <- rep(0,18)

tunes <- c( 2.5, #1
			0.0008, #2
			2.5, #3
			0.0008, #4
			2.5, #5
			0.0008, #6
			2.5, #7
			0.0008, #8
			0.0008, #9
			0.0008, #10
			2.5, #11
			0.0008, #12
			2.5, #13
			2.5, #14
			2.5, #15
			2.5, #16
			2.5, #17
			0.0008 #18
			)

for (i in 1:18) {
	hits <- penetrance.data.tabulated[i,2]
	trials <- penetrance.data.tabulated[i,3]
	cur.genotype <- as.character(penetrance.data.tabulated[i,1])
	cur.data <- penetrance.data[as.character(penetrance.data$genotype)==cur.genotype,]
	cur.model <- MCMClogit(ectopic_present ~ 1, data=cur.data, burnin=3000, mcmc=100000, thin=100, tune=tunes[i], verbose=10000)
	penetrance.data.tabulated$estimate[i] <- summary(cur.model)$statistics[1]
	penetrance.data.tabulated$ci.min[i] <- summary(cur.model)$quantiles[1]
	penetrance.data.tabulated$ci.max[i] <- summary(cur.model)$quantiles[5]
	penetrance.data.tabulated$p.min[i] <- plogis(penetrance.data.tabulated$ci.min[i])
	penetrance.data.tabulated$p.max[i] <- plogis(penetrance.data.tabulated$ci.max[i])
	penetrance.data.tabulated$p.est[i] <- penetrance.data.tabulated[i,2] / penetrance.data.tabulated[i,3]
}

penetrance.data.tabulated <- penetrance.data.tabulated[order(penetrance.data.tabulated$p.est),]
penetrance.data.tabulated$genotype <- as.factor(as.character(penetrance.data.tabulated$genotype))

pdf(file="../outputs/Figure3_RAL.pdf", width=6, height=10)
par(mfrow=c(2,1))

par(mar=c(5,5,1,1))
plotCI(x=1:18, y=penetrance.data.tabulated$p.est, 
    ui=penetrance.data.tabulated$p.max, 
    li=penetrance.data.tabulated$p.min, 
    cex.lab = 1.5, pch = 20, lwd = 2, cex = 1.5,
    xlim = c(0.5,18.5), ylim = c(0, 1.1), xaxt = "n", xlab = NA, 
    ylab = expression(paste(italic(scd)^1, " penetrance")))
    
xlabels <- as.character(penetrance.data.tabulated$genotype)
#xlabels[1] <- "scd[1] stock"
xlabels[1] <- expression(italic(scd)^1)

axis(side=1, at=1:18, labels=xlabels, las=2, cex.axis=0.8)
text(x=0.7, y=1.05, "A", cex=2)
								
#####################################################
#2. Expressivity of ts2 comb differs among lines
#####################################################
ts2.size.data <- scd.data[scd.data$t2_teeth>0,c("genotype", "t2_teeth")]
ts2.size.model <- glm(t2_teeth ~ genotype, data=ts2.size.data)

chi2.ts2 <- ts2.size.model$null.deviance - ts2.size.model$deviance
chi2.ts2
pchisq(chi2.ts2, df=1, lower.tail=F)
pchisq(chi2.ts2, df=length(ts2.size.model$coefficients)-1, lower.tail=F)



ts2.confints <- confint(ts2.size.model)
penetrance.data.tabulated$ts2.tooth.est <- rep(0,18)
penetrance.data.tabulated$ts2.tooth.ci.min <- rep(0,18)
penetrance.data.tabulated$ts2.tooth.ci.max <- rep(0,18)
for (i in 1:18) {
	cur.genotype <- as.character(penetrance.data.tabulated$genotype[i])
	cur.genotype <- ifelse(cur.genotype=="5070", "Intercept", cur.genotype)
	confint.row <- grep(cur.genotype, names(ts2.confints[,1]))
	penetrance.data.tabulated$ts2.tooth.ci.min[i] <- ts2.confints[confint.row,1]
	penetrance.data.tabulated$ts2.tooth.ci.max[i] <- ts2.confints[confint.row,2]
	penetrance.data.tabulated$ts2.tooth.est[i] <- coef(ts2.size.model)[confint.row]
	print('ok1')
	if (cur.genotype == "Intercept") {
		print('ok2')
		intercept <- penetrance.data.tabulated$ts2.tooth.est[i]
	}
	else {
		print('ok3')
		penetrance.data.tabulated$ts2.tooth.ci.min[i] <- penetrance.data.tabulated$ts2.tooth.ci.min[i] + intercept
		penetrance.data.tabulated$ts2.tooth.ci.max[i] <- penetrance.data.tabulated$ts2.tooth.ci.max[i] + intercept
		penetrance.data.tabulated$ts2.tooth.est[i] <- penetrance.data.tabulated$ts2.tooth.est[i] + intercept
	}
}

par(mar=c(5,5,1,1))

plotCI(x = 1:18, y = penetrance.data.tabulated$ts2.tooth.est, 
    ui = penetrance.data.tabulated$ts2.tooth.ci.max, 
    li = penetrance.data.tabulated$ts2.tooth.ci.min, 
    xlim = c(0.5,18.5), ylim = c(0, 5), xaxt = "n", xlab = NA, 
    cex.lab = 1.5, pch = 20, lwd = 2, cex = 1.5,
    #ylab="scd[1] expressivity \n (mean ectopic sex comb tooth number)")
    ylab = expression(paste(italic(scd)^1, " expressivity")))
        
xlabels <- as.character(penetrance.data.tabulated$genotype)
xlabels[1] <- expression(italic(scd)^1)
axis(side=1, at=1:18, labels=xlabels, las=2, cex.axis=1.1)

text(x=0.7, y=4.77, "B", cex=2)

dev.off()

#Plot distributions by genotype
#Are they normal-ish?
#If so, ok							

histogram( ~ t2_teeth | genotype, data=scd.data)


sessionInfo()