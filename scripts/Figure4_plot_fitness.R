
library(lme4)
library(glmmTMB)
library(sciplot)


fitness.data <- read.csv('../data/compensation_4_30.csv', na.strings="NA", stringsAsFactors = TRUE)

fitness.data$prop.scd <- fitness.data$scd / (fitness.data$scd + fitness.data$wt)

rep.colors <- c("#EE2222", "#55AA00", "#1133FF", "#992299")

fitness.data$plot.col <- rep.colors[fitness.data$rep]

#par(mar=c(4,4,1,1))

ci.95 <- function(x) {
	mean.x <- mean(x)
	se.x <- se(x)
	return(c(mean.x-2*se.x, mean.x+2*se.x))	
}

pdf(file="../outputs/Figure4_fitness.pdf", width = 5, height = 4)

par(mar=c(5,5,1,1))

x.vals <- c(0,1,2,3,4,5,9)
y.vals <- rep(NA, length(x.vals))
y.mins <- rep(NA, length(x.vals))
y.maxs <- rep(NA, length(x.vals))
for (n in 1:length(x.vals)) {
	cur.gen <- x.vals[n]
	cur.data <- fitness.data[fitness.data$gen==cur.gen,]
	y.vals[n] <- mean(cur.data$prop.scd)
	cur.ci <- ci.95(cur.data$prop.scd)
	y.mins[n] <- cur.ci[1]
	y.maxs[n] <- cur.ci[2]
}

plot(x=x.vals, y=y.vals, xlab="Generation", pch=20, cex=2, lwd=2, 
    cex.lab=1.5, ylab=expression(paste("Frequency of ", italic(scd)^1, " males")),
    ylim=c(0, 1))
for (n in 1:length(x.vals)) {
	arrows(x.vals[n], y.mins[n], x.vals[n], y.maxs[n], angle=90, length=0.05, code=3, lwd=2)
	#points(x=rep(x.vals[n], 2), y=c(y.mins[n], y.maxs[n]), type="both")
}

#original model -- estimate the intercept
gen.model.orig <- glm( cbind(scd, wt) ~ gen, 
   family = "binomial", 
   data = fitness.data)


#fix the intercept at a probability of 0.7
fitness.data$dummy.offset <- qlogis(0.7) # vector for offset

gen.model <- glm( cbind(scd, wt) ~ 0 + gen + offset(dummy.offset),
   family = "binomial", 
   data = fitness.data)

summary(gen.model)

#should there be a random effect of replicate in here too?


gen.model.2 <- glmer( cbind(scd, wt) ~ 
    0 + gen + ( 0 + gen| rep) + offset(dummy.offset), 
    family = "binomial", data = fitness.data)
    
summary(gen.model.2)
plogis(-0.19691)

#in the paper, give the estimate of the generation effect and p-value

#plot the best fit line [with confidence intervals] from the GLM above
predict.data <- data.frame(gen=seq(from = 0, to = 10, by = 1))
predict.data$dummy.offset <- qlogis(0.7)

y.pred <- predict(gen.model, newdata=predict.data, type="response", se.fit=T)
y.pred.top <- y.pred$fit + 2*y.pred$se
y.pred.bottom <- y.pred$fit - 2*y.pred$se

polygon(x=c(predict.data$gen, rev(predict.data$gen)), y=c(y.pred.top, rev(y.pred.bottom)),
border="#00000000", col="#00000031")

dev.off()

sessionInfo()
