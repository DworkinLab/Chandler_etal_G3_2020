library(ggplot2)
require(ggeffects)
library(sciplot)
require(bbmle)
require(effects)
require(glmmTMB) # model at the bottom.
require(MASS)
require(emmeans)
require(lme4)
library(grid)
library(simr)
library(car)


hss.index <- 1
lss.index <- 2
lv.index <- 3
wtc.index <- 4

plot.colors <- c(	"#FF0000", #HSS
					"#1122CC", #LSS
					"#449911", #LV
					"#993399" #WTC
			)
fill.colors <- c(	"#FF000033", #HSS
					"#1122CC33", #LSS
					"#44991133", #LV
					"#99339933" #WTC
			)


text.pars <- gpar(cex=1.5, fontface="bold")

#Setting up data set and sex comb teeth plots

sexcomb <- read.csv("../data/sex_comb_data_uptogen24.csv", header = T)

str(sexcomb) #check that treatment is a factor, generation is an interger and change replicate to a factor
sexcomb$rep <- factor(sexcomb$rep)
sexcomb$treatment <- relevel(sexcomb$treatment, "WTC")
#sexcomb <- na.omit(sexcomb)
sexcomb$individualid <- factor(paste(sexcomb$treatment, sexcomb$rep, sexcomb$gen, sexcomb$id, sep="."))

sexcomb$gen0 <- sexcomb$gen - 1 # Since our first measurements are in generation one, we want to make sure we are not extrapolating intercept for generation effect

####################################
#Primary comb
#First: run models for primary comb using glmmTMB, only HSS and LSS


#NOTES HERE:
#The model runs just fine with only HSS and LSS treatments
primary_evolve_model_glmm <- glmmTMB(primary_teeth ~ gen0*treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = poisson, 
   #ziformula = ~gen, # no evidence of zero inflation
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])

summary(primary_evolve_model_glmm)
plot(allEffects(primary_evolve_model_glmm))

plot(emmeans(primary_evolve_model_glmm, "gen0", 
    by = "treatment"))
     
#there is clearly no effect of generation


#However, when we add the WTC or LV treatments, we get warnings about model convergence
# --either non-positive-definite Hessian matrix, or
# --extreme or very small eigenvalues detected;
# but in both cases, the actual effect estimates are about the same (2.406 for the intercept, essentially 0 for the others,
#  but the std errors for those estimates are NA)
primary_evolve_model_glmm <- glmmTMB(primary_teeth ~ gen0*treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = poisson, 
   #ziformula = ~gen, # no evidence of zero inflation
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS" | sexcomb$treatment=="LV" ,])

summary(primary_evolve_model_glmm)
plot(allEffects(primary_evolve_model_glmm))



# These warnings go away, and we get standard error estimates, when we remove the random effects from the model
#  so the convergence problem probably has to do with difficulty in estimating random effects?
# (Is it too difficult to estimate the rep effect if there are only 2 reps for a treatment, or if 
#   some of the reps are missing data? e.g., there are only
#   2 reps for WTC, and the 3rd rep for LV is also missing some of the generations)
primary_evolve_model_glmm <- glmmTMB(primary_teeth ~ gen0*treatment 
    + (1|individualid), 
    family = poisson, 
   #ziformula = ~gen, # no evidence of zero inflation
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS" | sexcomb$treatment=="WTC" | sexcomb$treatment=="LV",])

summary(primary_evolve_model_glmm)
plot(allEffects(primary_evolve_model_glmm))

#with just the other random effect, the model works too
primary_evolve_model_glmm <- glmmTMB(primary_teeth ~ gen0*treatment 
    + (1|treatment:rep), 
    family = poisson, 
   #ziformula = ~gen, # no evidence of zero inflation
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS" | sexcomb$treatment=="WTC" | sexcomb$treatment=="LV",])

summary(primary_evolve_model_glmm)
plot(allEffects(primary_evolve_model_glmm))

#in any case the numbers are almost identical no matter which random effect we keep,
# and there's clearly no real effect of anything for
# primary comb tooth number






#Plot

as.data.frame(effect(c("gen0", "treatment"), primary_evolve_model_glmm ))

new_dat <- data.frame(gen0 = c(0:23, 0:23), 
    treatment = c(rep("LSS", 24), rep("HSS", 24), rep("WTC", 24), rep("LV", 24)),
    individualid = NA, rep = NA)

#Do predictions based on the model fit
predictTrial_glmmTMB <- predict(primary_evolve_model_glmm, 
   new_dat, interval = "confidence",
   allow.new.levels=TRUE, type = "response", se.fit = TRUE)


length(predictTrial_glmmTMB$fit)
head(predictTrial_glmmTMB$fit)
cbind(new_dat, predictTrial_glmmTMB)

length(predictTrial_glmmTMB$se)

summary(primary_evolve_model_glmm)    
fixef(primary_evolve_model_glmm)

#Run similar model, but using lme4

primary_evolve_model_glmm_lme4 <- glmer(primary_teeth ~ gen0*treatment 
    + (1|individualid), 
    family = poisson, 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS" | sexcomb$treatment=="WTC" | sexcomb$treatment=="LV",])
##plot(allEffects(primary_evolve_model_glmm_lme4))

as.data.frame(Effect(c("treatment", "gen0"), primary_evolve_model_glmm_lme4))
    
summary(primary_evolve_model_glmm_lme4)

#take-home here is that they give almost identical results, in spite of possible warnings about not converging

#Make plots

predictions <- ggemmeans(primary_evolve_model_glmm, terms=c("gen0", "treatment")) # on link scale. transform back via exp(). Very similar to using effects
predictions$x <- predictions$x + 1 #move back to original x-values for plotting

myplot.1 <- ggplot(predictions, aes(x, predicted, color=group)) + 
 geom_line() + geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, color=NA) + 
 scale_color_manual(values=c(plot.colors[wtc.index], plot.colors[hss.index], plot.colors[lss.index], plot.colors[lv.index])) + 
 labs(x="Generation", y="Number of primary sex comb teeth") + 
 scale_fill_manual(values=c(fill.colors[wtc.index], fill.colors[hss.index], fill.colors[lss.index], fill.colors[lv.index])) +
 theme(legend.title=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(), panel.border=element_rect(color="black", fill=NA)) +
 theme(plot.margin=unit(c(1,1,1,1), "cm"))

print(myplot.1) 




## Added Feb 13th, 2020 - Reviewer really wants a common intercept model for these.
## Of note, we disagreed. While LSS and HSS could be fit sensibly with a common intercept, it makes no sense for the other two treatment levels.
## So I have fit a model, where there is an intercept fit, and a generation:treatment effect
## Only for the HSS and LSS treatments as there is no reason to assume a common intercept for LV or WTC
## Please also note that this will generate two slopes (one for HSS and LSS), not just a slope for HSS and the LSS contrast. This seems to behave better computationally (less convergence issues). These remain on the log link scale.
 
primary_evolve_model_glmm_lme4_int <- glmer(primary_teeth ~ 1 + gen0:treatment 
    + (1|individualid), 
    family = poisson, 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])
    
##plot(allEffects(primary_evolve_model_glmm_lme4_int))

as.data.frame(Effect(c("treatment", "gen0"), primary_evolve_model_glmm_lme4_int))
    
summary(primary_evolve_model_glmm_lme4_int)
Anova(primary_evolve_model_glmm_lme4_int)

## End added Feb 13th 2020


####################################
#Primary comb - gaps
#First: run models for primary comb gaps using glmmTMB, only HSS and LSS



sexcomb$primary_gaps_binom <- ifelse(sexcomb$primary_gaps >= 1, 1, 0)

#The model runs just fine with only HSS and LSS treatments [the one we will ultimately plot for the paper]
gaps_evolve_model_glmm_plot <- glmmTMB(primary_gaps_binom ~ gen0*treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = binomial(link="logit"), 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])

summary(gaps_evolve_model_glmm_plot)
plot(allEffects(gaps_evolve_model_glmm_plot))
#according to the model there is an effect of generation; the probability of having a gap defect in the primary
# comb decreases a bit each generation
#the interaction is not significant; however, it LOOKS very different on the plot -- much steeper for HSS



#When we add the LV treatment...
gaps_evolve_model_glmm <- glmmTMB(primary_gaps_binom ~ gen0*treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = binomial(link="logit"), 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS" | sexcomb$treatment=="LV" ,])

summary(gaps_evolve_model_glmm)
plot(allEffects(gaps_evolve_model_glmm))
#largely the result is the same, but the LV main effect is significant
#it appears that the LV treatment has a lower probability of gaps overall


#And we add the WTC treatment... (even though defects are not expected here).

gaps_evolve_model_glmm <- glmmTMB(primary_gaps_binom ~ gen0*treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = binomial(link="logit"), 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS" | sexcomb$treatment=="LV" | sexcomb$treatment=="WTC",])

summary(gaps_evolve_model_glmm)
plot(allEffects(gaps_evolve_model_glmm))
#the results are different, not surprisingly (with wild type control after all). the generation main effect is gone, but this is because it is estimating the slope of the response for the WT which really has no defects.

#the HSS effect is significant (higher probability of having gaps overall) and the
# LSS effect is marginally significant (p=0.07)
#again, it LOOKS like probability of gaps decreases in the HSS treatment, but the interaction
# term is not significant






#Plot

as.data.frame(effect(c("gen0", "treatment"), gaps_evolve_model_glmm_plot ))

new_dat <- data.frame(gen0 = c(0:23, 0:23), 
    treatment = c(rep("LSS", 24), rep("HSS", 24)),
    individualid = NA, rep = NA)

#Do predictions based on the model fit
predictTrial_glmmTMB <- predict(gaps_evolve_model_glmm_plot, 
   new_dat, interval = "confidence",
   allow.new.levels=TRUE, type = "response", se.fit = TRUE)


length(predictTrial_glmmTMB$fit)
head(predictTrial_glmmTMB$fit)
cbind(new_dat, predictTrial_glmmTMB)

length(predictTrial_glmmTMB$se)

summary(gaps_evolve_model_glmm_plot)    
fixef(gaps_evolve_model_glmm_plot)

#Run similar model, but using lme4

gaps_evolve_model_glmm_lme4 <- glmer(primary_gaps_binom ~ gen0*treatment 
    + (1|individualid), 
    family = binomial(link="logit"), 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS" | sexcomb$treatment=="WTC" | sexcomb$treatment=="LV",])
plot(allEffects(gaps_evolve_model_glmm_lme4))

as.data.frame(Effect(c("treatment", "gen0"), gaps_evolve_model_glmm_lme4))
    
summary(gaps_evolve_model_glmm_lme4)

#Make plots

predictions <- ggemmeans(gaps_evolve_model_glmm_plot, terms=c("gen0", "treatment")) # are these predictions on link scale?
#they look like they are already probabilities
predictions$x <- predictions$x + 1 #move back to original x-values for plotting

myplot.2 <- ggplot(predictions, aes(x, predicted, color=group)) + 
 geom_line() + geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, color=NA) + 
 scale_color_manual(values=c(plot.colors[hss.index], plot.colors[lss.index])) + 
 labs(x="Generation", y="Proportion with primary comb defects") + 
 scale_fill_manual(values=c(fill.colors[hss.index], fill.colors[lss.index])) + 
 theme(legend.title=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(), panel.border=element_rect(color="black", fill=NA)) +
 theme(plot.margin=unit(c(1,1,1,1), "cm"))
print(myplot.2)
 



## Added February 13th, 2020. Common intercept model (fit for just HSS and LSS, see comments above)


gaps_evolve_model_glmm_int <- glmmTMB(primary_gaps_binom ~ 1 + gen0:treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = binomial(link="logit"), 
#   ziformula = ~gen, # no evidence of zero inflation
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])

summary(gaps_evolve_model_glmm_int)
plot(allEffects(gaps_evolve_model_glmm_int))


gaps_evolve_model_glmm_lme4_int <- glmer(primary_gaps_binom ~ 1 + gen0:treatment 
    + (1|individualid), 
    family = binomial(link="logit"), 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])
    
plot(allEffects(gaps_evolve_model_glmm_lme4_int))

as.data.frame(Effect(c("treatment", "gen0"), gaps_evolve_model_glmm_lme4))
    
summary(gaps_evolve_model_glmm_lme4_int)
# same as glmmTMB

## Plot for common intercept model
predictions_int <- ggemmeans(gaps_evolve_model_glmm_int, terms=c("gen0", "treatment")) # are these predictions on link scale?
#they look like they are already probabilities
predictions$x <- predictions$x + 1 #move back to original x-values for plotting

myplot.2b <- ggplot(predictions_int, aes(x, predicted, color=group)) + 
 geom_line() + geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, color=NA) + 
 scale_color_manual(values=c(plot.colors[hss.index], plot.colors[lss.index])) + 
 labs(x="Generation", y="Proportion with primary comb defects") + 
 scale_fill_manual(values=c(fill.colors[hss.index], fill.colors[lss.index])) + 
 theme(legend.title=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(), panel.border=element_rect(color="black", fill=NA)) +
 theme(plot.margin=unit(c(1,1,1,1), "cm"))
print(myplot.2b)


## end new bit, for V3, common intercept for defects in primary comb.



####################################
#Ectopic comb
#First: run models for ectopic comb using glmmTMB, only HSS and LSS

secondary_evolve_model_glmm <- glmmTMB(secondary_teeth ~ gen0*treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = poisson, 
   #ziformula = ~gen, # no evidence of zero inflation
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])



summary(secondary_evolve_model_glmm)
plot(allEffects(secondary_evolve_model_glmm))
#take-home message here is:
# there is an effect of generation (i.e., ectopic sex comb tooth number decreases)
# but slope does not differ between HSS and LSS



secondary_LV_evolve_model_glmm <- glmmTMB(secondary_teeth ~ gen0 + (1|individualid), 
    family = poisson, 
   #ziformula = ~gen, # no evidence of zero inflation
    data = sexcomb[sexcomb$treatment=="LV",])

summary(secondary_LV_evolve_model_glmm)
plot(allEffects(secondary_LV_evolve_model_glmm))
#take-home message here is:
# no effect of gen when LV is considered on its own


as.data.frame(effect(c("gen0", "treatment"), secondary_evolve_model_glmm ))

new_dat <- data.frame(gen0 = c(0:23, 0:23), 
    treatment = c(rep("LSS", 24), rep("HSS", 24)),
    individualid = NA, rep = NA)

#Do predictions based on the model fit
predictTrial_glmmTMB <- predict(secondary_evolve_model_glmm, 
   new_dat, interval = "confidence",
   allow.new.levels=TRUE, type = "response", se.fit = TRUE)


length(predictTrial_glmmTMB$fit)
head(predictTrial_glmmTMB$fit)
cbind(new_dat, predictTrial_glmmTMB)

length(predictTrial_glmmTMB$se)

summary(secondary_evolve_model_glmm)    
fixef(secondary_evolve_model_glmm)

#Run similar model, but using lme4 - also power analysis.

secondary_evolve_model_glmm_lme4 <- glmer(secondary_teeth ~ gen*treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = poisson, 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])
plot(allEffects(secondary_evolve_model_glmm_lme4))

as.data.frame(Effect(c("treatment", "gen"), secondary_evolve_model_glmm_lme4))
    
summary(secondary_evolve_model_glmm_lme4)

#take home here is that they give almost identical results, in spite of possible warnings about not converging



#Make plots

predictions <- ggemmeans(secondary_evolve_model_glmm, terms=c("gen0", "treatment")) # on link scale. transform back via exp(). Very similar to using effects
predictions$x <- predictions$x + 1 #move back to original x-values for plotting


myplot.3 <- ggplot(predictions, aes(x, predicted, color=group)) + 
 geom_line() + geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, color=NA) + 
 scale_color_manual(values=c(plot.colors[hss.index], plot.colors[lss.index])) + 
 labs(x="Generation", y="Number of ectopic sex comb teeth") + 
 scale_fill_manual(values=c(fill.colors[hss.index], fill.colors[lss.index])) +
 theme(legend.title=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(), panel.border=element_rect(color="black", fill=NA)) +
 theme(plot.margin=unit(c(1,1,1,1), "cm"))

print(myplot.3)


## Note: this part added Feb 13th, 2020. Fitting a common intercept for HSS and LSS as a reviewer wanted this.
secondary_evolve_model_glmm_int <- glmmTMB(secondary_teeth ~ 1 + gen0:treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = poisson, 
   #ziformula = ~gen, # no evidence of zero inflation
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])
    
# model.matrix(secondary_evolve_model_glmm_int) #just checking it is producing what I want 

summary(secondary_evolve_model_glmm_int)
confint(secondary_evolve_model_glmm_int)
Anova(secondary_evolve_model_glmm_int)
plot(allEffects(secondary_evolve_model_glmm_int))



secondary_evolve_model_glmm_lme4_int <- glmer(secondary_teeth ~ 1 + gen0:treatment 
    + (1|treatment:rep) + (1|individualid), 
    family = poisson, 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])
    
plot(allEffects(secondary_evolve_model_glmm_lme4_int))

as.data.frame(Effect(c("treatment", "gen0"), secondary_evolve_model_glmm_lme4_int))
    
summary(secondary_evolve_model_glmm_lme4_int)

#Plot for the common intercept model
predictions <- ggemmeans(secondary_evolve_model_glmm_int, terms=c("gen0", "treatment")) # on link scale. transform back via exp(). Very similar to using effects
predictions$x <- predictions$x + 1 #move back to original x-values for plotting


myplot.3b <- ggplot(predictions, aes(x, predicted, color=group)) + 
 geom_line() + geom_ribbon(aes(ymin=conf.low, ymax=conf.high, fill=group), alpha=0.2, color=NA) + 
 scale_color_manual(values=c(plot.colors[hss.index], plot.colors[lss.index])) + 
 labs(x="Generation", y="Number of ectopic sex comb teeth") + 
 scale_fill_manual(values=c(fill.colors[hss.index], fill.colors[lss.index])) +
 theme(legend.title=element_blank(), panel.background=element_blank(), panel.grid.major=element_blank(), panel.border=element_rect(color="black", fill=NA)) +
 theme(plot.margin=unit(c(1,1,1,1), "cm"))

print(myplot.3b)




##### End New part (Feb 13th, 2020)




####################################
#Make actual PDF 
pdf(file="../outputs/Figure5_EE.pdf", width=7, height=15)

vp.1 <- viewport(x=0, y=1, width=1, height=1/3, just=c("left", "top"))
print(myplot.1, vp=vp.1)
grid.text(label="A", x=0.01, y=0.98, just=c("left", "top"), gp=text.pars, vp=vp.1)

vp.2 <- viewport(x=0, y=2/3, width=1, height=1/3, just=c("left", "top"))
print(myplot.2b, vp=vp.2)
grid.text(label="B", x=0.01, y=0.98, just=c("left", "top"), gp=text.pars, vp=vp.2)

vp.3 <- viewport(x=0, y=1/3, width=1, height=1/3, just=c("left", "top"))
print(myplot.3b, vp=vp.3)
grid.text(label="C", x=0.01, y=0.98, just=c("left", "top"), gp=text.pars, vp=vp.3)

dev.off()

##### Power analysis ##### 

# For number of secondary teeth.
#For the power analysis I needed to use a slightly modified model or otherwise simr seems to have some difficulties (estimates a power of 0 for all effects, even when the magnitudes are 1000X greater than observed.)

# Please note that I have used a slightly different model specification for the power analysis so that I can directly contrast the slopes of the high vs low sexual selection treatments for the power analysis. 
# I have left the default factor levels on. So it fits the HSS and looks at LSS as a contrast to HSS. However, for our purposes the results (estimates, SE, CI and p) and only differ in the sign of the estimate.


secondary_evolve_model_glmm_lme4 <- glmer(secondary_teeth ~ 1 + gen0*treatment - treatment
    + (1|treatment:rep) , 
    family = poisson, 
    data = sexcomb[sexcomb$treatment=="HSS" | sexcomb$treatment=="LSS",])


summary(secondary_evolve_model_glmm_lme4)
plot(allEffects(secondary_evolve_model_glmm_lme4))

# Retrospective power analysis, What is the power to detect an effect, assuming our estimates are reasonable.

# observed power simulation for gen:treat.

fixef(secondary_evolve_model_glmm_lme4)["gen0"]
doTest(secondary_evolve_model_glmm_lme4, fixed("gen0", "z"))
doTest(secondary_evolve_model_glmm_lme4, fixed("gen0:treatmentLSS", "z"))


doTest(secondary_evolve_model_glmm_lme4, fcompare(~ 1 + gen)) # i.e. this checks whether a common slope model is a poorer fit than the seperate slopes model

# Making sure I can get it to work (this is a retrospective power analysis, so not particularly meaningful, but just confirms things are working)
sec_teeth_power <- powerSim(secondary_evolve_model_glmm_lme4,
    fixed("gen:treatmentLSS", "z"), nsim = 100)

print(sec_teeth_power)
summary(sec_teeth_power)


# testing with a large effect size, again just to confirm things are working for the power simulations.

power_mod_test <- secondary_evolve_model_glmm_lme4 # copy the model

fixef(power_mod_test)["gen0:treatmentLSS"] <- -100

sec_teeth_power2 <- powerSim(power_mod_test,
    fcompare(~ 1 + gen0), nsim = 100)

summary(sec_teeth_power2)
print(sec_teeth_power2) # pretty much 100% power as expected.


# try up to 30 reps. That is, keep the other effects at observed levels, but increase the number of replicate lineages.

# create a copy of the model to work with
power_mod_effect_rep <- secondary_evolve_model_glmm_lme4
fixef(power_mod_effect_rep)["gen0:treatmentLSS"] <- -0.00380 # keeping it approximately the same as observed ON PURPOSE!!! I only want to ater number of simulated replicate lineages.

power_more_reps <- extend(power_mod_effect_rep, along = "rep", n = 30)

pow_more_reps_out <- powerCurve(power_more_reps, along = "rep", 
    nsim = 100, fixed("gen:treatmentLSS", "z")) # run with nsim = 100 or more
 
summary(pow_more_reps_out)
plot(pow_more_reps_out)

# number of generations
power_more_gen <- secondary_evolve_model_glmm_lme4

fixef(power_more_gen)["gen0:treatmentLSS"] <- -0.00380 # again keeping it close to observed effect to understand what increasing number of generations does in terms of uncertainty of estimates.

power_gen <- extend(power_more_gen, along = "gen0", n = 36)

pow_gen_out <- powerCurve(power_gen, along = "gen0", 
    nsim = 200, fixed("gen0:treatmentLSS", "z")) # run with nsim = 1000
summary(pow_gen_out)
plot(pow_gen_out)

# How about if we alter the magnitude of the effect (i.e. LSS vs. HSS)
# For some reason extend was not working for this effect.
# Thus I did the power analysis at three magnitudes. Close to observed, close to the slope for the other treatment and then approximately double that. Not surprisingly it shows a steady increase in power.
power_more_mag1 <- secondary_evolve_model_glmm_lme4

fixef(power_more_mag1)["gen0:treatmentLSS"] <- -0.00380

sec_teeth_power1 <- powerSim(power_more_mag1,
    fcompare(~ 1 + gen0), nsim = 100)

    
summary(sec_teeth_power1)
print(sec_teeth_power1)



power_more_mag2 <- secondary_evolve_model_glmm_lme4

fixef(power_more_mag2)["gen0:treatmentLSS"] <- -0.0067 # of same magnitude as slope being contrasted

sec_teeth_power2 <- powerSim(power_more_mag2,
    fcompare(~ 1 + gen0), nsim = 100)

    
summary(sec_teeth_power2)
print(sec_teeth_power2)


power_more_mag3 <- secondary_evolve_model_glmm_lme4

fixef(power_more_mag3)["gen0:treatmentLSS"] <- -0.012 # of double magnitude as slope being contrasted

sec_teeth_power3 <- powerSim(power_more_mag3,
    fcompare(~ 1 + gen0), nsim = 100)

    
summary(sec_teeth_power3)
print(sec_teeth_power3)


sessionInfo()