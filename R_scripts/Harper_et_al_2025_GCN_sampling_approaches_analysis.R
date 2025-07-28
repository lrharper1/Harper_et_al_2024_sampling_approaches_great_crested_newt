#' ---
#' Title: "Extending sampling approaches for great crested newt (Triturus cristatus) eDNA monitoring"
#' Author: "Lynsey Harper"
#' Date: "28th July 2025"
#' ---
#' 
#' Paired water samples for ethanol precipitation and filtration were collected 
#' from ponds once a month from April to October by Atkins for NatureMetrics. 
#' The water sampling protocol set forth by Biggs et al. (2014) was used for 
#' both eDNA capture methods with minor modifications for filtration, where 
#' 20 x 125 mL subsamples were collected at equidistant intervals around the 
#' pond perimeter and pooled into a single sampling bag for homogenisation, 
#' following which as much water as possible was filtered. We aimed to include 
#' 17 positive and 3 negative ponds for great crested newt (GCN) each month, and 
#' field blanks (mineral water) were included from May onwards. Conventional 
#' population size assessments for each pond were also performed for comparison 
#' to eDNA results. 
#' 
#' The resulting data will be analysed to examine:
#' 
#' 1. The number of positive ponds produced each month by ethanol precipitation
#'    and filtration.
#' 2. The eDNA score (number of positive qPCR replicates) for each pond in each 
#'    month produced by ethanol precipitation and filtration.
#' 3. The number of positive ponds in-season (April to June) vs. out-of-season
#'    (July to October) produced by ethanol precipitation and filtration.
#' 4. The eDNA score for each pond in-season (April to June) vs. out-of-season
#'    (July to October) produced by ethanol precipitation and filtration
#' 5. Whether volume of water filtered influenced the eDNA scores for ponds 
#'    (and thus detection).
#' 6. Whether population size influenced the eDNA scores for ponds (and thus 
#'    detection).
#' 7. The probability of eDNA detection when imperfect detection is accounted
#'    for via hierachical occupancy modelling.
#'  


## Clear R memory.
rm(list=ls())

## Set working directory to location of R script.
# On first time running code: install.packages("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Check working directory.
getwd()

## Install and load required packages.
p <- c("plyr", "tidyverse", "ggpubr", "patchwork", "ggrepel", "scales", "ggh4x", 
       "sf", "mapview", "reshape2", "TMB", "glmmTMB", "DHARMa", "easystats", 
       "arm", "emmeans", "effects", "boot")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, 
                                          repos = "https://cran.ma.imperial.ac.uk/",
                                          dependencies = TRUE)
lapply(p, require, character.only = TRUE)

## To ensure reproducibility, print details about the version of R being used 
## for analysis.
sessionInfo()



###########
# DATASET #
###########

## Import data.
GCN.dat <- read.csv("Data/Field_study_data.csv",
                    header = TRUE)

## Remove redundant columns and field blanks. Make necessary columns factors.
## Order 'Month' factor chronologically and change 'Kit_type' to 'eDNA_capture'.
GCN.dat <- GCN.dat %>%
  dplyr::select(-c(Kit_ID, Sampled, Arrived, Inhibition, Degradation, Notes)) %>%
  drop_na(Pond_ID) %>%
  droplevels() %>%
  mutate_if(is.character, as.factor) %>%
  mutate_at("Pond_number", as.factor) %>%
  mutate_at(c("Volume", "eDNA_score", "Average_Cq"), as.numeric) %>%
  mutate(Month = fct_relevel(Month, 
                             "April", "May", "June", "July", "August", 
                             "September", "October")) %>%
  rename(eDNA_capture = Kit_type)


################
# MAP OF PONDS #
################

## Transform coordinates to geometric points for mapping.
ponds <- GCN.dat %>%
  dplyr::select(Pond_number, Previous_status, Latitude, Longitude) %>%
  distinct()

sites <- st_as_sf(ponds, coords = c("Longitude", "Latitude"), crs = 4326)

## Create and export interactive map of sites.
m1 <- mapview(sites, 
              zcol = "Previous_status", 
              col.regions = c("black", "orange", "grey"),
              alpha.regions = 0.7,
              cex = 5,
              layer.name = "Previous status",
              legend = FALSE)
m1

m2 <- mapview(sites, 
              zcol = "Previous_status", 
              col.regions = c("black", "orange", "grey"),
              alpha.regions = 0.7,
              cex = 5,
              layer.name = "Previous status",
              legend = FALSE)
m2

m3 <- mapview(sites, 
              zcol = "Previous_status", 
              col.regions = c("black", "orange", "grey"),
              alpha.regions = 0.7,
              cex = 5,
              layer.name = "Previous status",
              legend = FALSE)
m3

m4 <- mapview(sites, 
              zcol = "Previous_status", 
              col.regions = c("black", "orange", "grey"),
              alpha.regions = 0.7,
              cex = 5,
              layer.name = "Previous status",
              legend = FALSE)
m4

m <- leafsync::sync(m1, m2, m3, m4, sync = "none")
m

## Zoom to desired resolution on each panel to show all three sites at landscape
## scale and individual ponds at each site, then save as JPEG.


#############
# SUMMARIES #
#############

## Calculate number of positive samples for each eDNA capture method according
## to previous GCN status.
prev.survey <- GCN.dat %>% 
  group_by(Pond_number, Previous_status, eDNA_capture) %>%
  count(GCN_2022, sort = TRUE) %>%
  rename(Samples = n) %>%
  complete(GCN_2022)

## Plot number of positive samples for ponds in 2022 in relation to previous GCN
## status.
p1 <- ggplot(prev.survey, 
             aes(x = Pond_number, y = Samples, fill = eDNA_capture)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  scale_fill_manual(name = "eDNA capture method",
                    values = c("grey60", "white")) +
  scale_y_continuous(limit = c(0, 8), breaks = seq(0, 8, 1)) +
  labs(subtitle = "", x = "Pond", y = "Number of eDNA samples") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text.y = element_text(angle = 360),
        legend.position = "bottom",
        text = element_text(size = 20)) +
  facet_grid(Previous_status ~ GCN_2022)
p1

## Add titles to secondary x and y axes.
g1 <- ggarrange(p1)
g1 <- annotate_figure(g1, 
                      top = text_grob("2022 eDNA survey", color = "black", 
                                      size = 20, vjust = 1.5),
                      right = text_grob("Previous survey\n", color = "black", 
                                        size = 20, rot = 270, vjust = 0.9))

## Export plot.
#ggsave(filename = "Figures/GCN_previous_current_status.pdf", 
#       plot = g1, width = 25, height = 10, dpi = 300, units = "in")


## Add column for proportional eDNA score, then plot GCN detection with each 
## eDNA capture method for each pond over the course of the field study (April 
## to October).
GCN.dat.hm <- mutate(GCN.dat, prop_score = eDNA_score/12)
p2 <- ggplot(filter(GCN.dat.hm), 
             aes(x = eDNA_capture, y = fct_rev(Pond_number), 
                 fill = eDNA_score, group = eDNA_capture)) + 
  geom_tile(colour = "black", position = "dodge") +
  geom_text(aes(label = round(prop_score, 2)), colour = "white", size = 5,
            position = position_dodge(width = 1)) +
  scale_fill_gradientn(name = "eDNA score\n",
                       colours = c("black","black", "orange4","orange"),
                       values = c(0, 0.07, 0.08, 1), 
                       limits = c(0, 12),
                       breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
                       na.value = "grey",
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  labs(x = "eDNA capture method", y = "Pond") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"), 
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.background = element_blank(),
        ggh4x.facet.nestline = element_line(colour = "black"),
        text = element_text(size = 20),
        legend.key.size = unit(2, "lines"),
        legend.key.height = unit(2.5, "cm")) +
  facet_nested(. ~ Season + Month)
p2

## Export plot.
#ggsave(filename = "Figures/temporal_GCN_detection.pdf", 
#       plot = p2, width = 18, height = 13, dpi = 300, units = "in")


## Remove ponds that were previously negative for GCN, then remove column for
## previous status. Add column converting 'Positive/Negative' results to binary
## '1/0' results.
GCN2022 <- GCN.dat %>%
  filter(Previous_status == "Positive") %>%
  dplyr::select(-Previous_status) %>%
  mutate(GCN_binary = ifelse(GCN_2022 == "Positive", 1, 0)) %>%
  relocate(GCN_binary, .after = GCN_2022) %>%
  droplevels()



#################################
# eDNA CAPTURE METHOD AND MONTH #
#################################

#--------------------------#
# NUMBER OF POSITIVE PONDS #
#--------------------------#

## Test whether the difference in the number of positive ponds each month using
## ethanol precipitation and filtration is statistically significant using a 
## binomial GLMM.

## With interaction term between eDNA capture method and month:
glmm.pa.capture.month <- glmmTMB(GCN_binary ~ (1|Pond_number) + eDNA_capture 
                                 + Month + eDNA_capture:Month,
                                 family = binomial,
                                 data = GCN2022)
summary(glmm.pa.capture.month)
drop1(glmm.pa.capture.month, test = "Chi")


## The interaction term is not significant so run GLMM Without interaction term:
glmm.pa.capture.month <- glmmTMB(GCN_binary ~ (1|Pond_number) + eDNA_capture + Month,
                                 family = binomial,
                                 data = GCN2022)
summary(glmm.pa.capture.month)
drop1(glmm.pa.capture.month, test = "Chi")

## Run model diagnostics.
diagnose(glmm.pa.capture.month)                     # No issues detected by glmmTMB
check_model(glmm.pa.capture.month)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.pa.capture.month)           # Response follows Bernoulli distribution
check_singularity(glmm.pa.capture.month)            # Model fit is not singular
check_convergence(glmm.pa.capture.month)            # Model has converged
model_performance(glmm.pa.capture.month)            # Model has reasonable explanatory power
testDispersion(glmm.pa.capture.month, plot = T)     # Model is not overdispersed
testZeroInflation(glmm.pa.capture.month, plot = T)  # Data is not zero-inflated

## Get standardised residuals for model validation.
resid1 <- simulateResiduals(glmm.pa.capture.month, plot = T)
getResiduals(glmm.pa.capture.month)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid1, plot = F)
testQuantiles(resid1, plot = F)
testOutliers(resid1)
## Normally distributed.

## Check assumption of homogeneity of variances.
testCategorical(resid1, catPred = GCN2022$eDNA_capture)
testCategorical(resid1, catPred = GCN2022$Month)
## No heteroscedasticity present.

## Check dataset does not contain serial auto-correlation. This can arise 
## if there are unaccounted for relationships in data set or if there is 
## confounding effect of time or space.
recal.resid1 <- recalculateResiduals(resid1, group = GCN2022$Month)
testTemporalAutocorrelation(recal.resid1, time = unique(GCN2022$Month))

group.ponds <- aggregate(GCN2022[, 4:5], list(GCN2022$Pond_number), mean)
recal.resid2 <- recalculateResiduals(resid1, group = GCN2022$Pond_number)
testSpatialAutocorrelation(recal.resid2, group.ponds$Latitude, group.ponds$Longitude)
## Slight temporal autocorrelation and spatial autocorrelation detected.


## Post-hoc test for month.
glmm.pa.month.posthoc <- emmeans(glmm.pa.capture.month, 
                                 specs = "Month",
                                 type = "response")
multcomp::cld(glmm.pa.month.posthoc, Letters = letters)
pairs(glmm.pa.month.posthoc)
plot(glmm.pa.month.posthoc)

## Post-hoc test for eDNA capture method.
glmm.pa.capture.posthoc <- emmeans(glmm.pa.capture.month, 
                                   specs = "eDNA_capture",
                                   type = "response")
multcomp::cld(glmm.pa.capture.posthoc, Letters = letters)
pairs(glmm.pa.capture.posthoc)
plot(glmm.pa.capture.posthoc)


## Plot relationship between GCN presence, eDNA capture method and month using 
## predicted values from model.

## Calculate the number of GCN positive and negative ponds each month using 
## ethanol precipitation and filtration.
capture.ponds <- GCN2022 %>%
  group_by(eDNA_capture, Month) %>%
  count(GCN_2022) %>%
  rename(Ponds = n) %>%
  complete(GCN_2022)

## Plot number of GCN positive and negative ponds each month with each eDNA 
## capture method.
p3a <- ggplot(capture.ponds, 
             aes(x = eDNA_capture, y = Ponds, fill = GCN_2022)) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black") + 
  geom_text(aes(label = Ponds), 
            position = position_dodge(width = 0.9), 
            vjust = -1) +
  scale_fill_manual(name = "Great crested newt",
                    values = c("black", "orange")) +
  scale_y_continuous(limits = c(0, 20), breaks = seq(0, 20, 5)) +
  labs(title = "(a)",
       x = "eDNA capture method", y = "Number of ponds") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "bottom",
        text = element_text(size = 20)) +
  facet_grid(. ~ Month)
p3a

## Obtain predicted values from GLMM for the full data set.
pa.capture.month.fit <- data.frame(effect(c("Month", "eDNA_capture"), 
                                          glmm.pa.capture.month))

## Plot predicted against observed values for GCN detection in each month.
p3b <- ggplot() + 
  geom_errorbar(pa.capture.month.fit, 
                mapping = aes(x = eDNA_capture, ymin = lower, ymax = upper), 
                linewidth = 1, width = 0.3, color = "black") + 
  geom_point(pa.capture.month.fit, 
             mapping = aes(x = eDNA_capture, y = fit, fill = eDNA_capture), 
             size = 5, shape = 22) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  labs(title = "(b)",
       x = "eDNA capture method", y = "Probability of a positive eDNA result") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Month)
p3b


## Create multi-plot:
g2 <- ggarrange(p3a, p3b,
                nrow = 2, ncol = 1)

## Export plot.
#ggsave(filename = "Figures/GCN_detection_rate_month_capture.pdf", 
#       plot = g2, width = 18, height = 13, dpi = 300, units = "in")


#------------#
# eDNA SCORE #
#------------#

## Initial data exploration and modelling indicated that using eDNA score in its
## current form as the response variable resulted in problems with overdispersion,
## underdispersion or zero-inflation when Poisson, negative binomial,
## zero-inflated Poisson, and zero-inflated negative binomial GLMMs were employed.
## Even though eDNA score is integer count data, a Poisson or negative binomial 
## distribution is inappropriate because the data is bounded at 0 and 12 (i.e. a 
## sample cannot have an eDNA score lower or higher than this). This is because 
## the majority of eDNA scores are 0 and 12, resulting in data being inflated at 
## both boundaries. We will transform eDNA score to proportional eDNA score 
## (bounded from 0-1) to employ a binomial GLMM instead. Both zero-inflated and 
## regular binomial GLMMs were tested but the zero-inflated GLMM did not offer
## improved fit over the regular binomial GLMM, the AIC value for the regular
## binomial GLMM was lower, and the estimates produced by both GLMMs did not 
## differ. Therefore, results from the regular binomial GLMM are reported.

GCN2022 <- GCN2022 %>%
  mutate(prop_eDNA_score = round(eDNA_score/12, 2)) %>%
  relocate(prop_eDNA_score, .after = eDNA_score)

## Test whether the difference in proportional eDNA score for ponds each month 
## using ethanol precipitation and filtration is statistically significant using 
## a binomial GLMM.

## With interaction term between eDNA capture method and month:
glmm.score.capture.month <- glmmTMB(prop_eDNA_score ~ (1|Pond_number)
                                    + eDNA_capture + Month + eDNA_capture:Month,
                                    family = binomial(link="logit"),
                                    data = GCN2022)
summary(glmm.score.capture.month)
drop1(glmm.score.capture.month, test = "Chi")


## The interaction term is not significant so run GLMM Without interaction term:
glmm.score.capture.month <- glmmTMB(prop_eDNA_score ~ (1|Pond_number)
                                    + eDNA_capture + Month,
                                    family = binomial(link="logit"),
                                    data = GCN2022)
summary(glmm.score.capture.month)
drop1(glmm.score.capture.month, test = "Chi")

## Run model diagnostics.
diagnose(glmm.score.capture.month)                     # No issues detected by glmmTMB
check_model(glmm.score.capture.month)                  # Minor deviations in residuals
check_distribution(glmm.score.capture.month)           # Response follows binomial distribution
check_singularity(glmm.score.capture.month)            # Model fit is not singular
check_convergence(glmm.score.capture.month)            # Model has converged
model_performance(glmm.score.capture.month)            # Model has reasonable explanatory power
testDispersion(glmm.score.capture.month, plot = T)     # Model is underdispersed
testZeroInflation(glmm.score.capture.month, plot = T)  # Data is zero-inflated

## Get standardised residuals for model validation.
resid2 <- simulateResiduals(glmm.score.capture.month, plot = T)
getResiduals(glmm.score.capture.month)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid2, plot = F)
testQuantiles(resid2, plot = F)
testOutliers(resid2)
## Not normally distributed.

## Check assumption of homogeneity of variances.
testCategorical(resid2, catPred = GCN2022$eDNA_capture)
testCategorical(resid2, catPred = GCN2022$Month)
## Heteroscedasticity present for both eDNA capture method and month.

## Check dataset does not contain serial auto-correlation. This can arise 
## if there are unaccounted for relationships in data set or if there is 
## confounding effect of time or space.
recal.resid3 <- recalculateResiduals(resid2, group = GCN2022$Month)
testTemporalAutocorrelation(recal.resid3, time = unique(GCN2022$Month))

group.ponds <- aggregate(GCN2022[, 4:5], list(GCN2022$Pond_number), mean)
recal.resid4 <- recalculateResiduals(resid2, group = GCN2022$Pond_number)
testSpatialAutocorrelation(recal.resid4, group.ponds$Latitude, group.ponds$Longitude)
## No temporal or spatial autocorrelation detected.

## Unlike overdispersion, underdispersion is generally not considered to be a 
## major issue as estimates are more conservative as opposed to inflated with 
## overdispersed models. GLMMs are generally robust to non-normal residuals
## provided that the model predictions are a good fit to the observed data as is
## the case here. Therefore, we will proceed with post hoc analysis.


## Post-hoc test for month.
glmm.score.month.posthoc <- emmeans(glmm.score.capture.month,
                                    specs = "Month",
                                    type = "response")
multcomp::cld(glmm.score.month.posthoc, Letters = letters)
pairs(glmm.score.month.posthoc)
plot(glmm.score.month.posthoc)

## Post-hoc test for eDNA capture method.
glmm.score.capture.posthoc <- emmeans(glmm.score.capture.month,
                                      specs = "eDNA_capture",
                                      type = "response")
multcomp::cld(glmm.score.capture.posthoc, Letters = letters)
pairs(glmm.score.capture.posthoc)
plot(glmm.score.capture.posthoc)


## Plot eDNA score for ponds each month with each eDNA capture method.
p4a <- ggplot(GCN2022,
             aes(x = eDNA_capture, y = eDNA_score, fill = eDNA_capture)) +
  geom_jitter(pch = 21, cex = 2, height = 0.1, width = 0.2, show.legend = FALSE) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5, 
               colour = "black", fill = "gold") +
  scale_colour_manual(values = c("black", "white")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(breaks = seq(0, 12, 1)) +
  labs(title = "(a)",
       x = "eDNA capture method", y = "eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Month)
p4a

## Create new dataframe with eDNA scores for filtration and ethanol precipitation
## in wide format.
method.comp <- GCN2022 %>%
  dplyr::select(Pond_ID, Pond_number, eDNA_capture, Month, eDNA_score) %>%
  spread(eDNA_capture, eDNA_score) %>%
  distinct() %>%
  mutate(capture_diff = F - EP) %>%
  drop_na(capture_diff) %>%
  droplevels()

## Plot eDNA scores with ethanol precipitation against eDNA scores with filtration
## for each pond.
p4b <- ggplot(method.comp,
            aes(x = EP, y = F, label= Pond_number)) +
  geom_abline(intercept = 0, linetype = 2, linewidth = 1, colour = "red1") +
  geom_point(cex = 2) +
  geom_text_repel(aes(label = Pond_number),
                  size = 5,
                  min.segment.length = unit(0, "lines"),
                  colour = "black",
                  segment.color = "black",
                  max.overlaps = 25) +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  labs(title = "(b)",
       x = "eDNA score (EP)", y = "eDNA score (F)") +
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text.y = element_text(angle = 360),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Month)
p4b

## Obtain predicted values from GLMM for the full data set.
scores.capture.month.fit <- data.frame(effect(c("Month", "eDNA_capture"), 
                                              glmm.score.capture.month))

## Plot predicted against observed values for GCN detection in each month.
p4c <- ggplot() + 
  geom_errorbar(scores.capture.month.fit, 
                mapping = aes(x = eDNA_capture, ymin = lower, ymax = upper), 
                linewidth = 1, width = 0.3, color = "black") + 
  geom_point(scores.capture.month.fit, 
             mapping = aes(x = eDNA_capture, y = fit, fill = eDNA_capture), 
             size = 5, shape = 22) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "(c)",
       x = "eDNA capture method", y = "Proportional eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Month)
p4c


## Create multi-plot:
g3 <- ggarrange(p4a, p4b, p4c,
                nrow = 3, ncol = 1)

## Export plot.
#ggsave(filename = "Figures/eDNA_score_month_capture.pdf", 
#       plot = g3, width = 20, height = 20, dpi = 300, units = "in")



##################################
# eDNA CAPTURE METHOD AND SEASON #
##################################

## Subset GCN2022 dataframe for ponds with both in-season and out-of-season data.
GCN.season <- GCN2022 %>%
  filter(!Pond_number %in% c("17", "18", "19", "21", "23", "24", "25")) %>%
  droplevels() 

## Subset further for in-season and limited out-of-season (July, August) data.
GCN.season.ltd <- GCN.season %>%
  filter(!Month %in% c("September", "October")) %>%
  droplevels()


#---------------------------------#
# INCLUDING SEPTEMBER AND OCTOBER #
#---------------------------------#

#============================#
# NUMBER OF POSITIVE SAMPLES #
#============================#

## Test whether the difference in the number of positive samples in-season and 
## out-of-season using ethanol precipitation and filtration is statistically 
## significant using a binomial GLM.

## With interaction term between eDNA capture method and season:
glmm.pa.capture.season <- glmmTMB(GCN_binary ~ (1|Pond_number) + eDNA_capture 
                                  + Season + eDNA_capture:Season,
                                  family = binomial,
                                  data = GCN.season)
summary(glmm.pa.capture.season)
drop1(glmm.pa.capture.season, test = "Chi")


## Interaction term was not significant so run GLMM without:
glmm.pa.capture.season <- glmmTMB(GCN_binary ~ (1|Pond_number) + eDNA_capture 
                                  + Season,
                                  family = binomial,
                                  data = GCN.season)
summary(glmm.pa.capture.season)
drop1(glmm.pa.capture.season, test = "Chi")

## Run model diagnostics.
diagnose(glmm.pa.capture.season)                     # No issues detected by glmmTMB
check_model(glmm.pa.capture.season)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.pa.capture.season)           # Response follows Bernoulli distribution
check_singularity(glmm.pa.capture.season)            # Model fit is not singular
check_convergence(glmm.pa.capture.season)            # Model has converged
model_performance(glmm.pa.capture.season)            # Model has low explanatory power
testDispersion(glmm.pa.capture.season, plot = T)     # Model is not under or overdispersed
testZeroInflation(glmm.pa.capture.season, plot = T)  # Data is not zero-inflated

## Get standardised residuals for model validation.
resid3 <- simulateResiduals(glmm.pa.capture.season, plot = T)
getResiduals(glmm.pa.capture.season)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid3, plot = F)
testQuantiles(resid3, plot = F)
testOutliers(resid3)
## Normally distributed.

## Check assumption of homogeneity of variances.
testCategorical(resid3, catPred = GCN.season$eDNA_capture)
testCategorical(resid3, catPred = GCN.season$Month)
## Heteroscedascity present for month.

## Check dataset does not contain serial auto-correlation. This can arise 
## if there are unaccounted for relationships in data set or if there is 
## confounding effect of time or space.
recal.resid5 <- recalculateResiduals(resid3, group = GCN.season$Month)
testTemporalAutocorrelation(recal.resid5, time = unique(GCN.season$Month))

group.ponds <- aggregate(GCN.season[, 4:5], list(GCN.season$Pond_number), mean)
recal.resid6 <- recalculateResiduals(resid3, group = GCN.season$Pond_number)
testSpatialAutocorrelation(recal.resid6, group.ponds$Latitude, group.ponds$Longitude)
## No temporal autocorrelation but spatial autocorrelation detected.


## Post-hoc test for season.
glmm.pa.season.posthoc <- emmeans(glmm.pa.capture.season, 
                                  specs = "Season",
                                  type = "response")
multcomp::cld(glmm.pa.season.posthoc, Letters = letters)
pairs(glmm.pa.season.posthoc)
plot(glmm.pa.season.posthoc)

## Post-hoc test for eDNA capture method.
glmm.pa.season.posthoc <- emmeans(glmm.pa.capture.season, 
                                  specs = "eDNA_capture",
                                  type = "response")
multcomp::cld(glmm.pa.season.posthoc, Letters = letters)
pairs(glmm.pa.season.posthoc)
plot(glmm.pa.season.posthoc)


## Calculate the number of GCN positive and negative samples in-season and 
## out-of-season using ethanol precipitation and filtration.
capture.samples.season <- GCN.season %>%
  group_by(eDNA_capture, Season) %>%
  count(GCN_2022) %>%
  rename(Samples = n) %>%
  complete(GCN_2022)

## Plot number of GCN positive and negative samples in-season and out-of-season
## with each capture method.
p5a <- ggplot(capture.samples.season, 
             aes(x = eDNA_capture, y = Samples, fill = GCN_2022)) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_text(aes(label = Samples), 
            position = position_dodge(width = 0.9), 
            vjust = -1) +
  scale_fill_manual(name = "Great crested newt",
                    values = c("#262D38", "orange")) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  labs(title = "(a)",
       subtitle = "(i)",
       x = "eDNA capture method", y = "Number of eDNA samples") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "top",
        text = element_text(size = 20)) +
  facet_grid(. ~ Season)
p5a

## Obtain predicted values from GLMM for the full data set.
pa.capture.season.fit <- data.frame(effect(c("Season", "eDNA_capture"), 
                                           glmm.pa.capture.season))

## Plot predicted against observed values for GCN detection in each season.
p5b <- ggplot() + 
  geom_errorbar(pa.capture.season.fit, 
                mapping = aes(x = eDNA_capture, ymin = lower, ymax = upper), 
                linewidth = 1, width = 0.3, color = "black") + 
  geom_point(pa.capture.season.fit, 
             mapping = aes(x = eDNA_capture, y = fit, fill = eDNA_capture), 
             size = 5, shape = 22) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "(b)",
       subtitle = "(i)",
       x = "eDNA capture method", y = "Probability of a positive eDNA result") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Season)
p5b


#============#
# eDNA SCORE #
#============#

## Test whether the difference in proportional eDNA score for samples in-season 
## vs. out-of-season using ethanol precipitation and filtration is statistically 
## significant using a binomial GLMM.

## With interaction term between eDNA capture method and season:
glmm.score.capture.season <- glmmTMB(prop_eDNA_score ~ (1|Pond_number) 
                                     + eDNA_capture + Season + eDNA_capture:Season, 
                                     family = binomial(link = "logit"),
                                     data = GCN.season)
summary(glmm.score.capture.season)
drop1(glmm.score.capture.season, test = "Chi")


## The interaction term was not significant so run GLMM without With interaction term:
glmm.score.capture.season <- glmmTMB(prop_eDNA_score ~ (1|Pond_number) 
                                     + eDNA_capture + Season, 
                                     family = binomial(link = "logit"),
                                     data = GCN.season)
summary(glmm.score.capture.season)
drop1(glmm.score.capture.season, test = "Chi")


## Run model diagnostics.
diagnose(glmm.score.capture.season)                     # No issues detected by glmmTMB
check_model(glmm.score.capture.season)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.score.capture.season)           # Response follows binomial distribution
check_singularity(glmm.score.capture.season)            # Model fit is singular
check_convergence(glmm.score.capture.season)            # Model has converged
model_performance(glmm.score.capture.season)            # Model has low explanatory power
testDispersion(glmm.score.capture.season, plot = T)     # Model is underdispersed
testZeroInflation(glmm.score.capture.season, plot = T)  # Data is zero-inflated

## Get standardised residuals for model validation.
resid4 <- simulateResiduals(glmm.score.capture.season, plot = T)
getResiduals(glmm.score.capture.season)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid4, plot = F)
testQuantiles(resid4, plot = F)
testOutliers(resid4)
## Not normally distributed.

## Check assumption of homogeneity of variances.
testCategorical(resid4, catPred = GCN.season$eDNA_capture)
testCategorical(resid4, catPred = GCN.season$Month)
## Heteroscedascity present for eDNA capture method and month.

## Check dataset does not contain serial auto-correlation. This can arise 
## if there are unaccounted for relationships in data set or if there is 
## confounding effect of time or space.
recal.resid7 <- recalculateResiduals(resid4, group = GCN.season$Month)
testTemporalAutocorrelation(recal.resid7, time = unique(GCN.season$Month))

group.ponds <- aggregate(GCN.season[, 4:5], list(GCN.season$Pond_number), mean)
recal.resid8 <- recalculateResiduals(resid4, group = GCN.season$Pond_number)
testSpatialAutocorrelation(recal.resid8, group.ponds$Latitude, group.ponds$Longitude)
## No temporal autocorrelation but spatial autocorrelation detected.


## Post-hoc test for season.
glmm.score.season.posthoc <- emmeans(glmm.score.capture.season,
                                     specs = "Season",
                                     type = "response")
multcomp::cld(glmm.score.season.posthoc, Letters = letters)
pairs(glmm.score.season.posthoc)
plot(glmm.score.season.posthoc)

## Post-hoc test for eDNA capture method.
glmm.score.capture.posthoc <- emmeans(glmm.score.capture.season,
                                      specs = "eDNA_capture",
                                      type = "response")
multcomp::cld(glmm.score.capture.posthoc, Letters = letters)
pairs(glmm.score.capture.posthoc)
plot(glmm.score.capture.posthoc)


## Plot eDNA scores for samples in-season and out-of-season with each capture 
## method.
p6a <- ggplot(GCN.season,
             aes(x = eDNA_capture, y = eDNA_score, fill = eDNA_capture)) +
  geom_jitter(pch = 21, cex = 2, height = 0.1, width = 0.2, show.legend = FALSE) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5, 
               colour = "black", fill = "gold") +
  scale_colour_manual(values = c("black", "white")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(breaks = seq(0, 12, 1)) +
  labs(title = "(a)",
       subtitle = "(i)",
       x = "eDNA capture method", y = "eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Season)
p6a

## Obtain predicted values from GLMM for the full data set.
scores.capture.season.fit <- data.frame(effect(c("Season", "eDNA_capture"), 
                                               glmm.score.capture.season))

## Plot predicted against observed values for eDNA score in each season.
p6b <- ggplot() + 
  geom_errorbar(scores.capture.season.fit, 
                mapping = aes(x = eDNA_capture, ymin = lower, ymax = upper), 
                linewidth = 1, width = 0.3, color = "black") + 
  geom_point(scores.capture.season.fit, 
             mapping = aes(x = eDNA_capture, y = fit, fill = eDNA_capture), 
             size = 5, shape = 22) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "(b)",
       subtitle = "(i)",
       x = "eDNA capture method", y = "Proportional eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Season)
p6b


#---------------------------------#
# EXCLUDING SEPTEMBER AND OCTOBER #
#---------------------------------#

#============================#
# NUMBER OF POSITIVE SAMPLES #
#============================#

## Test whether the difference in the number of positive samples in-season and 
## out-of-season using ethanol precipitation and filtration is statistically 
## significant using a binomial GLM.

## With interaction term between eDNA capture method and season:
glmm.pa.capture.season.ltd <- glmmTMB(GCN_binary ~ (1|Pond_number) 
                                      + eDNA_capture + Season + eDNA_capture:Season,
                                      family = binomial,
                                      data = GCN.season.ltd)
summary(glmm.pa.capture.season.ltd)
drop1(glmm.pa.capture.season.ltd, test = "Chi")


## Interaction term was not significant to run GLM without:
glmm.pa.capture.season.ltd <- glmmTMB(GCN_binary ~ (1|Pond_number) 
                                      + eDNA_capture + Season,
                                      family = binomial,
                                      data = GCN.season.ltd)
summary(glmm.pa.capture.season.ltd)
drop1(glmm.pa.capture.season.ltd, test = "Chi")

## Run model diagnostics.
diagnose(glmm.pa.capture.season.ltd)                     # No issues detected by glmmTMB
check_model(glmm.pa.capture.season.ltd)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.pa.capture.season.ltd)           # Response follows Bernoulli distribution
check_singularity(glmm.pa.capture.season.ltd)            # Model fit is not singular
check_convergence(glmm.pa.capture.season.ltd)            # Model has converged
model_performance(glmm.pa.capture.season.ltd)            # Model has low explanatory power
testDispersion(glmm.pa.capture.season.ltd, plot = T)     # Model is not under or overdispersed
testZeroInflation(glmm.pa.capture.season.ltd, plot = T)  # Data is not zero-inflated

## Get standardised residuals for model validation.
resid5 <- simulateResiduals(glmm.pa.capture.season.ltd, plot = T)
getResiduals(glmm.pa.capture.season.ltd)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid5, plot = F)
testQuantiles(resid5, plot = F)
testOutliers(resid5)
## Normally distributed.

## Check assumption of homogeneity of variances.
testCategorical(resid5, catPred = GCN.season.ltd$eDNA_capture)
testCategorical(resid5, catPred = GCN.season.ltd$Month)
## No heteroscedascity present.

## Check dataset does not contain serial auto-correlation. This can arise 
## if there are unaccounted for relationships in data set or if there is 
## confounding effect of time or space.
recal.resid9 <- recalculateResiduals(resid5, group = GCN.season.ltd$Month)
testTemporalAutocorrelation(recal.resid9, time = unique(GCN.season.ltd$Month))

group.ponds <- aggregate(GCN.season.ltd[, 4:5], list(GCN.season.ltd$Pond_number), mean)
recal.resid10 <- recalculateResiduals(resid5, group = GCN.season.ltd$Pond_number)
testSpatialAutocorrelation(recal.resid10, group.ponds$Latitude, group.ponds$Longitude)
## No temporal autocorrelation but spatial autocorrelation detected.


## Post-hoc test for season.
glmm.pa.season.ltd.posthoc <- emmeans(glmm.pa.capture.season.ltd, 
                                      specs = "Season",
                                      type = "response")
multcomp::cld(glmm.pa.season.ltd.posthoc, Letters = letters)
pairs(glmm.pa.season.ltd.posthoc)
plot(glmm.pa.season.ltd.posthoc)

## Post-hoc test for eDNA capture method.
glmm.pa.season.ltd.posthoc <- emmeans(glmm.pa.capture.season.ltd, 
                                      specs = "eDNA_capture",
                                      type = "response")
multcomp::cld(glmm.pa.season.ltd.posthoc, Letters = letters)
pairs(glmm.pa.season.ltd.posthoc)
plot(glmm.pa.season.ltd.posthoc)


## Calculate the number of GCN positive and negative samples in-season and 
## out-of-season using ethanol precipitation and filtration.
capture.samples.season.ltd <- GCN.season.ltd %>%
  group_by(eDNA_capture, Season) %>%
  count(GCN_2022) %>%
  rename(Samples = n) %>%
  complete(GCN_2022)

## Plot number of GCN positive and negative samples in-season and out-of-season
## with each capture method.
p7a <- ggplot(capture.samples.season.ltd, 
              aes(x = eDNA_capture, y = Samples, fill = GCN_2022)) + 
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  geom_text(aes(label = Samples), 
            position = position_dodge(width = 0.9), 
            vjust = -1) +
  scale_fill_manual(name = "Great crested newt",
                    values = c("#262D38", "orange")) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  labs(title = "",
       subtitle = "(ii)",
       x = "eDNA capture method", y = "Number of eDNA samples") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "top",
        text = element_text(size = 20)) +
  facet_grid(. ~ Season)
p7a

## Obtain predicted values from GLMM for the full data set.
pa.capture.season.ltd.fit <- data.frame(effect(c("Season", "eDNA_capture"), 
                                               glmm.pa.capture.season.ltd))

## Plot predicted against observed values for GCN detection in each season.
p7b <- ggplot() + 
  geom_errorbar(pa.capture.season.ltd.fit, 
                mapping = aes(x = eDNA_capture, ymin = lower, ymax = upper), 
                linewidth = 1, width = 0.3, color = "black") + 
  geom_point(pa.capture.season.ltd.fit, 
             mapping = aes(x = eDNA_capture, y = fit, fill = eDNA_capture), 
             size = 5, shape = 22) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "",
       subtitle = "(ii)",
       x = "eDNA capture method", y = "Probability of a positive eDNA result") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Season)
p7b


#============#
# eDNA SCORE #
#============#

## Test whether the difference in proportional eDNA score for samples in-season 
## vs. out-of-season using ethanol precipitation and filtration is statistically 
## significant using a binomial GLMM.

## With interaction term between eDNA capture method and season:
glmm.score.capture.season.ltd <- glmmTMB(prop_eDNA_score ~ (1|Pond_number) 
                                         + eDNA_capture + Season + eDNA_capture:Season, 
                                         family = binomial(link = "logit"),
                                         data = GCN.season.ltd)
summary(glmm.score.capture.season.ltd)
drop1(glmm.score.capture.season.ltd, test = "Chi")


## The interaction term not significant so run GLMM without With interaction term:
glmm.score.capture.season.ltd <- glmmTMB(prop_eDNA_score ~ (1|Pond_number) 
                                         + eDNA_capture + Season, 
                                         family = binomial(link = "logit"),
                                         data = GCN.season.ltd)
summary(glmm.score.capture.season.ltd)
drop1(glmm.score.capture.season.ltd, test = "Chi")

## Run model diagnostics.
diagnose(glmm.score.capture.season.ltd)                     # No issues detected by glmmTMB
check_model(glmm.score.capture.season.ltd)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.score.capture.season.ltd)           # Response follows binomial distribution
check_singularity(glmm.score.capture.season.ltd)            # Model fit is not singular
check_convergence(glmm.score.capture.season.ltd)            # Model has converged
model_performance(glmm.score.capture.season.ltd)            # Model has low explanatory power
testDispersion(glmm.score.capture.season.ltd, plot = T)     # Model is underdispersed
testZeroInflation(glmm.score.capture.season.ltd, plot = T)  # Data is zero-inflated

## Get standardised residuals for model validation.
resid6 <- simulateResiduals(glmm.score.capture.season.ltd, plot = T)
getResiduals(glmm.score.capture.season.ltd)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid6, plot = F)
testQuantiles(resid6, plot = F)
testOutliers(resid6)
## Not normally distributed.

## Check assumption of homogeneity of variances.
testCategorical(resid6, catPred = GCN.season.ltd$eDNA_capture)
testCategorical(resid6, catPred = GCN.season.ltd$Season)
## Heteroscedascity present for eDNA capture method.

## Check dataset does not contain serial auto-correlation. This can arise 
## if there are unaccounted for relationships in data set or if there is 
## confounding effect of time or space.
recal.resid11 <- recalculateResiduals(resid6, group = GCN.season.ltd$Season)
testTemporalAutocorrelation(recal.resid11, time = unique(GCN.season.ltd$Season))

group.ponds <- aggregate(GCN.season.ltd[, 4:5], list(GCN.season.ltd$Pond_number), mean)
recal.resid12 <- recalculateResiduals(resid6, group = GCN.season.ltd$Pond_number)
testSpatialAutocorrelation(recal.resid12, group.ponds$Latitude, group.ponds$Longitude)
## No temporal autocorrelation but spatial autocorrelation detected.


## Post-hoc test for season.
glmm.score.season.ltd.posthoc <- emmeans(glmm.score.capture.season.ltd,
                                         specs = "Season",
                                         type = "response")
multcomp::cld(glmm.score.season.ltd.posthoc, Letters = letters)
pairs(glmm.score.season.ltd.posthoc)
plot(glmm.score.season.ltd.posthoc)

## Post-hoc test for eDNA capture method.
glmm.score.capture.ltd.posthoc <- emmeans(glmm.score.capture.season.ltd,
                                          specs = "eDNA_capture",
                                          type = "response")
multcomp::cld(glmm.score.capture.ltd.posthoc, Letters = letters)
pairs(glmm.score.capture.ltd.posthoc)
plot(glmm.score.capture.ltd.posthoc)


## Plot eDNA scores for samples in-season and out-of-season with each capture 
## method.
p8a <- ggplot(GCN.season.ltd,
              aes(x = eDNA_capture, y = eDNA_score, fill = eDNA_capture)) +
  geom_jitter(pch = 21, cex = 2, height = 0.1, width = 0.2, show.legend = FALSE) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) + 
  stat_summary(fun = mean, geom = "point", shape = 23, size = 5, 
               colour = "black", fill = "gold") +
  scale_colour_manual(values = c("black", "white")) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(breaks = seq(0, 12, 1)) +
  labs(title = "",
       subtitle = "(ii)",
       x = "eDNA capture method", y = "eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Season)
p8a

## Obtain predicted values from GLMM for the full data set.
scores.capture.season.ltd.fit <- data.frame(effect(c("Season", "eDNA_capture"), 
                                                   glmm.score.capture.season.ltd))

## Plot predicted against observed values for eDNA score in each season.
p8b <- ggplot() + 
  geom_errorbar(scores.capture.season.ltd.fit, 
                mapping = aes(x = eDNA_capture, ymin = lower, ymax = upper), 
                linewidth = 1, width = 0.3, color = "black") + 
  geom_point(scores.capture.season.ltd.fit, 
             mapping = aes(x = eDNA_capture, y = fit, fill = eDNA_capture), 
             size = 5, shape = 22) +
  scale_fill_manual(values = c("black", "white")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "",
       subtitle = "(ii)",
       x = "eDNA capture method", y = "Proportional eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.position = "none",
        text = element_text(size = 20)) +
  facet_grid(. ~ Season)
p8b


#---------------#
# SUMMARY PLOTS #
#---------------#

## Number of positive samples in-season vs. out-of-season:
g4 <- ggarrange(ggarrange(p5a, p7a, nrow = 1, ncol = 2, 
                          common.legend = TRUE, legend = "bottom"),
                ggarrange(p5b, p7b, nrow = 1, ncol = 2),
                nrow = 2, ncol = 1)

#ggsave(filename = "Figures/GCN_detection_rate_season_capture.pdf", 
#       plot = g4, width = 13, height = 13, dpi = 300, units = "in")


## eDNA score for samples in-season vs. out-of-season:
g5 <- ggarrange(p6a, p8a, p6b, p8b,
                nrow = 2, ncol = 2)

#ggsave(filename = "Figures/eDNA_score_season_capture.pdf", 
#       plot = g5, width = 13, height = 13, dpi = 300, units = "in")



############################
# VOLUME OF WATER FILTERED #
############################

## Create dataframe of samples processed using filtration.
filtration <- GCN2022 %>%
  filter(eDNA_capture == "F") %>%
  drop_na(Volume) %>%
  droplevels()

## Test whether the difference in eDNA score for samples of different volumes
## is statistically significant using a negative binomial GLM as a Poisson GLM 
## was overdispersed.
glmm.volume <- glmmTMB(prop_eDNA_score ~ (1|Pond_number) + Volume, 
                       family = binomial(link = "logit"),
                       data = filtration)
summary(glmm.volume)
drop1(glmm.volume, test = "Chi")

## Run model diagnostics.
diagnose(glmm.volume)                     # Large coefficients
check_model(glmm.volume)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.volume)           # Response follows binomial distribution
check_singularity(glmm.volume)            # Model fit is singular
check_convergence(glmm.volume)            # Model has converged
model_performance(glmm.volume)            # Model has low explanatory power
testDispersion(glmm.volume, plot = T)     # Model is not underdispersed
testZeroInflation(glmm.volume, plot = T)  # Data is not zero-inflated

## Get standardised residuals for model validation.
resid7 <- simulateResiduals(glmm.volume, plot = T)
getResiduals(glmm.volume)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid7, plot = F)
testQuantiles(resid7, plot = F)
testOutliers(resid7)
## Not normally distributed.

## Check dataset does not contain serial auto-correlation. This can arise 
## if there are unaccounted for relationships in data set or if there is 
## confounding effect of time or space.
recal.resid13 <- recalculateResiduals(resid7, group = filtration$Month)
testTemporalAutocorrelation(recal.resid13, time = unique(filtration$Month))

group.ponds <- aggregate(filtration[, 4:5], list(filtration$Pond_number), mean)
recal.resid14 <- recalculateResiduals(resid7, group = filtration$Pond_number)
testSpatialAutocorrelation(recal.resid14, group.ponds$Latitude, group.ponds$Longitude)
## No spatial autocorrelation but temporal autocorrelation detected.


## Obtain predicted values from GLMM for the full data set.
volume.fit <- data.frame(effect(c("Volume"), glmm.volume))

## Plot volume filtered against eDNA score for ponds.
p9 <- ggplot() +
  geom_jitter(filtration, mapping = aes(x = Volume, y = prop_eDNA_score, fill = GCN_2022),
              pch = 21, cex = 3.5, colour = "black") +
  geom_line(volume.fit, mapping = aes (x = Volume, y= fit),
            linewidth = 1) +
  geom_ribbon(volume.fit, mapping = aes(x = Volume, ymin = lower, ymax = upper),
              alpha = 0.25) +
  scale_fill_manual(name = "Great crested newt",
                    values = c("black", "orange")) +
  scale_x_continuous(limits = c(0, 2600), 
                     breaks = seq(0, 2500, 500)) +
  scale_y_continuous(limits = c(-0.05, 1.05), 
                     breaks = seq(0, 1, 0.1)) +
  labs(x = "Volume filtered (mL)", y = "Proportional eDNA score") +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", linewidth = 0.5, linetype = "solid"),
        axis.line.y = element_line(colour = "black", linewidth = 0.5, linetype = "solid"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        legend.key = element_rect(fill = NA),
        text = element_text(size = 20))
p9

## Export plot.
#ggsave(filename = "Figures/volume_vs_eDNA_score.pdf", 
#       plot = p9, width = 10, height = 8, dpi = 300, units = "in")



#########################
# POPULATION SIZE CLASS #
#########################

## No GCN were detected out-of-season using conventional methods, thus only 
## in-season population and eDNA data were compared.

## Import population size class and associated eDNA score data.
pop.dat <- read.csv("Data/Population_size_data.csv", header = TRUE)

## Refine dataset:
## - Filter dataframe for in-season samples only (as no GCN were found with
##   conventional methods out-of-season so no population data was available).
## - Remove in-season samples with no associated population data (missing or
##   no population assessment undertaken as pond was obmitted from eDNA surveys 
##   beyond May).
## - Reorder factor levels for plotting.
## - Convert Pond_number column to factor.
## - Remove Notes column.
pop.dat <- pop.dat %>%
  filter(Season == "In-season") %>%
  filter(!Kit_type == "Both") %>%
  drop_na(GCN_found) %>%
  mutate(Pond_number = as.factor(Pond_number),
         Population_class = fct_relevel(Population_class, 
                                        "None", "Small", "Medium")) %>%
  dplyr::select(-Notes) %>%
  droplevels()


#------------------#
# DATA EXPLORATION #
#------------------#

## Examine relationship between peak adult counts and eDNA score corresponding
## to the month most survey visits for population size estimation were undertaken
## (May).
p10a <- ggplot(pop.dat,
               aes(x = Peak_count, y = Corresponding_eDNA_score)) + 
  geom_jitter(pch = 16, cex = 3.5, colour = "black", alpha = 0.7) +
  scale_x_continuous(limits = c(-0.1, 26), 
                     breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-0.1, 12.5), 
                     breaks = seq(0, 12, 2)) +
  labs(title = "(a)", 
       subtitle = "(i)",
       x = "Peak adult count", y = "Corresponding eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text.y = element_text(angle = 360),
        legend.key = element_rect(fill = NA),
        text = element_text(size = 20)) +
  facet_grid(. ~ Kit_type)
p10a

## Examine relationship between population class and eDNA score corresponding
## to the month most survey visits for population size estimation were undertaken
## (May).
p10b <- ggplot(pop.dat,
               aes(x = Population_class, y = Corresponding_eDNA_score)) + 
  geom_jitter(width = 0.3, pch = 16, cex = 3.5, colour = "black", alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_y_continuous(limits = c(-0.1, 12.5), 
                     breaks = seq(0, 12, 2)) +
  labs(title = "", 
       subtitle = "(ii)",
       x = "Population size class", y = "Corresponding eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text.y = element_text(angle = 360),
        legend.key = element_rect(fill = NA),
        text = element_text(size = 20)) +
  facet_grid(. ~ Kit_type)
p10b

## Examine relationship between peak adult counts and average eDNA score
## derived from all eDNA scores for in-season months.
p10c <- ggplot(pop.dat,
               aes(x = Peak_count, y = Average_eDNA_score)) + 
  geom_jitter(pch = 16, cex = 3.5, colour = "black", alpha = 0.7) +
  scale_x_continuous(limits = c(-0.1, 26), 
                     breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-0.1, 12.5), 
                     breaks = seq(0, 12, 2)) +
  labs(title = "(b)", 
       subtitle = "(i)",
       x = "Peak adult count", y = "Average eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text.y = element_text(angle = 360),
        legend.key = element_rect(fill = NA),
        text = element_text(size = 20)) +
  facet_grid(. ~ Kit_type)
p10c

## Examine relationship between population class and and average eDNA score
## derived from all eDNA scores for in-season months.
p10d <- ggplot(pop.dat,
               aes(x = Population_class, y = Average_eDNA_score)) + 
  geom_jitter(width = 0.3, pch = 16, cex = 3.5, colour = "black", alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_y_continuous(limits = c(-0.1, 12.5), 
                     breaks = seq(0, 12, 2)) +
  labs(title = "", 
       subtitle = "(ii)",
       x = "Population size class", y = "Average eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text.y = element_text(angle = 360),
        legend.key = element_rect(fill = NA),
        text = element_text(size = 20)) +
  facet_grid(. ~ Kit_type)
p10d

## Examine relationship between peak adult counts and peak eDNA score
## derived from all eDNA scores for in-season months.
p10e <- ggplot(pop.dat,
               aes(x = Peak_count, y = Peak_eDNA_score)) + 
  geom_jitter(pch = 16, cex = 3.5, colour = "black", alpha = 0.7) +
  scale_x_continuous(limits = c(-0.1, 26), 
                     breaks = seq(0, 25, 5)) +
  scale_y_continuous(limits = c(-0.1, 12.5), 
                     breaks = seq(0, 12, 2)) +
  labs(title = "(c)", 
       subtitle = "(i)",
       x = "Peak adult count", y = "Peak eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text.y = element_text(angle = 360),
        legend.key = element_rect(fill = NA),
        text = element_text(size = 20)) +
  facet_grid(. ~ Kit_type)
p10e

## Examine relationship between population class and and peak eDNA score
## derived from all eDNA scores for in-season months.
p10f <- ggplot(pop.dat,
               aes(x = Population_class, y = Peak_eDNA_score)) + 
  geom_jitter(width = 0.3, pch = 16, cex = 3.5, colour = "black", alpha = 0.7) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  scale_y_continuous(limits = c(-0.1, 12.5), 
                     breaks = seq(0, 12, 2)) +
  labs(title = "",
       subtitle = "(ii)",
       x = "Population size class", y = "Peak eDNA score") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text.y = element_text(angle = 360),
        legend.key = element_rect(fill = NA),
        text = element_text(size = 20)) +
  facet_grid(. ~ Kit_type)
p10f

## Create multi-plot:
g6 <- ggarrange(p10a, p10b, p10c, p10d, p10e, p10f,
                nrow = 3, ncol = 2)

#ggsave(filename = "Figures/adult_counts_vs_eDNA_score.pdf", 
#       plot = g6, width = 15, height = 20, dpi = 300, units = "in")


## Calculate proportional eDNA scores for modelling and remove NAs.
pop.May <- pop.dat %>%
  drop_na(Corresponding_eDNA_score) %>%
  mutate(Corresponding_eDNA_score = round(Corresponding_eDNA_score/12, 3),
         Average_eDNA_score = round(Average_eDNA_score/12, 3),
         Peak_eDNA_score = round(Peak_eDNA_score/12, 3)) %>%
  droplevels()

## Subset data for ethanol precipitation and filtration.
pop.EP <- pop.May %>%
  filter(Kit_type == "EP") %>%
  droplevels()

pop.F <- pop.May %>%
  filter(Kit_type == "F") %>%
  droplevels()


#--------------------------#
# CORRESPONDING eDNA SCORE #
#--------------------------#

## Run GLMM for ethanol precipitation:
glmm.corr.EP <- glmmTMB(Corresponding_eDNA_score ~ (1|Pond_number) 
                        + Peak_count + Population_class,
                        family = binomial(link = "logit"),
                        data = pop.EP)
summary(glmm.corr.EP)
drop1(glmm.corr.EP, test = "Chi")

## Run model diagnostics.
diagnose(glmm.corr.EP)                     # Large coefficients
check_model(glmm.corr.EP)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.corr.EP)           # Response follows binomial distribution
check_singularity(glmm.corr.EP)            # Model fit is singular
check_convergence(glmm.corr.EP)            # Model has converged
model_performance(glmm.corr.EP)            # Model has low explanatory power
testDispersion(glmm.corr.EP, plot = T)     # Model is underdispersed
testZeroInflation(glmm.corr.EP, plot = T)  # Less 0s than expected

## Get standardised residuals for model validation.
resid8 <- simulateResiduals(glmm.corr.EP, plot = T)
getResiduals(glmm.corr.EP)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid8, plot = F)
testQuantiles(resid8, plot = F)
testOutliers(resid8)
## Normally distributed.


## Run GLMM for filtration:
glmm.corr.F <- glmmTMB(Corresponding_eDNA_score ~ (1|Pond_number)
                       + Peak_count + Population_class,
                       family = binomial(link = "logit"),
                       data = pop.F)
summary(glmm.corr.F)
drop1(glmm.corr.F, test = "Chi")

## Run model diagnostics.
diagnose(glmm.corr.F)                     # Large coefficients and convergence problems
check_model(glmm.corr.F)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.corr.F)           # Response follows binomial distribution
check_singularity(glmm.corr.F)            # Model fit is singular
check_convergence(glmm.corr.F)            # Model has converged
model_performance(glmm.corr.F)            # Model has low explanatory power
testDispersion(glmm.corr.F, plot = T)     # No problems with dispersion
testZeroInflation(glmm.corr.F, plot = T)  # Less 0s than expected

## Get standardised residuals for model validation.
resid9 <- simulateResiduals(glmm.corr.F, plot = T)
getResiduals(glmm.corr.F)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid9, plot = F)
testQuantiles(resid9, plot = F)
testOutliers(resid9)
## Some deviation from normality.


#--------------------#
# AVERAGE eDNA SCORE #
#--------------------#

## Run GLMM for ethanol precipitation:
glmm.avg.EP <- glmmTMB(Average_eDNA_score ~ (1|Pond_number) 
                       + Peak_count + Population_class,
                       family = binomial(link = "logit"),
                       data = pop.EP)
summary(glmm.avg.EP)
drop1(glmm.avg.EP, test = "Chi")

## Run model diagnostics.
diagnose(glmm.avg.EP)                     # Large coefficients
check_model(glmm.avg.EP)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.avg.EP)           # Response follows binomial distribution
check_singularity(glmm.avg.EP)            # Model fit is singular
check_convergence(glmm.avg.EP)            # Model has converged
model_performance(glmm.avg.EP)            # Model has low explanatory power
testDispersion(glmm.avg.EP, plot = T)     # Model is underdispersed
testZeroInflation(glmm.avg.EP, plot = T)  # Less 0s than expected

## Get standardised residuals for model validation.
resid10 <- simulateResiduals(glmm.avg.EP, plot = T)
getResiduals(glmm.avg.EP)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid10, plot = F)
testQuantiles(resid10, plot = F)
testOutliers(resid10)
## Not normally distributed.


## Run GLMM for filtration:
glmm.avg.F <- glmmTMB(Average_eDNA_score ~ (1|Pond_number)
                      + Peak_count + Population_class,
                      family = binomial(link = "logit"),
                      data = pop.F)
summary(glmm.avg.F)
drop1(glmm.avg.F, test = "Chi")

## Run model diagnostics.
diagnose(glmm.avg.F)                     # Large coefficient
check_model(glmm.avg.F)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.avg.F)           # Response follows binomial distribution
check_singularity(glmm.avg.F)            # Model fit is singular
check_convergence(glmm.avg.F)            # Model has converged
model_performance(glmm.avg.F)            # Model has low explanatory power
testDispersion(glmm.avg.F, plot = T)     # Data is underdispersed
testZeroInflation(glmm.avg.F, plot = T)  # Less 0s than expected

## Get standardised residuals for model validation.
resid11 <- simulateResiduals(glmm.avg.F, plot = T)
getResiduals(glmm.avg.F)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid11, plot = F)
testQuantiles(resid11, plot = F)
testOutliers(resid11)
## Some deviation from normality.


#-----------------#
# PEAK eDNA SCORE #
#-----------------#

## Run GLMM for ethanol precipitation:
glmm.peak.EP <- glmmTMB(Peak_eDNA_score ~ (1|Pond_number) 
                        + Peak_count + Population_class,
                        family = binomial(link = "logit"),
                        data = pop.EP)
summary(glmm.peak.EP)
drop1(glmm.peak.EP, test = "Chi")

## Run model diagnostics.
diagnose(glmm.peak.EP)                     # Large coefficients
check_model(glmm.peak.EP)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.peak.EP)           # Response follows binomial distribution
check_singularity(glmm.peak.EP)            # Model fit is singular
check_convergence(glmm.peak.EP)            # Model has converged
model_performance(glmm.peak.EP)            # Model has low explanatory power
testDispersion(glmm.peak.EP, plot = T)     # No issues with dispersion
testZeroInflation(glmm.peak.EP, plot = T)  # Less 0s than expected

## Get standardised residuals for model validation.
resid12 <- simulateResiduals(glmm.peak.EP, plot = T)
getResiduals(glmm.peak.EP)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid12, plot = F)
testQuantiles(resid12, plot = F)
testOutliers(resid12)
## Not normally distributed.


## Run GLMM for filtration:
glmm.peak.F <- glmmTMB(Peak_eDNA_score ~ (1|Pond_number)
                       + Peak_count + Population_class,
                       family = binomial(link = "logit"),
                       data = pop.F)
summary(glmm.peak.F)
drop1(glmm.peak.F, test = "Chi")

## Run model diagnostics.
diagnose(glmm.peak.F)                     # Large coefficient
check_model(glmm.peak.F)                  # Potentially problems with residuals but otherwise ok
check_distribution(glmm.peak.F)           # Response follows binomial distribution
check_singularity(glmm.peak.F)            # Model fit is singular
check_convergence(glmm.peak.F)            # Model has converged
model_performance(glmm.peak.F)            # Model has low explanatory power
testDispersion(glmm.peak.F, plot = T)     # Data is underdispersed
testZeroInflation(glmm.peak.F, plot = T)  # Less 0s than expected

## Get standardised residuals for model validation.
resid11 <- simulateResiduals(glmm.peak.F, plot = T)
getResiduals(glmm.peak.F)

## Check assumption of normal distribution. Most models are robust to slight 
## deviations from normality in the residuals.
testUniformity(resid11, plot = F)
testQuantiles(resid11, plot = F)
testOutliers(resid11)
## Some deviation from normality.

## Neither peak adult counts or population size class are related to corresponding
## eDNA score, average eDNA score, or peak eDNA score for each pond, thus 
## population size does not appear to influence the amount of GCN eDNA present.



#######################
# OCCUPANCY MODELLING #
#######################

## A potential flaw with GLMMs is they do not account for imperfect detection 
## via a hierarchical dependence structure like multi-scale occupancy models do.
## We will try applying eDNA occupancy models to our data, acknowledging that 
## repeat visits to the same pond will be treated as unique visits to different 
## ponds. 


#-----------------------------#
# GRIFFIN ET AL. (2020) MODEL #
#-----------------------------#

## The Griffin et al. (2020) model (https://doi.org/10.1111/rssc.12390), 
## implemented via the eDNAShinyApp package (Diana et al. 2021; 
## https://doi.org/10.1111/ecog.05718), accounts for false positive and 
## false negative error during sampling (stage 1) and qPCR (stage 2). However, 
## this package does not allow for different covariates to specified at each level 
## of the model. Both eDNA capture method and month are likely to influence 
## probability of eDNA being captured in a sample, but only eDNA capture method 
## is likely to influence probability of qPCR amplification as the DNA extraction 
## procedure differs for ethanol precipitation versus filtration.

## We will run the model using eDNA capture method and month as covariates at
## Stage 1 only.

## Format data for upload to eDNA Shiny App.
GCN.app <- GCN2022 %>%
  dplyr::select(eDNA_score, eDNA_capture, Month) %>%
  drop_na()
write.csv(GCN.app, "Data/GCN_data_for_eDNAShinyApp.csv", row.names = FALSE)

## Install and load eDNAShinyApp package.
#library(devtools)
#install_github("alexdiana1992/eDNAShinyApp", ref = "master")
library(eDNAShinyApp)
runeDNAShinyApp()

## All Stage 1 and Stage 2 observations were considered with default prior 
## settings from Griffin et al. (2020) of 0.9 for θ11 and p11 and 0.1 for θ10 
## and p10. The model was run using 2,000 burn-in iterations, 10,000 iterations, 
## 1 chain and 20 thinned iterations, with ‘probability of site occupancy’ set to 
## 0.5, ‘variance of probability of site occupancy’ set to 4, ‘variance of 
## coefficients of probability of site occupancy’ set to 0.25, and number of 
## significant covariates’ set to 2.

## Read in relevant results from the eDNAShinyApp.
app.theta11 <- read.csv("eDNAshinyapp/output/download_theta11.csv", header = TRUE) %>%
  mutate(Error = "True Positive")
app.theta10 <- read.csv("eDNAshinyapp/output/download_theta10.csv", header = TRUE) %>%
  mutate(Error = "False Positive")
app.beta.theta11 <- read.csv("eDNAshinyapp/output/download_beta_theta11.csv", header = TRUE) %>%
  mutate(Error = "True Positive")
app.beta.theta10 <- read.csv("eDNAshinyapp/output/download_beta_theta10.csv", header = TRUE) %>%
  mutate(Error = "False Positive")

## Combine dataframes for plotting.
app.theta <- rbind(app.theta11, app.theta10) %>%
  rename(lower = X2.5Credible.Interval, upper = X97.5Credible.Interval) %>%
  group_by(Error) %>%
  summarise(mean = mean(Posterior.Mean),
            median = median(Posterior.Mean),
            lower = min(lower),
            upper = max(upper)) %>%
  ungroup()

app.beta.theta <- rbind(app.beta.theta11, app.beta.theta10) %>%
  mutate(X = fct_relevel(X,
                         "(Intercept)", "eDNA_captureF", "MonthMay", "MonthJune",
                         "MonthJuly", "MonthAugust", "MonthSeptember", "MonthOctober"))

## Plot results.
p11a <- ggplot(app.theta, aes(x = Error, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  scale_x_discrete(labels = c("False Positive" = expression(theta[10]),
                              "True Positive" = expression(theta[11]))) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "(a) Griffin et al. (2020) model",
       subtitle = "(i)",
       x = "", 
       y = expression(paste(theta, " (median +/- 95% credible intervals)"))) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        text = element_text(size = 20))
p11a

app.beta.theta$Error <- factor(app.beta.theta$Error, 
                               levels = c("False Positive", "True Positive"),
                               ordered = TRUE, 
                               labels = c(expression(theta[10]), expression(theta[11])))

p11b <- ggplot(app.beta.theta, aes(x = X, y = Posterior.Mean))+
  geom_point() +
  geom_errorbar(aes(ymin = X2.5Credible.Interval, ymax = X97.5Credible.Interval), 
                width = 0) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "black", linewidth = 1) +
  scale_y_continuous(limits = c(-6, 5), breaks = seq(-6, 5, 1)) +
  labs(title = "",
       subtitle = "(ii)",
       x = "", 
       y = "\nCoefficient") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(colour = "black"),
        text = element_text(size = 20)) +
  facet_grid(. ~ Error, labeller = label_parsed)
p11b


#---------------------------------#
# DORAZIO & ERICKSON (2018) MODEL #
#---------------------------------#

## The Dorazio & Erickson (2018) model (https://doi.org/10.1111/1755-0998.12735)
## only accounts for false negative error but does allow different covariates to 
## specified at each level. This Bayesian model is computationally demanding. The 
## msocc package (Stratton et al. 2020; https://doi.org/10.1111/2041-210X.13442)
## provides an alternative method, implementing a Gibbs sampler that quickly 
## converges to the posterior distribution and features a web application capable 
## of investigating study design, fitting models and analysing the results.

## Install and load msocc package.
#devtools::install_github("StrattonCh/msocc")
library(msocc)

## Code modified from Stelzer et al. (2023). See https://doi.org/10.1002/edn3.501
## and https://github.com/jamesord/PKD_eDNA).
## Define custom functions. Modified lppd function to include burn-in option:
compute_lppd_burnin <- function(msocc_mod, burnin){
  # pull model info
  num.mcmc <- msocc_mod$model.info$num.mcmc-burnin
  J <- msocc_mod$model.info$J
  K <- msocc_mod$model.info$K
  z <- tail(msocc_mod$model.info$z.vec,num.mcmc)[, rep(1:ncol(msocc_mod$model.info$z.vec), K)]
  A <- tail(msocc_mod$model.info$A,num.mcmc)[, rep(1:ncol(msocc_mod$model.info$A), K)]
  y <- matrix(rep(msocc_mod$model.info$y, num.mcmc), nrow = num.mcmc, byrow = T)
  # convert MCMC samples to probability
  psi <- tail(psi_mcmc(msocc_mod), num.mcmc)
  theta <- tail(theta_mcmc(msocc_mod), num.mcmc)
  p <- tail(p_mcmc(msocc_mod), num.mcmc)
  site <- psi^z * (1 - psi)^(1 - z)
  sample <- (z*theta)^A * (1 - z*theta)^(1 - A)
  rep <- (A*p)^y * (1 - A*p)^(1 - y)
  tmp <- (site * sample * rep)[2:num.mcmc,]
  return(sum(log(colMeans(tmp))))
}

## Modified pwaic2 function to include burn-in option:
compute_pwaic2_burnin <- function(msocc_mod, burnin){
  # pull model info
  num.mcmc <- msocc_mod$model.info$num.mcmc-burnin
  J <- msocc_mod$model.info$J
  K <- msocc_mod$model.info$K
  z <- tail(msocc_mod$model.info$z.vec,num.mcmc)[, rep(1:ncol(msocc_mod$model.info$z.vec), K)]
  A <- tail(msocc_mod$model.info$A,num.mcmc)[, rep(1:ncol(msocc_mod$model.info$A), K)]
  y <- matrix(rep(msocc_mod$model.info$y, num.mcmc), nrow = num.mcmc, byrow = T)
  # first piece of the difference
  # convert MCMC samples to probability
  psi <- tail(psi_mcmc(msocc_mod), num.mcmc)
  theta <- tail(theta_mcmc(msocc_mod), num.mcmc)
  p <- tail(p_mcmc(msocc_mod), num.mcmc)
  site <- psi^z * (1 - psi)^(1 - z)
  sample <- (z*theta)^A * (1 - z*theta)^(1 - A)
  rep <- (A*p)^y * (1 - A*p)^(1 - y)
  tmp <- log(site * sample * rep)[2:num.mcmc,]
  pwaic2 <- sum(apply(tmp, 2, var))
  return(pwaic2)
}

## Modified waic function (with type 2 penalty score) to include burn-in option:
waic2_burnin <- function(msocc_mod, burnin){
  waic <- -2 * (compute_lppd_burnin(msocc_mod,burnin) - compute_pwaic2_burnin(msocc_mod,burnin))
  return(waic)
}

## Function to get a table of WAIC and associated parameters from model fit:
get_waic_table <- function(msocc_mod, burnin){
  formula <- paste("psi(",msocc_mod$model.info$site_mod[2],"),
                   theta(",msocc_mod$model.info$samp_mod[2],"),
                   p(",msocc_mod$model.info$rep_mod[2],")",sep="")
  waic <- waic2_burnin(msocc_mod, burnin)
  pwaic2 <- compute_pwaic2_burnin(msocc_mod, burnin)
  lppd <- compute_lppd_burnin(msocc_mod, burnin)
  num.mcmc <- msocc_mod$model.info$num.mcmc
  return(data.frame(formula = formula,
                    waic = waic,
                    waic2 = pwaic2,
                    lppd = lppd,
                    iter = num.mcmc,
                    burnin = burnin))
}

## Format data for occupancy modelling.
GCN.msocc <- GCN2022 %>%
  dplyr::select(Pond_number, Month, eDNA_capture, eDNA_score) %>%
  unite("sample", eDNA_capture:Month, remove = FALSE) %>%
  mutate(sample = as.factor(fct_relevel(sample,
                                        "EP_April", "F_April",
                                        "EP_May", "F_May",
                                        "EP_June", "F_June",
                                        "EP_July", "F_July",
                                        "EP_August", "F_August",
                                        "EP_September", "F_September",
                                        "EP_October", "F_October"))) %>%
  mutate(sample = as.numeric(sample)) %>%
  dplyr::select(-c(eDNA_capture, Month)) %>%
  mutate(pcr1 = ifelse(eDNA_score > 0, 1, 0),
         pcr2 = ifelse(eDNA_score > 1, 1, 0),
         pcr3 = ifelse(eDNA_score > 2, 1, 0),
         pcr4 = ifelse(eDNA_score > 3, 1, 0),
         pcr5 = ifelse(eDNA_score > 4, 1, 0),
         pcr6 = ifelse(eDNA_score > 5, 1, 0),
         pcr7 = ifelse(eDNA_score > 6, 1, 0),
         pcr8 = ifelse(eDNA_score > 7, 1, 0),
         pcr9 = ifelse(eDNA_score > 8, 1, 0),
         pcr10 = ifelse(eDNA_score > 9, 1, 0),
         pcr11 = ifelse(eDNA_score > 10, 1, 0),
         pcr12 = ifelse(eDNA_score > 11, 1, 0)) %>%
  dplyr::select(-eDNA_score) %>%
  rename(site = Pond_number)

## Create covariate tables for occupancy modelling.
GCN.msocc.site <- GCN.msocc %>%
  dplyr::select(-c(sample:pcr12)) %>%
  distinct(site, .keep_all = TRUE)

GCN.msocc.sample <- GCN2022 %>%
  dplyr::select(Pond_number, Month, eDNA_capture) %>%
  unite("sample", eDNA_capture:Month, remove = FALSE) %>%
  mutate(sample = as.factor(fct_relevel(sample,
                                        "EP_April", "F_April",
                                        "EP_May", "F_May",
                                        "EP_June", "F_June",
                                        "EP_July", "F_July",
                                        "EP_August", "F_August",
                                        "EP_September", "F_September",
                                        "EP_October", "F_October"))) %>%
  mutate(sample = as.numeric(sample)) %>%
  rename(site = Pond_number)


## Fit null occupancy model to data.
GCN.msocc.null <- msocc_mod(GCN.msocc, 
                            num.mcmc = 11000,
                            site = list(model = ~1, cov_tbl = GCN.msocc.site),
                            sample = list(model = ~1, cov_tbl = GCN.msocc.sample),
                            rep = list(model = ~1, cov_tbl = GCN.msocc.sample),
                            progress = T, 
                            print = NULL, 
                            seed = NULL, 
                            beta_bin = T)

## Obtain a one-row dataframe with the model formula, WAIC score and associated 
## metrics (used to compare all subsequent models).
null_score <- get_waic_table(GCN.msocc.null, burnin = 1000)
null_score

## If no covariates are specified in the model, the posterior is already 
## converted to the probability scale. For the trace plots, we want to see the 
## sampling distribution on the original logit scale, so use the logit() function.
par(mfrow = c(2,2))
plot(logit(GCN.msocc.null$psi), type= "l", main = "beta", ylab = "sampled value")
plot(logit(GCN.msocc.null$theta), type= "l", main = "alpha", ylab = "sampled value")
plot(logit(GCN.msocc.null$p), type= "l", main = "delta", ylab = "sampled value")

## Fit occupancy model to data with covariates at the sample and replicate levels.
## Month is likely to influence probability of eDNA capture because great crested
## newt abundance and activity will vary over time. eDNA capture method is likely
## to influence the probably of eDNA capture, where filtration is expected to
## capture more DNA present in a water sample. Only eDNA capture method is
## expected to influence probability of eDNA amplification as the extraction
## methods for ethanol precipitation and filtration differ.

## Model with no covariates at sample level:
GCN.msocc.m1 <- msocc_mod(GCN.msocc, 
                          num.mcmc = 11000,
                          site = list(model = ~1, cov_tbl = GCN.msocc.site),
                          sample = list(model = ~1, cov_tbl = GCN.msocc.sample),
                          rep = list(model = ~ eDNA_capture, cov_tbl = GCN.msocc.sample),
                          progress = T, 
                          print = NULL, 
                          seed = NULL, 
                          beta_bin = T)
m1_score <- get_waic_table(GCN.msocc.m1, burnin = 1000)

## Model with only month at sample level:
GCN.msocc.m2 <- msocc_mod(GCN.msocc, 
                          num.mcmc = 11000,
                          site = list(model = ~1, cov_tbl = GCN.msocc.site),
                          sample = list(model = ~ Month, cov_tbl = GCN.msocc.sample),
                          rep = list(model = ~1, cov_tbl = GCN.msocc.sample),
                          progress = T, 
                          print = NULL, 
                          seed = NULL, 
                          beta_bin = T)
m2_score <- get_waic_table(GCN.msocc.m2, burnin = 1000)

## Model with only eDNA_capture method at sample level:
GCN.msocc.m3 <- msocc_mod(GCN.msocc, 
                          num.mcmc = 11000,
                          site = list(model = ~1, cov_tbl = GCN.msocc.site),
                          sample = list(model = ~ eDNA_capture, cov_tbl = GCN.msocc.sample),
                          rep = list(model = ~1, cov_tbl = GCN.msocc.sample),
                          progress = T, 
                          print = NULL, 
                          seed = NULL, 
                          beta_bin = T)
m3_score <- get_waic_table(GCN.msocc.m3, burnin = 1000)

## Model with month and eDNA_capture method at sample level:
GCN.msocc.m4 <- msocc_mod(GCN.msocc, 
                          num.mcmc = 11000,
                          site = list(model = ~1, cov_tbl = GCN.msocc.site),
                          sample = list(model = ~ Month + eDNA_capture, cov_tbl = GCN.msocc.sample),
                          rep = list(model = ~1, cov_tbl = GCN.msocc.sample),
                          progress = T, 
                          print = NULL, 
                          seed = NULL, 
                          beta_bin = T)
m4_score <- get_waic_table(GCN.msocc.m4, burnin = 1000)

## Model with eDNA capture method at both levels:
GCN.msocc.m5 <- msocc_mod(GCN.msocc, 
                          num.mcmc = 11000,
                          site = list(model = ~1, cov_tbl = GCN.msocc.site),
                          sample = list(model = ~ eDNA_capture, cov_tbl = GCN.msocc.sample),
                          rep = list(model = ~ eDNA_capture, cov_tbl = GCN.msocc.sample),
                          progress = T, 
                          print = NULL, 
                          seed = NULL, 
                          beta_bin = T)
m5_score <- get_waic_table(GCN.msocc.m5, burnin = 1000)

## Model with eDNA capture method at replicate level and month at sample level:
GCN.msocc.m6 <- msocc_mod(GCN.msocc, 
                          num.mcmc = 11000,
                          site = list(model = ~1, cov_tbl = GCN.msocc.site),
                          sample = list(model = ~ Month, cov_tbl = GCN.msocc.sample),
                          rep = list(model = ~ eDNA_capture, cov_tbl = GCN.msocc.sample),
                          progress = T, 
                          print = NULL, 
                          seed = NULL, 
                          beta_bin = T)
m6_score <- get_waic_table(GCN.msocc.m6, burnin = 1000)

## Model with all covariates at relevant levels:
GCN.msocc.m7 <- msocc_mod(GCN.msocc, 
                          num.mcmc = 11000,
                          site = list(model = ~1, cov_tbl = GCN.msocc.site),
                          sample = list(model = ~ Month + eDNA_capture, cov_tbl = GCN.msocc.sample),
                          rep = list(model = ~ eDNA_capture, cov_tbl = GCN.msocc.sample),
                          progress = T, 
                          print = NULL, 
                          seed = NULL, 
                          beta_bin = T)
m7_score <- get_waic_table(GCN.msocc.m7, burnin = 1000)

## Model with all covariates at relevant levels and interaction term:
GCN.msocc.m8 <- msocc_mod(GCN.msocc, 
                          num.mcmc = 11000,
                          site = list(model = ~1, cov_tbl = GCN.msocc.site),
                          sample = list(model = ~ Month*eDNA_capture, cov_tbl = GCN.msocc.sample),
                          rep = list(model = ~ eDNA_capture, cov_tbl = GCN.msocc.sample),
                          progress = T, 
                          print = NULL, 
                          seed = NULL, 
                          beta_bin = T)
m8_score <- get_waic_table(GCN.msocc.m8, burnin = 1000)

## Use WAIC scores to compare the models. The WAIC is an 'information criterion' 
## reflecting a compromise between how well the model fits the data and how 
## generalisable it is. Too far in either direction leads to underfitting or 
## overfitting (underfitting = mdoel does not fit data very well; overfitting =
## model fits data too well and is less likely to reflect reality). Lower WAIC 
## scores reflect a better compromise between model fit and 'generalisability', 
## and therefore models with lower scores are considered 'better'.

## Compile all model formulae and WAIC scores into a dataframe.
model_scores <- rbind(null_score, m1_score, m2_score, m3_score, m4_score, m5_score, 
                      m6_score, m7_score, m8_score)
model_scores <- model_scores[order(model_scores$waic),] # sort by WAIC score in ascending order (lower = better)

## Models with scores of 1-2 unit difference of the 'best' model deserve 
## consideration, while those with 3-7 units of difference have less support.
## We could use this rule in the other direction to say that models with scores 
## at least 3 units lower than the null model deserve consideration.
## Either way, the model with all covariates at relevant levels, but no interaction
## term, has the lowest WAIC score.

## Overall summary of occupancy given by posterior medians:
head(posterior_summary(GCN.msocc.m7, level = "overall", burnin = 1000))

## For a more in-depth look at each of the site, sample, and replicate levels 
## and a description of uncertainty, we can specify the level in posterior_summary.
## Posterior summary of occupancy at the sample: 
head(posterior_summary(GCN.msocc.m7, level = "sample"))

## Summarise the sample level occurrence probabilities by plotting credibility 
## intervals. n controls how many samples are plotted on each plot.
cred.plots <- cred_plot(GCN.msocc.m7, level = "sample", n = 14)
gridExtra::grid.arrange(cred.plots[[1]], cred.plots[[2]], nrow = 2, ncol = 1)

cred.plots <- cred_plot(GCN.msocc.m7, level = "rep", n = 40)
cred.plots

## Examine the trace plots of the parameter estimates for best model before going 
## further with estimates of theta and p.
par(mfrow = c(2,3))
plot(logit(GCN.msocc.m7$psi), type = "l", main = "beta", ylab = "sampled value")
plot(as.data.frame(GCN.msocc.m7$alpha)$`(Intercept)`, 
     type = "l", main = "alpha.intercept", ylab = "sampled value")
plot(as.data.frame(GCN.msocc.m7$alpha)$MonthMay, 
     type = "l", main = "alpha.month", ylab = "sampled value")
plot(as.data.frame(GCN.msocc.m7$alpha)$eDNA_capture, 
     type = "l", main = "alpha.eDNAcapture", ylab = "sampled value")
plot(as.data.frame(GCN.msocc.m7$delta)$`(Intercept)`,
     type = "l", main = "delta.intercept", ylab = "sampled value")
plot(as.data.frame(GCN.msocc.m7$delta)$eDNA_capture,
     type = "l", main = "delta.eDNAcapture", ylab = "sampled value")

## Extract posterior samples of delta parameters from the model (last 10k samples
## = 1k burn-in), and convert to p estimates for ethanol precipitation and filtration.
alphas <- tail(as.data.frame(GCN.msocc.m7$alpha), 10000)
deltas <- tail(as.data.frame(GCN.msocc.m7$delta), 10000)

thetas <- data.frame(theta = c(inv.logit(alphas$`(Intercept)`),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthMay),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthJune),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthJuly),
                               inv.logit(alphas$`(Intercept)`+ alphas$MonthAugust),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthSeptember),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthOctober),
                               inv.logit(alphas$`(Intercept)` + alphas$eDNA_captureF),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthMay + alphas$eDNA_captureF),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthJune + alphas$eDNA_captureF),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthJuly + alphas$eDNA_captureF),
                               inv.logit(alphas$`(Intercept)`+ alphas$MonthAugust + alphas$eDNA_captureF),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthSeptember + alphas$eDNA_captureF),
                               inv.logit(alphas$`(Intercept)` + alphas$MonthOctober + alphas$eDNA_captureF)),
                     Month = c(rep("April", 10000),
                               rep("May", 10000),
                               rep("June", 10000),
                               rep("July", 10000),
                               rep("August", 10000),
                               rep("September", 10000),
                               rep("October", 10000),
                               rep("April", 10000),
                               rep("May", 10000),
                               rep("June", 10000),
                               rep("July", 10000),
                               rep("August", 10000),
                               rep("September", 10000),
                               rep("October", 10000)),
                     eDNA_capture = c(rep("EP", 70000), rep("F", 70000)))

ps <- data.frame(p = c(inv.logit(deltas$`(Intercept)`),
                       inv.logit(deltas$`(Intercept)` + deltas$eDNA_captureF)),
                 eDNA_capture = c(rep("EP", 10000), rep("F", 10000)))

## Obtain median and credible intervals for theta.
theta_summarised <- thetas %>%
  group_by(Month, eDNA_capture) %>%
  summarise(median = median(theta),
            lower = quantile(theta, probs = c(0.025, 0.975))[1],
            upper = quantile(theta, probs = c(0.025, 0.975))[2]) %>%
  ungroup() %>%
  mutate(Month = fct_relevel(Month,
                             "April", "May", "June", "July", "August",
                             "September", "October"))

## Obtain median and credible intervals for p.
p_summarised <- ps %>% 
  group_by(eDNA_capture) %>%
  summarise(median = median(p),
            lower = quantile(p, probs = c(0.025,0.975))[1],
            upper = quantile(p, probs = c(0.025,0.975))[2])

## Plot of theta as a function of month and eDNA capture method.
p12a <- ggplot(theta_summarised, aes(x = eDNA_capture, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "\n(b) Dorazio & Erickson (2018) model",
       subtitle = "(i)",
       x = "eDNA capture method", 
       y = expression(paste(theta, " (median +/- 95% credible intervals)"))) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        text = element_text(size = 20)) +
  facet_grid(. ~ Month)
p12a

## Plot of p as a function of eDNA capture method.
p12b <- ggplot(p_summarised, aes(x = eDNA_capture, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "\n",
       subtitle = "(ii)",
       x = "eDNA capture method", 
       y = "\np (median +/- 95% credible intervals)") +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        text = element_text(size = 20))
p12b


#-----------------------------#
# NICHOLS ET AL. (2008) MODEL #
#-----------------------------#

## The eDNA capture methods used in this study are co-located, i.e. the not truly 
## independent of each other as would be assumed by a single season occupancy model
## for each method. Additionally, there are three possible reasons why great crested
## newt eDNA may not be detected:
##
## - eDNA is truly absent.
## - eDNA is present but simply not captured by either eDNA capture method.
## - eDNA is present but simply “unavailable” to be detected by any of the methods. 
##   For instance, if newts do not inhabit a pond, it will be impossible for any 
##   eDNA release and subsequent capture to occur.
##
## Single season, multi-method occupancy models allow estimation of occupancy, 
## local availability, and detection rates for different monitoring methods, thus
## providing method-specific detection probabilities (Nichols et al. 2008; 
## https://doi.org/10.1111/j.1365-2664.2008.01509.x).

## Fit single season, multi-method occupancy model with RPresence (MacKenzie & 
## Hines 2025) to compare great crested newt detection probability with filtration 
## and ethanol precipitation.

## Install and load RPresence package.
#install.packages("RPresence", repo = 'https://eesc.usgs.gov/mbr/mbrCRAN')
library(RPresence)

## Format data for occupancy modelling.
GCN.occ <- GCN2022 %>%
  dplyr::select(Pond_number, eDNA_capture, Month, GCN_binary) %>%
  unite("Survey", eDNA_capture:Month, remove = FALSE) %>%
  dplyr::select(-c(eDNA_capture, Month)) %>%
  pivot_wider(names_from = Survey, values_from = GCN_binary) %>%
  column_to_rownames("Pond_number")

## Create RPresence input object using data.
pao <- createPao(data = GCN.occ, 
                 frq = rep(1, nrow(GCN.occ)),
                 methods = 2)
summary(pao)

## Single season, multi-method occupancy model. Fix site occupancy (psi) and local
## occupancy to 1 because newts are known to occur across all three sites where 
## ponds were sampled:
mmmod <- occMod(data = pao,      
                type = 'so.mm',                        # so.mm : static occupancy, multi-method model
                model = list(psi ~ 1,                  # occupancy constant for all sites
                             theta ~ 1,                # local occupancy constant for each survey (month)
                             p ~ METHOD * SURVEY))     # Detection different for each method and survey (month) 


## Model summary:
summary(mmmod)

## Print real estimates. Detection parameter, p, is printed with the method 
## number just after ‘p’, followed by the survey number in parentheses. The 
## survey number is referenced by season and survey-within-season, separated by 
## a ‘-’.
print(mmmod$real)
coef(object = mmmod, param = "p", prob = 0.95)
mmmod$dmat

## Make dataframe with results from RPresence.
occ.dat <- rbind(mmmod$real$psi,
                 mmmod$real$p,
                 mmmod$real$theta) %>%
  mutate(parameter = c("psi", "April F", "April EP", "May F", "May EP", "June F",
                       "June EP", "July F", "July EP", "August F", "August EP",
                       "September F", "September EP", "October F", "October EP",
                       "theta April", "theta May", "theta June", "theta July",
                       "theta August", "theta September", "theta October"))

## Plot all model estimates.
ggplot(data = occ.dat, aes(x = parameter, y = est)) +
  geom_col() +
  ylab("Estimate") +
  geom_errorbar(aes(ymin = ciLow, ymax = ciHi), width = 0.25) +
  ylim(c(0, 1))

## Make data frame with detection probability estimates.
det.dat <- mmmod$real$p %>% 
  mutate(method = c("F", "EP", "F", "EP", "F", "EP", "F", "EP", "F", "EP", "F", 
                    "EP", "F", "EP"),
         month = c("April", "April", "May", "May", "June", "June", "July", "July", 
                   "August", "August", "September", "September", "October", "October"),
         month = fct_relevel(month,
                             "April", "May", "June", "July", "August", "September", 
                             "October"))

## Plot detection probability estimates.
p13 <- ggplot(det.dat, aes(x = method, y = est, fill = method)) +
  geom_bar(stat = "identity", colour = "black") +
  geom_errorbar(aes(ymin = ciLow, ymax = ciHi), width = 0.2, position = "dodge") + 
  scale_fill_manual(values = c("grey60", "white")) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(title = "\n(c) Nichols et al. (2008) model",
       subtitle = "",
       x = "eDNA capture method", y = "Probability of eDNA detection") + 
  theme_bw() + 
  theme(panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"), 
        plot.title = element_text(face = "bold"),
        axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        text = element_text(size = 20),
        legend.position = "none") +
  facet_grid(. ~ month)
p13

## The problem with this type of occupancy model is it doesn't account for the
## hierarchical nature of eDNA analysis, i.e. probability that a species is
## present at a site, probability that a species' DNA is captured in a sample,
## and probability that a species' DNA is amplified in a qPCR replicate. However,
## models that do account for this are not equipped to handle multi-method, 
## repeated sampling of the same sites.


#---------#
# SUMMARY #
#---------#

## Create multi-plot with results of all occupancy models:
g7 <- (p11a | p11b) / (p12a + p12b + plot_layout(widths = c(2, 1), guides = "collect")) / p13

## Export plot.
#ggsave(filename = "Figures/occupancy_models.pdf", 
#       plot = g7, width = 15, height = 20, dpi = 300, units = "in")


#################
# END OF SCRIPT #
#################
