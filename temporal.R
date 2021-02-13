# Calculations for temporal niche manuscript
# 1. Data prep: Conversion to relative times (sunrise = 0600, sunset = 1800)
# 2. Data prep: Determine crepuscular captures
# 3. Temporal comparison: w/in species by year
# 4. Temporal comparison: btw species all years combined (main test)
# 5. Temporal comparison: single vs. groups of raccoons
# 6. Temporal comparison: single sp. vs. overlap
# 7. Distribution of capture times: diurnal, crepuscular, nocturnal
# 8. Rose diagrams for raccoons at sites with >10 opossum detections
# 9. Create an activity overlap plot


# LOAD LIBRARIES
library("circular")
library("ggplot2")
library("dplyr")
library("forcats")
library("camtrapR")

# CONVERT HH:MM:SS SINCE SUNRISE / SUNSET TO DECIMAL HOURS
# Upload data
captures <- read.csv("capturetimes.csv", header = T)

# Converts the timestamp data to decimal hours
dec_since_start <- sapply(strsplit(as.character(captures$since_start),":"),
       function(x) {
         x <- as.numeric(x)
         x[1]+x[2]/60+x[3]/3600
       }
)

# Append to data frame
captures <- cbind.data.frame(captures,dec_since_start)


# CONVERT TO RELATIVE (DECIMAL HOURS & CAPTURE HOUR) TIME ON 24H CLOCK
# Convert absolute number of hours since start to relative number of hours since start of window
rel_since_start <- captures$dec_since_start / captures$window_length * 12
captures <- cbind.data.frame(captures,rel_since_start)

# Convert to a relative time on 24h clock
rel_time <- numeric()
for (i in 1:length(captures$rel_since_start)) {
  if (captures$window_type[i] == "day") {
    rel_time <- c(rel_time,captures$rel_since_start[i] + 6)
  } else {
    new_time = captures$rel_since_start[i] + 18
    if (new_time >= 24) {
      rel_time <- c(rel_time, new_time - 24)
    } else {
      rel_time <- c(rel_time, new_time)
    }
  }
}
captures <- cbind.data.frame(captures,rel_time)

# Truncate to get capture hour
rel_hour <- floor(captures$rel_time)
captures <- cbind.data.frame(captures,rel_hour)


# DETERMINE WHICH CAPTURES WERE CREPUSCULAR
levels(captures$window_type) <- c("day", "night", "crepuscular")
for (i in 1:length(captures$rel_since_start)) {
  if (captures$rel_since_start[i] < 1.2 || captures$rel_since_start[i] > 10.8) {
    captures$window_type[i] <- as.factor("crepuscular")
  }
}


# SUBSET BY SPECIES AND BY YEAR WITHIN SPECIES
pl <- subset(captures, species == "Procyon lotor")
pl14 <- subset(pl, project == 'SC14')
pl15 <- subset(pl, project == 'SC15')
pl16 <- subset(pl, project == 'SC16')
pl17 <- subset(pl, project == 'SC17')

dv <- subset(captures, species == "Didelphis virginiana")
dv14 <- subset(dv, project == 'SC14')
dv15 <- subset(dv, project == 'SC15')
dv16 <- subset(dv, project == 'SC16')
dv17 <- subset(dv, project == 'SC17')


# CONVERT TO CIRCULAR VARIABLES ON RELATIVE TIME (DECIMAL HOUR)
# Also compute summary statistics for each
pl.circular <- circular(pl$rel_time, units = "hours", template = "clock24")
mean(pl.circular)
rho.circular(pl.circular)  # Measure of concentration (p. 602 in Zar)
angular.variance(pl.circular) # p. 604 in Zar

pl14.circular <- circular(pl14$rel_time, units = "hours", template = "clock24")
mean(pl14.circular) + 24 # b/c negative
rho.circular(pl14.circular)
angular.variance(pl14.circular)

pl15.circular <- circular(pl15$rel_time, units = "hours", template = "clock24")
mean(pl15.circular)
rho.circular(pl15.circular)
angular.variance(pl15.circular)

pl16.circular <- circular(pl16$rel_time, units = "hours", template = "clock24")
mean(pl16.circular) + 24
rho.circular(pl16.circular)
angular.variance(pl16.circular)

pl17.circular <- circular(pl17$rel_time, units = "hours", template = "clock24")
mean(pl17.circular)
rho.circular(pl17.circular)
angular.variance(pl17.circular)


dv.circular <- circular(dv$rel_time, units = "hours", template = "clock24")
mean(dv.circular) + 24
rho.circular(dv.circular)
angular.variance(dv.circular)

dv14.circular <- circular(dv14$rel_time, units = "hours", template = "clock24")
mean(dv14.circular) + 24
rho.circular(dv14.circular)
angular.variance(dv14.circular)

dv15.circular <- circular(dv15$rel_time, units = "hours", template = "clock24")
mean(dv15.circular) + 24
rho.circular(dv15.circular)
angular.variance(dv15.circular)

dv16.circular <- circular(dv16$rel_time, units = "hours", template = "clock24")
mean(dv16.circular) + 24
rho.circular(dv16.circular)
angular.variance(dv16.circular)

dv17.circular <- circular(dv17$rel_time, units = "hours", template = "clock24")
mean(dv17.circular) + 24
rho.circular(dv17.circular)
angular.variance(dv17.circular)

# Convert means to timestamp format
# Raccoon mean = 00:09:54
mean(pl.circular) # 0.1648872  0 hours
0.1648872*60 # 9.893232 minutes
0.893232*60 # 53.59392 seconds

# Opossum mean = 23:26:11
mean(dv.circular)+24 # 23.4363
0.4363*60 # 26.178
0.178*60 # 10.68
# Opossum 95% CI (in circular format): -2.272647 to 1.270454
mean(dv.circular)+1.96*sd(dv.circular)
mean(dv.circular)-1.96*sd(dv.circular)


# TEST DIFFERENCES WITHIN SPECIES BETWEEN YEARS
# Determine if any dataset follows a uniform circular distribution using the Rayleigh test
rayleigh.test(pl.circular)
rayleigh.test(pl14.circular)
rayleigh.test(pl15.circular)
rayleigh.test(pl16.circular)
rayleigh.test(pl17.circular)

rayleigh.test(dv.circular)
rayleigh.test(dv14.circular)
rayleigh.test(dv15.circular)
rayleigh.test(dv16.circular)
rayleigh.test(dv17.circular)

# Because they do not follow a uniform distribution, can compare among groups
# Use Watson-Wheeler (non-parametric)
watson.wheeler.test(list(pl14.circular, pl15.circular,pl16.circular,pl17.circular))
watson.wheeler.test(list(dv14.circular, dv15.circular,dv16.circular,dv17.circular))

# Post-hoc comparison between pairs of years for raccoons.
watson.wheeler.test(list(pl14.circular, pl15.circular))
watson.wheeler.test(list(pl14.circular, pl16.circular))
watson.wheeler.test(list(pl14.circular, pl17.circular))
watson.wheeler.test(list(pl15.circular,pl16.circular))
watson.wheeler.test(list(pl15.circular,pl17.circular))
watson.wheeler.test(list(pl16.circular,pl17.circular))


# TEST DIFFERENCES BETWEEN SPECIES COMBINING ALL YEARS
watson.two.test(pl.circular, dv.circular)


# Rose diagrams, separate, color
# png(filename = "analysis/figures/pl_circular.png", width = 5, height = 5, units = "in", res = 600, pointsize = 12)
rose.diag(pl.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 2.2, col = "darkolivegreen3")
arrows.circular(mean(pl.circular), lwd = 2)
# dev.off()

# png(filename = "analysis/figures/dv_circular.png", width = 5, height = 5, units = "in", res = 600, pointsize = 12)
rose.diag(dv.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 2.3, col = "darkolivegreen3")
arrows.circular(mean(dv.circular), lwd = 2)
# dev.off()

# Rose diagrams, separate, B&W
#png(filename = "analysis/figures/pl_circular_bw.png", width = 5, height = 5, units = "in", res = 600, pointsize = 12)
rose.diag(pl.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 2.2, col = "grey", shrink = 0.95, tol = 0)
arrows.circular(mean(pl.circular), lwd = 2)
# dev.off()

# png(filename = "analysis/figures/dv_circular_bw.png", width = 5, height = 5, units = "in", res = 600, pointsize = 12)
rose.diag(dv.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 2.3, col = "grey", shrink = 0.95, tol = 0)
arrows.circular(mean(dv.circular), lwd = 2)
# dev.off()


# COMPARE SINGLE VS GROUPS OF RACCOONS
pl_single <- subset(pl, animals == 1)
pl_single.circular <- circular(pl_single$rel_time, units = "hours", template = "clock24")

pl_group <- subset(pl, animals > 1)
pl_group.circular <- circular(pl_group$rel_time, units = "hours", template = "clock24")

rayleigh.test(pl_single.circular)
rayleigh.test(pl_group.circular)

watson.two.test(pl_single.circular,pl_group.circular)

# Rose diagrams of the comparison
# png(filename = "analysis/figures/pl_single_group.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(pl_single.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 2, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Single raccoons")
rose.diag(pl_group.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 2, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Groups")
par(op)
# dev.off()



# COMPARE OVERLAPPING VS. NON-OVERLAPPING LOCATIONS
# Sample sizes in each category
table(subset(captures, project == 'SC14')$overlap)
table(subset(captures, project == 'SC15')$overlap)
table(subset(captures, project == 'SC16')$overlap)
table(subset(captures, project == 'SC17')$overlap)

# Make data frames
plboth14 <- subset(pl14, overlap == 'B')
plonly14 <- subset(pl14, overlap == 'P')
dvboth14 <- subset(dv14, overlap == 'B')
dvonly14 <- subset(dv14, overlap == 'D')

plboth14.circular <- circular(plboth14$rel_time, units = "hours", template = "clock24")
plonly14.circular <- circular(plonly14$rel_time, units = "hours", template = "clock24")
dvboth14.circular <- circular(dvboth14$rel_time, units = "hours", template = "clock24")
dvonly14.circular <- circular(dvonly14$rel_time, units = "hours", template = "clock24")

plboth15 <- subset(pl15, overlap == 'B')
plonly15 <- subset(pl15, overlap == 'P')
dvboth15 <- subset(dv15, overlap == 'B')
dvonly15 <- subset(dv15, overlap == 'D')

plboth15.circular <- circular(plboth15$rel_time, units = "hours", template = "clock24")
plonly15.circular <- circular(plonly15$rel_time, units = "hours", template = "clock24")
dvboth15.circular <- circular(dvboth15$rel_time, units = "hours", template = "clock24")
dvonly15.circular <- circular(dvonly15$rel_time, units = "hours", template = "clock24")

plboth16 <- subset(pl16, overlap == 'B')
plonly16 <- subset(pl16, overlap == 'P')
dvboth16 <- subset(dv16, overlap == 'B')
dvonly16 <- subset(dv16, overlap == 'D')

plboth16.circular <- circular(plboth16$rel_time, units = "hours", template = "clock24")
plonly16.circular <- circular(plonly16$rel_time, units = "hours", template = "clock24")
dvboth16.circular <- circular(dvboth16$rel_time, units = "hours", template = "clock24")
dvonly16.circular <- circular(dvonly16$rel_time, units = "hours", template = "clock24")

plboth17 <- subset(pl17, overlap == 'B')
plonly17 <- subset(pl17, overlap == 'P')
dvboth17 <- subset(dv17, overlap == 'B')
dvonly17 <- subset(dv17, overlap == 'D')

plboth17.circular <- circular(plboth17$rel_time, units = "hours", template = "clock24")
plonly17.circular <- circular(plonly17$rel_time, units = "hours", template = "clock24")
dvboth17.circular <- circular(dvboth17$rel_time, units = "hours", template = "clock24")
dvonly17.circular <- circular(dvonly17$rel_time, units = "hours", template = "clock24")

rayleigh.test(plboth14.circular)
rayleigh.test(plonly14.circular)
rayleigh.test(dvboth14.circular)
rayleigh.test(dvonly14.circular)
rayleigh.test(plboth15.circular)
rayleigh.test(plonly15.circular)
rayleigh.test(dvboth15.circular)
rayleigh.test(dvonly15.circular)
rayleigh.test(plboth16.circular)
rayleigh.test(plonly16.circular)
rayleigh.test(dvboth16.circular)
rayleigh.test(dvonly16.circular)
rayleigh.test(plboth17.circular)
rayleigh.test(plonly17.circular)
rayleigh.test(dvboth17.circular)
rayleigh.test(dvonly17.circular)

watson.two.test(plboth14.circular,plonly14.circular)
watson.two.test(plboth15.circular,plonly15.circular)
watson.two.test(plboth16.circular,plonly16.circular)
watson.two.test(plboth17.circular,plonly17.circular)

watson.two.test(dvboth14.circular,dvonly14.circular)
watson.two.test(dvboth15.circular,dvonly15.circular)
watson.two.test(dvboth16.circular,dvonly16.circular)
watson.two.test(dvboth17.circular,dvonly17.circular)


# Raccoon figures showing overlapping and non-overlapping sites
# png(filename = "analysis/figures/pl14_overlap.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(plboth14.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.9, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Both spp")
rose.diag(plonly14.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.9, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Raccoon only")
# title(main = "Raccoon 2014", adj = 0, line = -10) # for screen
title(main = "Raccoon 2014", adj = 0, line = -3) # for saving to file
par(op)
# dev.off()

# png(filename = "analysis/figures/pl15_overlap.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(plboth15.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Both spp")
rose.diag(plonly15.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.9, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Raccoon only")
# title(main = "Raccoon 2015", adj = 0, line = -10) # for screen
title(main = "Raccoon 2015", adj = 0, line = -3) # for saving to file
par(op)
# dev.off()

# png(filename = "analysis/figures/pl16_overlap.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(plboth16.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.9, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Both spp")
rose.diag(plonly16.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.9, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Raccoon only")
# title(main = "Raccoon 2016", adj = 0, line = -10) # for screen
title(main = "Raccoon 2016", adj = 0, line = -3) # for saving to file
par(op)
# dev.off()

# png(filename = "analysis/figures/pl17_overlap.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(plboth17.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Both spp")
rose.diag(plonly17.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.9, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Raccoon only")
# title(main = "Raccoon 2017", adj = 0, line = -10) # for screen
title(main = "Raccoon 2017", adj = 0, line = -3) # for saving to file
par(op)
# dev.off()


# Opossum figures showing overlapping and non-overlapping sites
# png(filename = "analysis/figures/dv14_overlap.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(dvboth14.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Both spp")
rose.diag(dvonly14.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Opossum only")
# title(main = "Opossum 2014", adj = 0, line = -10) # for screen
title(main = "Opossum 2014", adj = 0, line = -3) # for saving to file
par(op)
# dev.off()

# png(filename = "analysis/figures/dv15_overlap.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(dvboth15.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Both spp")
rose.diag(dvonly15.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Opossum only")
# title(main = "Opossum 2015", adj = 0, line = -10) # for screen
title(main = "Opossum 2015", adj = 0, line = -3) # for saving to file
par(op)
# dev.off()

# png(filename = "analysis/figures/dv16_overlap.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(dvboth16.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Both spp")
rose.diag(dvonly16.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Opossum only")
# title(main = "Opossum 2016", adj = 0, line = -10) # for screen
title(main = "Opossum 2016", adj = 0, line = -3) # for saving to file
par(op)
# dev.off()

# png(filename = "analysis/figures/dv17_overlap.png", width = 8, height = 5, units = "in", res = 600, pointsize = 12)
op <- par(mfrow = c(1,2), mar = c(0,0,0,0), oma = c(0,0,0,0))
rose.diag(dvboth17.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.8, -1.2, adj = c(0,0), labels = "Both spp")
rose.diag(dvonly17.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.4, col = "darkolivegreen3")
text(-0.5, -1.2, adj = c(0,0), labels = "Opossum only")
# title(main = "Opossum 2017", adj = 0, line = -10) # for screen
title(main = "Opossum 2017", adj = 0, line = -3) # for saving to file
par(op)
# dev.off()

# Number of captures of each capture type (to determine sample sizes of above tests)
nrow(plonly14)
nrow(dvonly14)
nrow(plboth14)
nrow(dvboth14)
nrow(plonly15)
nrow(dvonly15)
nrow(plboth15)
nrow(dvboth15)
nrow(plonly16)
nrow(dvonly16)
nrow(plboth16)
nrow(dvboth16)
nrow(plonly17)
nrow(dvonly17)
nrow(plboth17)
nrow(dvboth17)


# DISTRIBUTION OF CAPTURE TIMES: DIURNAL VS. CREPUSCULAR VS. NOCTURNAL
# Set probabilities, which aren't equal because the crepuscular window is only 20% of the day
probs <- c(0.4, 0.4, 0.2)

# Get capture counts for each type
pl_distribution <- table(pl$window_type)
dv_distribution <- table(dv$window_type)

# Chi square tests
chisq.test(pl_distribution,p = probs)
chisq.test(dv_distribution,p = probs)

# Capture counts to data frame for both species
capture_df_dv <- as.data.frame(dv_distribution)
capture_df_dv <- cbind.data.frame(capture_df_dv,rep("Virginia opossum",3))
names(capture_df_dv) <- c("time.category","captures","species")

capture_df_pl <- as.data.frame(pl_distribution)
capture_df_pl <- cbind.data.frame(capture_df_pl,rep("Raccoon",3))
names(capture_df_pl) <- c("time.category","captures","species")

capture_df <- rbind.data.frame(capture_df_dv,capture_df_pl)

# Add diurnal and nocturnal as levels before substituting them
levels(capture_df$time.category) <- c(levels(capture_df$time.category),"diurnal","nocturnal")
capture_df$time.category[capture_df$time.category == "day"] <- "diurnal"
capture_df$time.category[capture_df$time.category == "night"] <- "nocturnal"

# Reorder the categories with fact_relevel from the forcats package then pass to
# grouped barplot for captures for each species in different time categories
#png(filename = "analysis/figures/fig4.png", width = 5, height = 5, units = "in", res = 600, pointsize=12)
capture_df %>%
  mutate(time.category = fct_relevel(time.category, "diurnal",
                            "crepuscular", "nocturnal")) %>%
  ggplot(aes(x=time.category, y=captures, fill = species)) + 
    geom_bar(stat="identity", position = "dodge") + 
    ylab("Number of Detections") +
    xlab("Time Category") +
    scale_fill_manual(values=c("gray42", "gray84")) +
    theme_classic() +
    theme(legend.position = "top", legend.title = element_blank()) +
    theme(legend.text = element_text(size = 12)) +
    theme(axis.text = element_text(size = 14)) +
    theme(axis.title.x = element_text(size = 16, margin = margin(10,0,0,0))) +
    theme(axis.title.y = element_text(size = 16, margin = margin(0,6,0,0)))
#dev.off()



# DISTRIBUTION OF CAPTURE TIMES: CREPUSCULAR VS. NOCTURNAL
# Set probabilities. Nocturnal + crepuscular is 14.4 hours, so the proportions are 67% and 33%
probs <- c(2/3, 1/3)

# Get capture counts for each type
pl_distribution <- table(pl$window_type)[2-3]
dv_distribution <- table(dv$window_type)[2-3]

# Chi square tests
chisq.test(pl_distribution,p = probs)
chisq.test(dv_distribution,p = probs)


# SUBSET RACCOONS AT SAMPLE UNITS WITH MANY OPOSSUM DETECTIONS AND MAKE ROSE DIAGRAMS
pl14_highdv <- subset (pl14, (studyarea == 'Colman' & unit == 2 ) | (studyarea == 'Seward Park' & unit == 6) | (studyarea == 'West Duwamish' & unit == 4))
pl14_highdv.circular <- circular(pl14_highdv$rel_time, units = "hours", template = "clock24")

pl15_highdv <- subset (pl15, (studyarea == 'Seward Park' & unit == 1 ) | (studyarea == 'Seward Park' & unit == 4) | (studyarea == 'Delridge' & unit == 5) | (studyarea == 'Delridge' & unit == 1) | (studyarea == 'Lincoln Park' & unit == 2))
pl15_highdv.circular <- circular(pl15_highdv$rel_time, units = "hours", template = "clock24")

# Rose diagrams, separate, color
# rose.diag(pl14_highdv.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.5, col = "darkolivegreen3")
# arrows.circular(mean(pl14_highdv.circular), lwd = 2)

# rose.diag(pl15_highdv.circular, axes = TRUE, bins = 24, upper = FALSE, prop = 1.3, col = "darkolivegreen3")
# arrows.circular(mean(pl15_highdv.circular), lwd = 2)

# Watson two-sample tests comparing subset to entire set within each year
watson.two.test(pl14_highdv.circular,pl14.circular)
watson.two.test(pl15_highdv.circular,pl15.circular)


#9. Create an activity overlap plot

library("camtrapR")

#Load data file with capture times
records <- read.csv("rel_time.csv", header=T)

#Ensure Date/Time column is properly formatted
records$DateTimeOriginal <- strptime(paste(records$Date, records$Time, sep = " "), format = "%d-%b-%y %H:%M") 
records$DateTimeOriginal <- as.character(records$DateTimeOriginal)

#Remove NAs
records <- na.omit(records)

#Load Stations file
stations <- read.csv("stations.csv", header=T)

#Define species A (raccoon) and B (opossum) for plotting
speciesA <- "Pl"
speciesB <- "Dv"

#Create activity overlap plot
activityOverlap(records,
                speciesA,
                speciesB,
                speciesCol = "Species",
                recordDateTimeCol = "DateTimeOriginal",
                recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                plotR = TRUE, 
                writePNG = FALSE, 
                addLegend = TRUE,
                legendPosition = "topleft",
                plotDirectory, 
                createDir = FALSE, 
                pngMaxPix = 1000,
                add.rug = TRUE,
                overlapEstimator = c("Dhat4")  # if more than 75 observations for each species, use dhat4
)


# Using the Overlap package to format the activity overlap plot

#Load packages
library(overlap)
library(dplyr)

#Load data for each species
relpl <- read.csv("rel_time_pl.csv", header=T)

reldv <- read.csv("rel_time_dv.csv", header=T)

#Convert date time to radians (radians are needed for analyses in the overlap package)
ClocktimeToTimeRad <- function(Clocktime,
                               timeformat = "%Y-%m-%d %H:%M:%S"){
  
  DateTime2 <- strptime(as.character(Clocktime), format = timeformat, tz = "UTC")
  Time2     <- format(DateTime2, format = "%H:%M:%S", usetz = FALSE)
  Time.rad  <- (as.numeric(as.POSIXct(strptime(Time2, format = "%H:%M:%S", tz = "UTC"))) -
                  as.numeric(as.POSIXct(strptime("0", format = "%S", tz = "UTC")))) / 3600 * (pi/12)
  return(Time.rad)
}

#Convert DateTimeOriginal column to radians, with the output as Values
activitytime_radians <- ClocktimeToTimeRad(records$DateTimeOriginal)

#Append radian activity time to original activity time dataset
records <- cbind(records, activitytime_radians)

#Subset radian activity times by species (raccoon)
pl_radians <- records[records$Species == "Pl", ]
pl_radians <- pl_radians %>% dplyr::select("activitytime_radians")

#subset radian activity times by species (opossum)
dv_radians <- records[records$Species == "Dv", ]
dv_radians <- dv_radians %>% dplyr::select("activitytime_radians")

#Assign radian times to each species name
Pl <- pl_radians
Dv <- dv_radians

#Convert radian times to vector for each species
pl_rad_vect<- unlist(pl_radians, use.names=FALSE)
dv_rad_vect <- unlist(dv_radians, use.names=FALSE)

#Assign radian times to each species name
Pl <- pl_rad_vect
Dv <- dv_rad_vect

# Bootstrap analysis
bootstrap <- bootstrap(Pl, Dv, 1000, type="Dhat4")
mean(bootstrap)

bootstrap2 <- bootstrap(Pl, Dv, 10000, type="Dhat4")

bootCI(0.79, bootstrap, conf=0.95)

#Create png image of plot
png("/yourpath", width=3.25,
    height=3.25,
    units="in",
    res=1200,
    pointsize=6)         
overlapPlot(Pl, Dv, main="", xscale = 24, extend=NULL, linet=c(1,2), linecol=c(1,"blue"), xlab="Time of Day")
legend('topleft', c("Pl", "Dv"), lty=c(1,2), col=c(1,"blue"), bty='n')
dev.off()






