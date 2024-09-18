# Re-analysis of Delavaux et al. 2024
# Changes in the original Delavaux scripts include:
# - Myc_analyses.R: Save data that was used to produce the world map (Figure 1b in the original manuscript)
# - Nfix_Analyses.R, Myc_Analyses.R, and Poll_Analyses.R: Add Random Forest estimator for expected mutualism ratios
# - Nfix_Analyses.R, Myc_Analyses.R, and Poll_Analyses.R: Return observed mutualism ratio on islands
# - Joint_Analyses.R:  Calculate observed mutualism ratio on islands (new variable called biotic.ml_obs)
# - Joint_Analyses.R:  Calculate predicted/corrected mutualism ratio on islands (new variable called biotic.ml_rf)
# - Joint_Analyses.R:  Save data that was used to calculate the effects (e.g. Figure 2c in the original manuscript)


library(mgcv); library(gridExtra); library(betareg); library(MASS); library(lme4); library(lmerTest); library(lsmeans); library(ggeffects); library(spdep); library(ggplot2); library(ncf); library(ape); library(sjPlot); library(gridExtra); library(MuMIn); library(tidyverse); library(maps); library(sf); library(tidyverse); library(relaimpo);library(spdep);library(randomForest);library(ggtext)
options(na.action = "na.fail")

# generated in the "Joint_Analyses.R" script
dat = readRDS("data/reanalysis_data.RDS")

# remove "scale" attributes 
dat$biotic.ml = as.numeric(dat$biotic.ml) 
dat$biotic.ml_rf = as.numeric(dat$biotic.ml_rf) 

# Rsquareds for RF and the original GAMs 
r2_RF_Myc = readRDS("data/r2_RF_Myc.RDS")
r2_RF_Nfix = readRDS("data/r2_RF_Nfix.RDS")
r2_RF_Poll = readRDS("data/r2_RF_Poll.RDS")
results_random_forest = 
  rbind(r2_RF_Myc, r2_RF_Nfix, r2_RF_Poll)

# Average R2 for RF and GAM for the mutualism predictions, FULL == Random Forest
colMeans(results_random_forest[c(1, 5, 7),c(4, 5)])


# RAC function
Spat.cor <- function(mod,dat, dist) {
  coords <- cbind(dat$longitude, dat$latitude)
  matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat = T)
  return(rac)
}

# RAC function when locations repeat (shift latlon)
Spat.cor.rep <- function(mod,dat, dist) {
  coords <- cbind(dat$longitude, dat$latitude) + matrix(runif(2*nrow(dat), 0, 0.00001), nrow = nrow(dat), ncol = 2)
  matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat = T)
  return(rac)
}



#### Plots

### Figure 1

# Models
## Original model 
sprichdiff.mod <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + biotic.ml, data = dat) 
rac <- Spat.cor.rep(sprichdiff.mod, dat, 2000)
sprichdiff.mod.rac <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + rac + biotic.ml, data = dat) 
summary(sprichdiff.mod.rac)

## Our model with the new biotic.ml_rf variable
sprichdiff.mod_rf <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + biotic.ml_rf, data = dat) 
rac <- Spat.cor.rep(sprichdiff.mod_rf, dat, 2000)
sprichdiff.mod_rf_rac <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + rac + biotic.ml_rf, data = dat) 
summary(sprichdiff.mod_rf_rac)

mod.dat <- summary(sprichdiff.mod.rac)$coefficients %>% as.data.frame() %>%
  mutate(variable = rownames(.)) %>%
  filter(! variable=="(Intercept)") %>%
  filter(! variable=="rac") %>%
  mutate(variable = case_when(variable=="abslatitude" ~ "Absolute latitude",
                              variable=="area" ~ "Area", 
                              variable=="biotic.ml" ~ "Mutualism filter strength",
                              variable=="dist" ~ "Distance",
                              variable=="prec" ~ "Precipitation",
                              variable=="elev_range" ~ "Elevation range")) %>%
  filter(!is.na(variable)) 

colnames(mod.dat)<-c("est", "std.err","tval","pval", "variable")

mod.dat.ordered <- mod.dat
mod.dat.ordered$variable = fct_reorder(mod.dat$variable, mod.dat$est)
order_effects = fct_reorder(mod.dat$variable, mod.dat$est)

forest.plot_original <-
  ggplot(data=mod.dat.ordered, aes(x=variable, y=est, ymin=est-std.err, ymax=est+std.err), fill = "coral3") +
  geom_pointrange(alpha = 0.8, size = 2, color="coral3") + 
  geom_hline(yintercept=0, lty=2, color='darkgrey') +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab(" ") + ylab("Effect of 1 SD of predictor \non species deficit") +
  theme_classic(base_size = 40) +
  theme(legend.position = 'none') +
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(size = 30), axis.text.y = element_text(angle = 0, size = 30)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))  + labs(tag = "B")


mod.dat <- summary(sprichdiff.mod_rf_rac)$coefficients %>% as.data.frame() %>%
  mutate(variable = rownames(.)) %>%
  filter(! variable=="(Intercept)") %>%
  filter(! variable=="rac") %>%
  mutate(variable = case_when(variable=="abslatitude" ~ "Absolute latitude",
                              variable=="area" ~ "Area", 
                              variable=="biotic.ml_rf" ~ "Mutualism filter strength",
                              variable=="dist" ~ "Distance",
                              variable=="prec" ~ "Precipitation",
                              variable=="elev_range" ~ "Elevation range")) %>%
  filter(!is.na(variable)) 

colnames(mod.dat)<-c("est", "std.err","tval","pval", "variable")
mod.dat.ordered <- mod.dat
mod.dat.ordered$variable = order_effects

forest.plot_original_rf <-
  ggplot(data=mod.dat.ordered, aes(x=variable, y=est, ymin=est-std.err, ymax=est+std.err), fill = "coral3") +
  geom_pointrange(alpha = 0.8, size = 2, color="coral3") + 
  geom_hline(yintercept=0, lty=2, color='darkgrey') +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab(" ") + ylab("Effect of 1 SD of predictor \non species deficit") +
  theme_classic(base_size = 40) +
  theme(legend.position = 'none') +
  #theme(legend.justification=c(1,1), legend.position=c(1,1))+
  guides(fill = FALSE) +
  theme(axis.text.x = element_text(size = 30), axis.text.y = element_text(angle = 0, size = 30)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 15))   + labs(tag = "D")


## Maps
df_map = readRDS("data/data_for_map.RDS")

rf = randomForest(AM~latitude+longitude, data = df_map$dat.ml, mtry = 2L, ntree = 1000L)
gam_lat = gam(AM~s(latitude), data = df_map$dat.ml)
dat.pred_rf = df_map$dat.ml
dat.pred_rf$AM = predict(rf, newdata =df_map$dat.ml)
dat.pred_gam = df_map$dat.ml
dat.pred_gam$AM = predict(gam_lat, newdata = df_map$dat.ml)

map_gam = 
  ggplot()+
  geom_polygon(data = df_map$world, aes(x = long, y = lat, group = group), fill = "gray10", alpha = 0.5) +
  geom_point(data = dat.pred_gam, aes(x = longitude, y = latitude, color = AM, size = AM),  
             pch = 19, alpha = 0.8, stroke = 0) +
  scale_color_viridis(option = "D", begin = 0.15, end = 1, alpha = 0.5) +
  scale_size_continuous(range = c(2, 7)) +
  geom_text(aes(x = -168, y = 5, label = paste0("R^2 == ", round(r2_RF_Myc[1,]$GAM, 3))), parse = TRUE, fontface = "bold", size = 10) +
  xlab("") + ylab ("") +
  theme_minimal() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), 
        legend.position = "bottom", legend.key.width = unit(1,"cm"),
        axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(color = "Predicted AM richness") +
  guides(size = "none")  +
  labs(tag = "A")+
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(angle = 45, size = 20), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        text = element_text(size=30)) 


map_rf = 
  ggplot()+
  geom_polygon(data = df_map$world, aes(x = long, y = lat, group = group), fill = "gray10", alpha = 0.5) +
  geom_point(data = dat.pred_rf, aes(x = longitude, y = latitude, color = AM, size = AM),  
             pch = 19, alpha = 0.8, stroke = 0) +
  scale_color_viridis(option = "D", begin = 0.15, end = 1, alpha = 0.5) +
  scale_size_continuous(range = c(2, 7)) +
  geom_text(aes(x = -168, y = 5, label = paste0("R^2 == ", round(r2_RF_Myc[1,]$Full, 3))), parse = TRUE, fontface = "bold", size = 10) +
  xlab("") + ylab ("") +
  theme_minimal() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), 
        legend.position = "bottom", legend.key.width = unit(1,"cm"),
        axis.text.x = element_blank(), axis.title.x = element_blank()) +
  labs(color = "Predicted AM richness") +
  guides(size = "none")  + labs(tag = "C") + 
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(angle = 45, size = 20), 
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        text = element_text(size=30)) 
grob <- arrangeGrob(map_gam, forest.plot_original, map_rf, forest.plot_original_rf, ncol = 2)
ggsave("figures/Comparison_models.png",grob, width = 26, height = 16)



### Figure 2
# Set mean to zero (because of the effect plots), no impact on results/effects
dat = 
  dat |> mutate(abslatitude = abslatitude - mean(abslatitude),
                area = area - mean(area),
                dist = dist - mean(dist),
                elev_range = elev_range - mean(elev_range),
                prec = prec - mean(prec),
                biotic.ml = biotic.ml - mean(biotic.ml),
                biotic.ml_rf = biotic.ml_rf - mean(biotic.ml_rf))

# Model from Delavaux et al
## Linear effect on biotic.ml
sprichdiff.mod <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + biotic.ml, data = dat) 
rac <- Spat.cor.rep(sprichdiff.mod, dat, 2000)
rac = rac 
sprichdiff.mod.rac <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + rac + biotic.ml, data = dat) 
summary(sprichdiff.mod.rac)

## Splines on all terms
sprichdiff.mod_gam <- gam(sprichdiff ~ s(abslatitude) + s(area) + s(dist) + s(elev_range) +  s(prec) + s(biotic.ml), data = dat) 
rac_gam <- Spat.cor.rep(sprichdiff.mod_gam, dat, 2000)
rac_gam = rac_gam 
sprichdiff.mod.rac_gam <- gam(sprichdiff ~ s(abslatitude) + s(area) + s(dist) + s(elev_range) +  s(prec) + s(rac_gam) + s(biotic.ml), data = dat) 
summary(sprichdiff.mod.rac_gam)

## Spline on all terms except for biotic.ml
sprichdiff.mod_gam_linear <- gam(sprichdiff ~ s(abslatitude) + s(area) + s(dist) + s(elev_range) +  s(prec) + biotic.ml, data = dat) 
rac_gam_linear <- Spat.cor.rep(sprichdiff.mod_gam_linear, dat, 2000)
sprichdiff.mod.rac_gam_linear <- gam(sprichdiff ~ s(abslatitude) + s(area) + s(dist) + s(elev_range) +  s(prec) + s(rac_gam_linear) + biotic.ml, data = dat) 
summary(sprichdiff.mod.rac_gam_linear)


# New models / new biotic.ml variable (based on random forest)
## Linear effect on biotic.ml_rf
sprichdiff.mod_rf <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + biotic.ml_rf, data = dat) 
rac_rf <- Spat.cor.rep(sprichdiff.mod_rf, dat, 2000)
sprichdiff.mod.rac_rf <- glm(sprichdiff ~ abslatitude + area + dist + elev_range +  prec + rac_rf + biotic.ml_rf, data = dat) 
summary(sprichdiff.mod.rac_rf)

## Splines on all terms
sprichdiff.mod_rf_gam <- gam(sprichdiff ~ s(abslatitude) + s(area) + s(dist) + s(elev_range) +  s(prec) + s(biotic.ml_rf), data = dat) 
rac_rf_gam <- Spat.cor.rep(sprichdiff.mod_rf_gam, dat, 2000)
sprichdiff.mod.rac_rf_gam <- gam(sprichdiff ~ s(abslatitude) + s(area) + s(dist) + s(elev_range) +  s(prec) + s(rac_rf_gam) + s(biotic.ml_rf), data = dat) 

## Spline on all terms except for biotic.ml_rf
sprichdiff.mod_rf_gam_linear <- gam(sprichdiff ~ s(abslatitude) + s(area) + s(dist) + s(elev_range) +  s(prec) + biotic.ml_rf, data = dat) 
rac_rf_gam_linear <- Spat.cor.rep(sprichdiff.mod_rf_gam_linear, dat, 2000)
sprichdiff.mod.rac_rf_gam_linear <- gam(sprichdiff ~ s(abslatitude) + s(area) + s(dist) + s(elev_range) +  s(prec) + s(rac_rf_gam_linear) + biotic.ml_rf, data = dat) 
summary(sprichdiff.mod.rac_rf_gam_linear)


### Figure 2 
pdf("figures/figure_2.pdf", width = 14, height = 5)
par(mfrow = c(1, 3))
lwd = 1.2
plot(dat$abslatitude, dat$biotic.ml, pch = 16, col = "#00000055", xlab = "Mutualism Filter Strength", ylab = "Absolute Latitude", las = 1, ylim = c(-3.3, 1.8))
points(dat$abslatitude, dat$biotic.ml_rf, pch = 16, col = "#FF000075")
legend("bottomright", col = c("black", "red"), pch = 16,, lty = 1, bty = "n",legend = c("Delavaux et. al Mutualism filter", "Our Mutualism filter") )
text(x = -2, y = 2.1, pos = 3, labels = "A", font = 2, xpd = NA, cex = 1.5)


plot(dat$biotic.ml, dat$sprichdiff, pch = 16, col = "#00000055", xlab = "Mutualism Filter Strength", ylab = "Species Difference", las = 1, xlim = c(-5.2, 2))
text(x = -5.9, y = 3300, pos = 3, labels = "B", font = 2, xpd = NA, cex = 1.5)

points(dat$biotic.ml_rf, dat$sprichdiff, pch = 16, col = "#FF000075")

pred = predict(sprichdiff.mod, newdata = dat, type = "terms", se.fit = TRUE)

polygon(c(dat$biotic.ml[order(dat$biotic.ml)], 
          dat$biotic.ml[order(dat$biotic.ml, decreasing = TRUE)] ), 
        c((pred$fit[,6]-pred$se.fit[,6]+coef(sprichdiff.mod)[1])[order(dat$biotic.ml)],
          (pred$fit[,6]+pred$se.fit[,6]+coef(sprichdiff.mod)[1])[order(dat$biotic.ml, decreasing = TRUE)]  ), col = "#22222222", border = NA )
points(dat$biotic.ml[order(dat$biotic.ml)], (pred$fit[,6]+coef(sprichdiff.mod)[1])[order(dat$biotic.ml)], type = "l", col = "black", lwd = lwd, , lty = 1)

pred = predict(sprichdiff.mod_gam_linear, newdata = dat, type = "terms", se.fit = TRUE)
polygon(c(dat$biotic.ml[order(dat$biotic.ml)], 
          dat$biotic.ml[order(dat$biotic.ml, decreasing = TRUE)] ), 
        c((pred$fit[,1]-pred$se.fit[,1]+coef(sprichdiff.mod_gam_linear)[1])[order(dat$biotic.ml)],
          (pred$fit[,1]+pred$se.fit[,1]+coef(sprichdiff.mod_gam_linear)[1])[order(dat$biotic.ml, decreasing = TRUE)]  ), col = "#22222222", border = NA )
points(dat$biotic.ml[order(dat$biotic.ml)], (pred$fit[,1]+coef(sprichdiff.mod_gam_linear)[1])[order(dat$biotic.ml)], type = "l", col = "black", lwd = lwd, , lty = 2)


pred = predict(sprichdiff.mod_rf, newdata = dat, type = "terms", se.fit = TRUE)

polygon(c(dat$biotic.ml_rf[order(dat$biotic.ml_rf)], 
          dat$biotic.ml_rf[order(dat$biotic.ml_rf, decreasing = TRUE)] ), 
        c((pred$fit[,6]-pred$se.fit[,6]+coef(sprichdiff.mod_rf)[1])[order(dat$biotic.ml_rf)],
          (pred$fit[,6]+pred$se.fit[,6]+coef(sprichdiff.mod_rf)[1])[order(dat$biotic.ml_rf, decreasing = TRUE)]  ), col = "#FF000033", border = NA )
points(dat$biotic.ml_rf[order(dat$biotic.ml_rf)], (pred$fit[,6]+coef(sprichdiff.mod_rf)[1])[order(dat$biotic.ml_rf)], type = "l", col = "red", lwd = lwd, , lty = 1)

pred = predict(sprichdiff.mod_rf_gam_linear, newdata = dat, type = "terms", se.fit = TRUE)
polygon(c(dat$biotic.ml_rf[order(dat$biotic.ml_rf)], 
          dat$biotic.ml_rf[order(dat$biotic.ml_rf, decreasing = TRUE)] ), 
        c((pred$fit[,1]-pred$se.fit[,1]+coef(sprichdiff.mod_rf_gam_linear)[1])[order(dat$biotic.ml_rf)],
          (pred$fit[,1]+pred$se.fit[,1]+coef(sprichdiff.mod_rf_gam_linear)[1])[order(dat$biotic.ml_rf, decreasing = TRUE)]  ), col = "#FF000033", border = NA )
points(dat$biotic.ml_rf[order(dat$biotic.ml_rf)], (pred$fit[,1]+coef(sprichdiff.mod_rf_gam_linear)[1])[order(dat$biotic.ml_rf)], type = "l", col = "red", lwd = lwd, , lty = 2)

legend("bottomright", col = c("black", "red"), pch = 16,, lty = 1, bty = "n",legend = c("Delavaux et. al Mutualism filter", "Our Mutualism filter") )
legend("topleft", bty = "n", legend = c("LM", "GAM"), lty = c(1, 2), col = c("black", "black"), lwd = lwd)



gam_predicted = gam(biotic.ml_rf~s(abslatitude, k = 5), data = dat)
gam_obs = gam(biotic.ml_obs~s(abslatitude, k = 5), data = dat)
plot(dat$abslatitude, dat$biotic.ml_obs, pch = 16, col = "#00000055", ylab = "Mutualism Filter Strength", xlab = "Absolute Latitude", las = 1, ylim = c(-3.6, 1.8))
points(dat$abslatitude, dat$biotic.ml_rf, pch = 16, col = "#FF000075")

text(x = -2, y = 2.1, pos = 3, labels = "C", font = 2, xpd = NA, cex = 1.5)


pred = predict(gam_predicted, newdata = dat, type = "terms", se.fit = TRUE)

polygon(c(dat$abslatitude[order(dat$abslatitude)], 
          dat$abslatitude[order(dat$abslatitude, decreasing = TRUE)] ), 
        c((pred$fit[,1]-pred$se.fit[,1]+coef(gam_predicted)[1])[order(dat$abslatitude)],
          (pred$fit[,1]+pred$se.fit[,1]+coef(gam_predicted)[1])[order(dat$abslatitude, decreasing = TRUE)]  ), col = "#FF000033", border = NA )
points(dat$abslatitude[order(dat$abslatitude)], (pred$fit[,1]+coef(gam_predicted)[1])[order(dat$abslatitude)], type = "l", col = "red", lwd = lwd, , lty = 1)

pred = predict(gam_obs, newdata = dat, type = "terms", se.fit = TRUE) # [,6]+coef(sprichdiff.mod)[1])[order(dat$biotic.ml_rf)]
polygon(c(dat$abslatitude[order(dat$abslatitude)], 
          dat$abslatitude[order(dat$abslatitude, decreasing = TRUE)] ), 
        c((pred$fit[,1]-pred$se.fit[,1]+coef(gam_obs)[1])[order(dat$abslatitude)],
          (pred$fit[,1]+pred$se.fit[,1]+coef(gam_obs)[1])[order(dat$abslatitude, decreasing = TRUE)]  ), col = "#22222222", border = NA )
points(dat$abslatitude[order(dat$abslatitude)], (pred$fit[,1]+coef(gam_obs)[1])[order(dat$abslatitude)], type = "l", col = "black", lwd = lwd, , lty = 1)
legend("bottomleft", col = c("black", "red"), pch = 16,, lty = 1, bty = "n",legend = c("Observed Mutualism filter", "Corrected Mutualism filter") )
dev.off()




# Figure Appendix
pdf(file = "figures/S1.pdf", width = 12, height = 3) 
par(mfrow = c(1,3),mar = c(4, 4, 2, 1))
lwd = 2.0
# Area
plot(dat$area, dat$sprichdiff, pch = 16, col = "grey", xlab = "Area", ylab = "Species Richness", las = 1)
points(dat$area[order(dat$area)], (predict(sprichdiff.mod_rf, newdata = dat, type = "terms")[,2]+coef(sprichdiff.mod_rf)[1])[order(dat$area)], type = "l", col = "black", lwd = lwd, lty = 1)
points(dat$area[order(dat$area)], (predict(sprichdiff.mod_rf_gam, newdata = dat, type = "terms")[,2]+coef(sprichdiff.mod_rf_gam)[1])[order(dat$area)], type = "l", col = "black", lwd = lwd, lty = 2)
legend("bottomleft", bty = "n", legend = c("LM", "GAM"), lty = c(1, 2), col = c("black", "black"), lwd = lwd)
text(-1, y = 3400, xpd = NA, pos = 2, labels = "A", fon = 2, cex = 1.5)

# Distance
plot(dat$dist, dat$sprichdiff, pch = 16, col = "grey", xlab = "Distance", ylab = "", las = 1)
points(dat$dist[order(dat$dist)], (predict(sprichdiff.mod_rf, newdata = dat, type = "terms")[,3]+coef(sprichdiff.mod_rf)[1])[order(dat$dist)], type = "l", col = "black", lwd = lwd, lty = 1)
points(dat$dist[order(dat$dist)], (predict(sprichdiff.mod_rf_gam, newdata = dat, type = "terms")[,3]+coef(sprichdiff.mod_rf_gam)[1])[order(dat$dist)], type = "l", col = "black", lwd = lwd, lty = 2)
legend("bottomleft", bty = "n", legend = c("LM", "GAM"), lty = c(1, 2), col = c("black", "black"), lwd = lwd)
text(-1.2, y = 3400, xpd = NA, pos = 2, labels = "B", fon = 2, cex = 1.5)

# Abslatitude
plot(dat$abslatitude, dat$sprichdiff, pch = 16, col = "grey", xlab = "Abs Latitute", ylab = "", las = 1)
points(dat$abslatitude[order(dat$abslatitude)], (predict(sprichdiff.mod_rf, newdata = dat, type = "terms")[,1]+coef(sprichdiff.mod_rf)[1])[order(dat$abslatitude)], type = "l", col = "black", lwd = lwd, lty = 1)
points(dat$abslatitude[order(dat$abslatitude)], (predict(sprichdiff.mod_rf_gam, newdata = dat, type = "terms")[,1]+coef(sprichdiff.mod_rf_gam)[1])[order(dat$abslatitude)], type = "l", col = "black", lwd = lwd, lty = 2)
legend("bottomleft", bty = "n", legend = c("LM", "GAM"), lty = c(1, 2), col = c("black", "black"), lwd = lwd)
text(-1.7, y = 3400, xpd = NA, pos = 2, labels = "C", fon = 2, cex = 1.5)
dev.off()

