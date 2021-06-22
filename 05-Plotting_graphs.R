# NRI plots for paper
library(ggplot2)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(broom)
library(ggtext)


# Arrange data ------------------------------------------------------------

# To subset to years with more than 10 data
subsample_data <- function(p, d, smin){
  nd <- subset(d, year==p)
  if(nrow(nd) > smin){return(nd)}
}

data.cels.subset <- pcs.data

data.cels.subset <- data.cels.subset[, c("nri", "nri_p", "nri_significant",  "nti", "nti_p", "nti_significant", "year", "Tmin", "Tmax", "Pmin", "Pmax", "AET", "ETR", "WDI", "Deglac.", "x", "y")]
data.cels.subset <- lapply(paste(PERIODS, "BP", sep=""), FUN=subsample_data, data.cels.subset, 10)
data.cels.subset <- do.call(rbind, data.cels.subset)



# Scatter - regression plots -----------------------------------------------


plot_scatter_regression <- function(data, model, var, pcs_metric, years=NULL){
  # data <- data.cels.subset
  # model <- pcs.ols
  # var <- "Tmin"
  # pcs_metric <- "nri"

  if(var == "Tmin") x.label <- "Minimum temperature (ºC)"  
  if(var == "Tmax") x.label <- "Maximum temperature (ºC)"  
  if(var == "Pmin") x.label <- "Minimum precipitation (mm)"  
  if(var == "Pmax") x.label <- "Maximum precipitation (mm)"  
  if(var == "AET") x.label <- "Actual Evapotranspiration (mm)"  
  if(var == "ETR") x.label <- "Evapotranspiration Ratio"  
  if(var == "WDI") x.label <- "Water Deficit Index"  
  if(var == "Deglac."){
    x.label <- "Time since deglaciation (ka)"  
    data$Deglac. <- data$Deglac / 1000
  } 
  
  if(pcs_metric == "nri") y.label <- "NRI"
  if(pcs_metric == "nti") y.label <- "NTI"
  
  m1.r2 <- signif(broom::glance(model[[pcs_metric]][[var]][["m1"]])$r.squared, 2)
  m2.r2 <- signif(broom::glance(model[[pcs_metric]][[var]][["m2"]])$r.squared, 2)
  m3.r2 <- signif(broom::glance(model[[pcs_metric]][[var]][["m3"]])$r.squared, 2)
  
  m1.p <- signif(broom::glance(model[[pcs_metric]][[var]][["m1"]])$p.value[[1]], 3)
  m2.p <- signif(broom::glance(model[[pcs_metric]][[var]][["m2"]])$p.value[[1]], 3)
  m3.p <- signif(broom::glance(model[[pcs_metric]][[var]][["m3"]])$p.value[[1]], 3)
  
  m1.interc <- coef(model[[pcs_metric]][[var]][["m1"]])[1]
  m1.slope <- coef(pcs.ols[[pcs_metric]][[var]][["m1"]])[2]
  m2.coef <- data.frame(year=unique(data$year))
  m2.coef$estimate <- summary(model[[pcs_metric]][[var]][["m2"]])$coefficients[2:20, 1]
  m2.coef$intercept <- summary(model[[pcs_metric]][[var]][["m2"]])$coefficients[1, 1]
  m2.coef$var.intercepts <- m2.coef$intercept + m2.coef$estimate
  m2.coef$var.intercepts[1] <- summary(model[[pcs_metric]][[var]][["m2"]])$coefficients[1, 1]
  m2.coef$slope <- summary(model[[pcs_metric]][[var]][["m2"]])$coefficients[2, 1]

  fig.a <- ggplot(data, aes_string(x=var, y=pcs_metric)) +
    stat_binhex(bins=30) + 
    geom_smooth(method="lm", color="#E69F00", fill="#E69F00", size=1, fullrange=TRUE) +
    theme_bw() +
    labs(x=x.label, y=y.label) +
    ggtitle(paste(expression("Adj R² = "), m1.r2, ", p =", m1.p)) +
    scale_fill_gradient(low="grey70", high="black", name="# of points", breaks=c(1, 50, 100, 150), limits=c(1,150)) +
    coord_cartesian(ylim=c(-3.5, 3.5)) +
    theme(panel.grid.minor=element_blank(), plot.title=element_text(color="#E69F00", size=9))
  
  if(!is.null(years)){
    data <- subset(data, year %in% years)
    m2.coef <- subset(m2.coef, year %in% years)
  }
  
  fig.b <- ggplot(data, aes_string(x=var, y=pcs_metric)) +
    stat_binhex(bins=30) +
    geom_abline(intercept=m1.interc, slope=m1.slope, colour="#E69F00", size=1) +
    geom_abline(data=m2.coef, aes(intercept=intercept, slope=slope), color = "#009E73", size=1) +
    geom_smooth(method="lm", color="#56B4E9", fill="royalblue4", size=1, fullrange=TRUE) +
    facet_wrap(~year) +
    theme_bw() +
    labs(x=x.label,
         y=y.label,
         title=paste(
           expression("<span style = 'color: #009E73;'>Adj R² = "),
           m2.r2,
           ", p =",
           m2.p,
           "</span>- <span style = 'color: #56B4E9;'>",
           expression("Adj R² = "),
           m3.r2,
           ", p =",
           m3.p,
           "</span>")) +
    scale_fill_gradient(low="grey70", high="black", name="# of points", breaks=c(1, 50, 100, 150), limits=c(1,150)) +
    coord_cartesian(ylim=c(-3.5, 3.5)) +
    theme(panel.grid.minor=element_blank(), plot.title = element_markdown(size=9))
  
  fig <- list(all=fig.a, panels=fig.b)
  return(fig)
}

arrange_scatterplots <- function(figs, arrange = "vertical"){
  if(arrange == "vertical"){
    ncol <- 1
    nrow <- 2
    legend <- "bottom"
    widths <- c(1, 1)
  }else{
    ncol <- 2
    nrow <- 1
    legend <- "bottom"
    widths <- c(1, 1.5)
  } 
  fig <- ggarrange(figs$all, figs$panels, ncol=ncol, nrow=nrow, labels=c("a)", "b)"), common.legend=TRUE, legend=legend, widths=widths)
  return(fig)
} 


fig.nri <- lapply(VARS_NEW_NAMES, FUN=plot_scatter_regression, data=data.cels.subset, model=pcs.ols, pcs_metric="nri")

fig.nri <- lapply(fig.nri, FUN=arrange_scatterplots)

mapply(FUN=ggexport, fig.nri, filename = paste0("Results/Figures/Fig-scatter_regression/Fig-NRI_", VARS_NEW_NAMES, ".pdf"), MoreArgs = list(width=6, height=8))

fig.nti <- lapply(VARS_NEW_NAMES, FUN=plot_scatter_regression, data=data.cels.subset, model=pcs.ols, pcs_metric="nti")

fig.nti <- lapply(fig.nti, FUN=arrange_scatterplots)

mapply(FUN=ggexport, fig.nti, filename = paste0("Results/Figures/Fig-scatter_regression/Fig-NTI_", VARS_NEW_NAMES, ".pdf"), MoreArgs = list(width=6, height=8))

nri.figs.list <- paste0("Results/Figures/Fig-scatter_regression/Fig-NRI_", VARS_NEW_NAMES, ".pdf")
system(paste0("qpdf --empty --pages ", paste(nri.figs.list, collapse=" "), " -- Results/Manuscript_figures_and_tables/Supplement_A3.pdf"))

nti.figs.list <- paste0("Results/Figures/Fig-scatter_regression/Fig-NTI_", VARS_NEW_NAMES, ".pdf")
system(paste0("qpdf --empty --pages ", paste(nti.figs.list, collapse=" "), " -- Results/Manuscript_figures_and_tables/Supplement_A4.pdf"))



years <- c("0BP", "5000BP", "12000BP", "16000BP", "17000BP", "18000BP")

fig_2 <- lapply(c("Tmin", "Pmin", "ETR"), FUN=plot_scatter_regression, data=data.cels.subset, model=pcs.ols, years=years, pcs_metric="nri")

fig_2 <- ggarrange(fig_2[[1]]$all, fig_2[[1]]$panels, fig_2[[2]]$all, fig_2[[2]]$panels, fig_2[[3]]$all, fig_2[[3]]$panels, ncol=2, nrow=3, labels=c("a)", "b)", "c)", "d)", "e)", "f)"), common.legend=TRUE, legend="bottom", widths=c(1, 1.5))

fig_2

ggexport(fig_2, filename = "Results/Manuscript_figures_and_tables/Figure_2.pdf", width=7, height=8)


# Deglaciation ------------------------------------------------------------

# REGRESSION PLOTS

# NRI
ggplot(data.cels.subset, aes(x=as.factor(Deglac./1000), y=nri)) +
  geom_point(alpha=0.25, colour="red") +
  geom_smooth(method=lm, color="black", size=1) +
  facet_wrap(~year) +
  labs(x="Time since deglaciation (ka)", y="NRI") +
  theme_bw() +
  theme(panel.grid.minor=element_blank())


# NTI 
ggplot(data.cels.subset, aes(x=as.factor(Deglac./1000), y=nti)) +
  geom_point(alpha=0.25, colour="red") +
  geom_smooth(method=lm, color="black", size=1) +
  facet_wrap(~year) +
  labs(x="Time since deglaciation (ka)", y="NTI") +
  theme_bw() +
  theme(panel.grid.minor=element_blank())


# HEATMAPS 

# NRI 
fig.a <- ggplot(data.cels.subset, aes(x=year, y=as.factor(Deglac./1000), z=nri)) +
  stat_summary_2d() +
  # scale_fill_gradient2(limits=c(-4, 4), low="red", mid="#ffffbf", high="blue", name="mean\nNRI") + 
  scale_fill_gradient2(low="red", mid="#ffffbf", high="blue", name="mean", limits=c(-1.7, 1.7), breaks=c(-1.5, 0, 1.5)) + 
  coord_equal() +
  labs(x="Time period", y="Time since deglaciation (ka)") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_rect(fill="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank())


# NTI 
fig.b <- ggplot(data.cels.subset, aes(x=year, y=as.factor(Deglac./1000), z=nti)) +
  stat_summary_2d() +
  # scale_fill_gradient2(limits=c(-4, 4), low="red", mid="#ffffbf", high="blue", name="mean\nNTI") + 
  scale_fill_gradient2(low="red", mid="#ffffbf", high="blue", name="mean", limits=c(-1.7, 1.7), breaks=c(-1.5, 0, 1.5)) + 
  coord_equal() +
  labs(x="Time period", y="Time since deglaciation (ka)") +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_rect(fill="black"),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())


# NRI ~ climate residuals 
nri.model.start <- lm(nri ~ 1, data=data.cels.subset)
nri.form.step <- formula(paste0("nri ~ ", paste(CLIMVARS_NEW_NAMES, collapse=" + ")))
nri.model.back <- step(nri.model.start, nri.form.step, direction="both")
nri.residuals <- residuals(nri.model.back)

fig.c <- ggplot(data.cels.subset, aes(x=year, y=as.factor(Deglac./1000), z=nri.residuals)) +
  stat_summary_2d() +
  # scale_fill_gradient2(limits=c(-1.8, 1.8), low="red", mid="#ffffbf", high="blue", name="mean\nNRI res.") + 
  scale_fill_gradient2(low="red", mid="#ffffbf", high="blue", name="mean", limits=c(-1.7, 1.7), breaks=c(-1.5, 0, 1.5)) + 
  coord_equal() +
  labs(x="Time period", y="Time since deglaciation (ka)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_rect(fill="black"))


# NTI ~ climate residuals 
nti.model.start <- lm(nti ~ 1, data=data.cels.subset)
nti.form.step <- formula(paste0("nti ~ ", paste(CLIMVARS_NEW_NAMES, collapse=" + ")))
nti.model.back <- step(nti.model.start, nti.form.step, direction="both")
nti.residuals <- residuals(nti.model.back)

fig.d <- ggplot(data.cels.subset, aes(x=year, y=as.factor(Deglac./1000), z=nti.residuals)) +
  stat_summary_2d() +
  # scale_fill_gradient2(limits=c(-1.8, 1.8), low="red", mid="#ffffbf", high="blue", name="mean\nNTI res.") + 
  scale_fill_gradient2(low="red", mid="#ffffbf", high="blue", name="mean", limits=c(-1.7, 1.7), breaks=c(-1.5, 0, 1.5)) + 
  coord_equal() +
  labs(x="Time period", y="Time since deglaciation (ka)") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_rect(fill="black"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank())

fig <- ggarrange(fig.a, fig.b, fig.c, fig.d, ncol=2, nrow=2, labels=c("a)", "b)", "c)", "d)"), common.legend=TRUE, legend="right", align="hv")
ggexport(fig, filename = "Results/Figures/Deglaciation/Heatmaps.pdf", width=7.5, height=7.5)
ggexport(fig, filename = "Results/Manuscript_figures_and_tables/Figure_5.pdf", width=7.5, height=7.5)


# Plot intercepts ---------------------------------------------------------


# Load data for global temperature from ice cores
global.temp <- read.csv("../Data/Published-Paleoclimate/gicc05-60ka-20yr-short.csv")
global.temp <- global.temp[which(global.temp$Age < 18000), c(1, 3)]
colnames(global.temp) <- c("Time","delta18O")
global.temp$Time <- global.temp$Time/1000
index <- .bincode(global.temp$Time, c(0, seq(0.5, 17.5, by=1), 18), FALSE)
global.temp.avg <- aggregate(global.temp$delta18O, list(index), FUN=mean)
global.temp.avg <- data.frame(Time=seq(0, 18, by=1), delta18O=global.temp.avg$x)
global.temp.avg$dif <- c(rev(diff(rev(global.temp.avg$delta18O))), 0)
global.temp.avg$delta18O_mod <- (global.temp.avg$delta18O + 45) / 2.5
global.temp.avg$dif_mod <- (global.temp.avg$dif / 2) - 40
global.temp.avg <- melt(global.temp.avg, id.vars="Time")
global.temp.avg$linetype <- as.factor(c(rep(1, 19), rep(2, 19), rep(1, 19), rep(2, 19)))

extract_intercepts_ols <- function(mods, periods){

  years <- paste0("year", periods, "BP")
  
  mod <- mods$m1
  modint <- coef(mod)[which(names(coef(mod)) == "(Intercept)")]
  result.m1 <- data.frame(intercept=modint, years=NA, model_type="ols")

  mod <- mods$m2
  int.0 <- coef(mod)[which(names(coef(mod)) == "(Intercept)")]
  int.part <- coef(mod)[which(names(coef(mod)) %in% years)]
  modint <- c(int.0, int.0 + int.part)
  result.m2 <- data.frame(intercept=modint, years=periods, model_type="ols")

  mod <- mods$m3
  int.0 <- coef(mod)[which(names(coef(mod)) == "(Intercept)")]
  int.part <- coef(mod)[which(names(coef(mod)) %in% years)]
  modint <- c(int.0, int.0 + int.part)
  result.m3 <- data.frame(intercept=modint, years=periods, model_type="ols")
  
  result <- list(result.m1, result.m2, result.m3)
  names(result) <- c("m1", "m2", "m3")
  return(result)
}

extract_intercepts_sar <- function(mods, periods){
  
  years <- paste0("year", periods, "BP")

  mod <- mods$m1
  modint <- coef(mod)[which(names(coef(mod)) == "(Intercept)")]
  result.m1 <- data.frame(intercept=modint, years=NA, model_type="sar")

  mod <- mods$m2
  int.0 <- coef(mod)[which(names(coef(mod)) == "(Intercept)")]
  int.part <- coef(mod)[which(names(coef(mod)) %in% years)]
  modint <- c(int.0, int.0 + int.part)
  result.m2 <- data.frame(intercept=modint, years=periods, model_type="sar")

  mod <- mods$m3
  int.0 <- coef(mod)[which(names(coef(mod)) == "(Intercept)")]
  int.part <- coef(mod)[which(names(coef(mod)) %in% years)]
  modint <- c(int.0, int.0 + int.part)
  result.m3 <- data.frame(intercept=modint, years=periods, model_type="sar")
  
  result <- list(result.m1, result.m2, result.m3)
  names(result) <- c("m1", "m2", "m3")
  return(result)
}

ols.interc <- list()
ols.interc$nri <- lapply(pcs.ols$nri, FUN=extract_intercepts_ols, PERIODS[1:19])
ols.interc$nti <- lapply(pcs.ols$nti, FUN=extract_intercepts_ols, PERIODS[1:19])

sar.lag.interc <- list()
sar.lag.interc$nri <- list()
sar.lag.interc$nri$"120" <- lapply(pcs.sar.lag$nri$"120", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.lag.interc$nri$"360" <- lapply(pcs.sar.lag$nri$"360", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.lag.interc$nri$"480" <- lapply(pcs.sar.lag$nri$"480", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.lag.interc$nti <- list()
sar.lag.interc$nti$"120" <- lapply(pcs.sar.lag$nti$"120", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.lag.interc$nti$"360" <- lapply(pcs.sar.lag$nti$"360", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.lag.interc$nti$"480" <- lapply(pcs.sar.lag$nti$"480", FUN=extract_intercepts_sar, PERIODS[1:19])

sar.err.interc <- list()
sar.err.interc$nri <- list()
sar.err.interc$nri$"120" <- lapply(pcs.sar.err$nri$"120", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.err.interc$nri$"360" <- lapply(pcs.sar.err$nri$"360", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.err.interc$nri$"480" <- lapply(pcs.sar.err$nri$"480", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.err.interc$nti <- list()
sar.err.interc$nti$"120" <- lapply(pcs.sar.err$nti$"120", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.err.interc$nti$"360" <- lapply(pcs.sar.err$nti$"360", FUN=extract_intercepts_sar, PERIODS[1:19])
sar.err.interc$nti$"480" <- lapply(pcs.sar.err$nti$"480", FUN=extract_intercepts_sar, PERIODS[1:19])

ols.interc <- melt(ols.interc, id.vars=c("intercept", "years", "model_type"))
sar.lag.interc <- melt(sar.lag.interc, id.vars=c("intercept", "years", "model_type"))
sar.err.interc <- melt(sar.err.interc, id.vars=c("intercept", "years", "model_type"))

sar.lag.interc["model_type"] <- "sar.lag"
sar.err.interc["model_type"] <- "sar.err"

names(ols.interc) <- c("intercept", "years", "model_type", "model", "variable", "pcs_metric")
names(sar.lag.interc) <- names(sar.err.interc) <- c("intercept", "years", "model_type", "model", "variable", "distance", "pcs_metric")

ols.interc$distance <- NA

intercepts <- rbind(ols.interc, sar.lag.interc, sar.err.interc)
intercepts$inter <- interaction(intercepts$model, intercepts$variable)
intercepts$linetype <- "solid"

epoch.scale <- data.frame(name=c("Holocene","Pleistocene"), start=c(11.7, 18), end=c(0,11.7))
# age.scale <- data.frame(name=c("YD","B-A"), start=c(12.7,14.7), end=c(11.5,12.8))
epoch.scale$center <- (epoch.scale$start + epoch.scale$end)/2
# age.scale$center <- (age.scale$start + age.scale$end)/2

intercepts$model_type <- factor(intercepts$model_type, levels=c("ols", "sar.lag", "sar.err"), labels=c("OLS", "SARlag", "SARerr"))
intercepts$pcs_metric <- factor(intercepts$pcs_metric, levels=c("nri", "nti"), labels=c("NRI", "NTI"))
intercepts$variable <- factor(intercepts$variable, levels=c("x", "y", "Tmin", "Tmax", "Pmin", "Pmax", "AET", "ETR", "WDI", "Deglac."))

a <- ggplot(subset(intercepts, model == "m2" &
                     model_type %in% c("OLS", "SARerr") &
                     variable %in% VARS_NEW_NAMES &
                     ((pcs_metric == "NRI" & distance %in% c(NA, 480))|(pcs_metric == "NTI" & distance %in% c(NA, 360))))) +
  
  # geom_rect(data=age.scale, aes(xmin=start, xmax=end, ymin=-1.5, ymax=1.5), colour="gray20", fill="black", alpha=0.05, linetype=2) +
  geom_rect(data=epoch.scale, aes(xmin=start, xmax=end, ymin=-1.5, ymax=-1.1), colour="black", fill="white") +
  geom_text(data=epoch.scale, aes(x=center, y=-1.3, label=name), colour="black", size=4) +
  
  geom_line(aes(x=years/1000, y=intercept, colour=variable, group=inter), size=0.8) +
  
  scale_x_reverse(expand=c(0,0), name="", labels=NULL) +
  scale_y_continuous(expand=c(0,0), name="Intercepts (model 2)") +
  scale_colour_discrete(name="") +
  facet_grid(pcs_metric ~ model_type) +
  coord_cartesian() +
  expand_limits(y=c(-1.5, 1.5)) +
  theme_bw() 

b <- ggplot(subset(intercepts, model == "m3" &
                     model_type %in% c("OLS", "SARerr") &
                     variable %in% VARS_NEW_NAMES &
                     ((pcs_metric == "NRI" & distance %in% c(NA, 480))|(pcs_metric == "NTI" & distance %in% c(NA, 360))))) +
  
  # geom_rect(data=age.scale, aes(xmin=start, xmax=end, ymin=-3.5, ymax=7.5), colour="gray20", fill="black", alpha=0.05, linetype=2) +
  geom_rect(data=epoch.scale, aes(xmin=start, xmax=end, ymin=-5, ymax=-3.3), colour="black", fill="white") +
  geom_text(data=epoch.scale, aes(x=center, y=-4.15, label=name), colour="black", size=4) +
  
  geom_line(aes(x=years/1000, y=intercept, colour=variable, group=inter), size=0.8) +
  
  scale_x_reverse(expand=c(0,0), name="Time (ka BP)") +
  scale_y_continuous(expand=c(0,0), name="Intercepts (model 3)") +
  scale_colour_discrete(name="") +
  facet_grid(pcs_metric ~ model_type) +
  expand_limits(y=c(-5, 8)) +
  theme_bw() 

fig.intercepts <- ggarrange(a, b, ncol=1, nrow=2, labels=c("a)", "b)"), align = "v", common.legend=TRUE, legend="top")
ggexport(fig.intercepts, filename = "Results/Figures/Regression_through_time/Intercepts.pdf", height=7.5, width=6.5)
ggexport(fig.intercepts, filename = "Results/Manuscript_figures_and_tables/Figure_3.pdf", height=7.5, width=6.5)



# Plot slopes ---------------------------------------------------------

extract_slopes_ols <- function(var, mods, periods){
  
  years <- paste0("year", periods, "BP")
  
  mod <- mods[[var]]$m1
  modslp <- coef(mod)[which(names(coef(mod)) == var)]
  result.m1 <- data.frame(slope=modslp, years=NA, model_type="ols")
  
  mod <- mods[[var]]$m2
  modslp <- coef(mod)[which(names(coef(mod)) == var)]
  result.m2 <- data.frame(slope=modslp, years=NA, model_type="ols")
  
  mod <- mods[[var]]$m3
  slp.0 <- coef(mod)[which(names(coef(mod)) == var)]
  slp.part <- coef(mod)[which(names(coef(mod)) %in% paste0(var, ":", years))]
  modslp <- c(slp.0, slp.0 + slp.part)
  result.m3 <- data.frame(slope=modslp, years=periods, model_type="ols")
  
  result <- list(result.m1, result.m2, result.m3)
  names(result) <- c("m1", "m2", "m3")
  return(result)
}

extract_slopes_sar <- function(var, mods, periods){
  
  years <- paste0("year", periods, "BP")
  
  mod <- mods[[var]]$m1
  modslp <- coef(mod)[which(names(coef(mod)) == var)]
  result.m1 <- data.frame(slope=modslp, years=NA, model_type="sar")

  mod <- mods[[var]]$m2
  modslp <- coef(mod)[which(names(coef(mod)) == var)]
  result.m2 <- data.frame(slope=modslp, years=NA, model_type="sar")
  
  mod <- mods[[var]]$m3
  slp.0 <- coef(mod)[which(names(coef(mod)) == var)]
  slp.part <- coef(mod)[which(names(coef(mod)) %in% paste0(var, ":", years))]
  modslp <- c(slp.0, slp.0 + slp.part)
  result.m3 <- data.frame(slope=modslp, years=periods, model_type="sar")

  result <- list(result.m1, result.m2, result.m3)
  names(result) <- c("m1", "m2", "m3")
  return(result)
}


ols.slope <- list()
ols.slope$nri <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_ols, mods=pcs.ols$nri, periods=PERIODS[1:19])
ols.slope$nti <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_ols, mods=pcs.ols$nti, periods=PERIODS[1:19])
names(ols.slope$nri) <- names(ols.slope$nti) <- VARS_NEW_NAMES

sar.lag.slope <- list()
sar.lag.slope$nri <- list()
sar.lag.slope$nri$"120" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.lag$nri$"120", periods=PERIODS[1:19])
sar.lag.slope$nri$"360" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.lag$nri$"360", periods=PERIODS[1:19])
sar.lag.slope$nri$"480" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.lag$nri$"480", periods=PERIODS[1:19])
names(sar.lag.slope$nri$"120") <- names(sar.lag.slope$nri$"360") <- names(sar.lag.slope$nri$"480") <- VARS_NEW_NAMES
sar.lag.slope$nti <- list()
sar.lag.slope$nti$"120" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.lag$nti$"120", periods=PERIODS[1:19])
sar.lag.slope$nti$"360" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.lag$nti$"360", periods=PERIODS[1:19])
sar.lag.slope$nti$"480" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.lag$nti$"480", periods=PERIODS[1:19])
names(sar.lag.slope$nti$"120") <- names(sar.lag.slope$nti$"360") <- names(sar.lag.slope$nti$"480") <- VARS_NEW_NAMES

sar.err.slope <- list()
sar.err.slope$nri <- list()
sar.err.slope$nri$"120" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.err$nri$"120", periods=PERIODS[1:19])
sar.err.slope$nri$"360" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.err$nri$"360", periods=PERIODS[1:19])
sar.err.slope$nri$"480" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.err$nri$"480", periods=PERIODS[1:19])
names(sar.err.slope$nri$"120") <- names(sar.err.slope$nri$"360") <- names(sar.err.slope$nri$"480") <- VARS_NEW_NAMES
sar.err.slope$nti <- list()
sar.err.slope$nti$"120" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.err$nti$"120", periods=PERIODS[1:19])
sar.err.slope$nti$"360" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.err$nti$"360", periods=PERIODS[1:19])
sar.err.slope$nti$"480" <- lapply(VARS_NEW_NAMES, FUN=extract_slopes_sar, mods=pcs.sar.err$nti$"480", periods=PERIODS[1:19])
names(sar.err.slope$nti$"120") <- names(sar.err.slope$nti$"360") <- names(sar.err.slope$nti$"480") <- VARS_NEW_NAMES


ols.slope <- melt(ols.slope, id.vars=c("slope", "years", "model_type"))
sar.lag.slope <- melt(sar.lag.slope, id.vars=c("slope", "years", "model_type"))
sar.err.slope <- melt(sar.err.slope, id.vars=c("slope", "years", "model_type"))

sar.lag.slope["model_type"] <- "sar.lag"
sar.err.slope["model_type"] <- "sar.err"

names(ols.slope) <-c("slope", "years", "model_type", "model", "variable", "pcs_metric")
names(sar.lag.slope) <- names(sar.err.slope) <- c("slope", "years", "model_type", "model", "variable", "distance", "pcs_metric")

ols.slope$distance <- NA

slopes <- rbind(ols.slope, sar.lag.slope, sar.err.slope)
slopes$inter <- interaction(slopes$model, slopes$variable)
slopes$linetype <- "solid"
slopes$model_type <- factor(slopes$model_type, levels=c("ols", "sar.lag", "sar.err"), labels=c("OLS", "SARlag", "SARerr"))
slopes$pcs_metric <- factor(slopes$pcs_metric, levels=c("nri", "nti"), labels=c("NRI", "NTI"))
slopes$variable <- factor(slopes$variable, levels=VARS_NEW_NAMES)

c <- ggplot(subset(slopes, model == "m3" &
                     model_type %in% c("OLS", "SARerr") &
                     variable %in% VARS_NEW_NAMES &
                     ((pcs_metric == "NRI" & distance %in% c(NA, 480))|(pcs_metric == "NTI" & distance %in% c(NA, 360))))) +
  
  # geom_rect(data=age.scale, aes(xmin=start, xmax=end, ymin=-5.5, ymax=5.5), colour="gray20", fill="black", alpha=0.05, linetype=2) +
  geom_rect(data=epoch.scale, aes(xmin=start, xmax=end, ymin=-9.5, ymax=-7.7), colour="black", fill="white") +
  geom_text(data=epoch.scale, aes(x=center, y=-8.6, label=name), colour="black", size=4) +
  
  geom_line(aes(x=years/1000, y=slope, colour=variable, group=inter), size=0.8) +
  
  scale_x_reverse(expand=c(0,0), name="", labels=NULL) +
  scale_y_continuous(expand=c(0,0), name="Slopes (model 3)") +
  scale_colour_discrete(name="") +
  facet_grid(pcs_metric ~ model_type) +
  expand_limits(y=c(-9.5, 5.5)) +
  theme_bw() 

d <- ggplot(subset(slopes, model == "m3" &
                     model_type %in% c("OLS", "SARerr") &
                     variable %in% VARS_NEW_NAMES[-6] &
                     ((pcs_metric == "NRI" & distance %in% c(NA, 480))|(pcs_metric == "NTI" & distance %in% c(NA, 360))))) +
  
  # geom_rect(data=age.scale, aes(xmin=start, xmax=end, ymin=-.1, ymax=.05), colour="gray20", fill="black", alpha=0.05, linetype=2) +
  geom_rect(data=epoch.scale, aes(xmin=start, xmax=end, ymin=-.1, ymax=-.08), colour="black", fill="white") +
  geom_text(data=epoch.scale, aes(x=center, y=-.09, label=name), colour="black", size=4) +
  
  geom_line(aes(x=years/1000, y=slope, colour=variable, group=inter), size=0.8) +
  
  scale_x_reverse(expand=c(0,0), name="Time (ka BP)") +
  scale_y_continuous(expand=c(0,0), name="Slopes (model 3)") +
  scale_colour_discrete(name="") +
  facet_grid(pcs_metric ~ model_type) +
  expand_limits(y = c(-.1, .05)) + 
  theme_bw() 

fig.slope <- ggarrange(c, d, ncol=1, nrow=2, align = "v", common.legend=TRUE, legend="top", labels=c("a)", "b)"))
ggexport(fig.slope, filename = "Results/Figures/Regression_through_time/Slopes.pdf", width=6.5, height=7.5)
ggexport(fig.slope, filename = "Results/Manuscript_figures_and_tables/Figure_4.pdf", width=6.5, height=7.5)

