library(ape)
library(raster)
library(spdep)
library(spatialreg)

# Arrange data ------------------------------------------------------------

# Subsetting to periods with more than 10 data
subsample_data <- function(p, d, smin){
  nd <- subset(d, year==p)
  if(nrow(nd) > smin){return(nd)}
}

data.cels.subset <- pcs.data
data.cels.subset <- data.cels.subset[, c("nri", "nri_p", "nri_significant",  "nti", "nti_p", "nti_significant", "year", "Tmin", "Tmax", "Pmin", "Pmax", "AET", "ETR", "WDI", "Deglac.", "x", "y")]
data.cels.subset <- lapply(paste(PERIODS, "BP", sep=""), FUN=subsample_data, data.cels.subset, 10)
data.cels.subset <- do.call(rbind, data.cels.subset)



# Neighbour definition for Spatial Autocorrelation tests ------------------

# define neighbors with distance d
# damax=1
neigh.d1 <- dnearneigh(cbind(data.cels.subset$x, data.cels.subset$y), d1=0, d2=120, longlat=TRUE)
w.d1 <- nb2listw(neigh.d1, style = "W", zero.policy =TRUE)

# damax=2
neigh.d2 <- dnearneigh(cbind(data.cels.subset$x, data.cels.subset$y), d1=0, d2=360, longlat=TRUE)
w.d2 <- nb2listw(neigh.d2, style = "W", zero.policy = TRUE)

# damax=3
neigh.d3 <- dnearneigh(cbind(data.cels.subset$x, data.cels.subset$y), d1=0, d2=480, longlat=TRUE)
w.d3 <- nb2listw(neigh.d3, style = "W", zero.policy = TRUE)

rm(neigh.d1, neigh.d2, neigh.d3)



# Raw Spatial Autocorrelation tests ---------------------------------------

get_raw_moran <- function(moran.test){
  output <- data.frame(moran.i=double(), p.value=double(), expectation=double(),
                       variance=double(), std_dev=double(),
                       alternative=as.character(), stringsAsFactors=FALSE)
  output[1, "moran.i"] <- moran.test$estimate[1]
  output[1, "p.value"] <- moran.test$p.value
  output[1, "std_dev"] <- moran.test$statistic
  output[1, "expectation"] <- moran.test$estimate[2]
  output[1, "variance"] <- moran.test$estimate[3]
  output[1, "alternative"] <- moran.test$alternative
  
  return(output)
}

# test for spatial autoregression in raw NRI and NTi
pcs.moran <- list() # list for pcs metrics
pcs.moran$nri <- list() # list for spatial autocorrelation distances
pcs.moran$nri$"120" <- moran.test(data.cels.subset$nri, listw=w.d1, zero.policy= TRUE, na.action=na.exclude)
pcs.moran$nri$"360" <- moran.test(data.cels.subset$nri, listw=w.d2, zero.policy= TRUE, na.action=na.exclude)
pcs.moran$nri$"480" <- moran.test(data.cels.subset$nri, listw=w.d3, zero.policy= TRUE, na.action=na.exclude)

pcs.moran$nti <- list() # list for spatial autocorrelation distances
pcs.moran$nti$"120" <- moran.test(data.cels.subset$nti, listw=w.d1, zero.policy= TRUE, na.action=na.exclude)
pcs.moran$nti$"360" <- moran.test(data.cels.subset$nti, listw=w.d2, zero.policy= TRUE, na.action=na.exclude)
pcs.moran$nti$"480" <- moran.test(data.cels.subset$nti, listw=w.d3, zero.policy= TRUE, na.action=na.exclude)

# extract parameters
pcs.moran$nri <- lapply(pcs.moran$nri, FUN=get_raw_moran)
pcs.moran$nti <- lapply(pcs.moran$nti, FUN=get_raw_moran)

# melt parameters
pcs.moran <- melt(pcs.moran, id.vars=colnames(pcs.moran$nri$"120"))
colnames(pcs.moran) <- c("moran.i", "p.value", "expectation", "variance", "std_dev", "alternative", "distance", "pcs_metrics")

# write to the hard disk drive
write.csv(pcs.moran, file="Results/Tables/RAW_moran.csv")



# Fit OLS models ----------------------------------------------------------

fit_ols <- function(var, pcs_metric, data){

  formula.1 <- formula(paste0(pcs_metric, "~", var))
  formula.2 <- formula(paste0(pcs_metric, "~", var, " + year"))
  formula.3 <- formula(paste0(pcs_metric, "~", var, " * year"))
  
  ols.1 <- eval(bquote(lm(.(formula.1), data)))
  ols.2 <- eval(bquote(lm(.(formula.2), data)))
  ols.3 <- eval(bquote(lm(.(formula.3), data)))
  
  ols <- list(ols.1, ols.2, ols.3)
  names(ols) <- c("m1", "m2", "m3")
  return(ols)
}

get_ols_param <- function(ols.list){
  require(broom)
  
  m1.output <- cbind(coef(ols.list$m1)[1], coef(ols.list$m1)[2], glance(ols.list$m1))
  m2.output <- cbind(coef(ols.list$m2)[2], glance(ols.list$m2))
  m3.output <- glance(ols.list$m3)
  
  col.names <- colnames(m3.output)
  colnames(m1.output) <- paste0("m1.", c("intercept", "slope", col.names))
  colnames(m2.output) <- paste0("m2.", c("slope", col.names))
  colnames(m3.output) <- paste0("m3.", col.names)
  
  output <- cbind(m1.output, m2.output, m3.output)
  return(output)
}

pcs.ols <- list()
pcs.ols$nri <- lapply(c(COORD_COLNAMES, VARS_NEW_NAMES), FUN=fit_ols, "nri", data.cels.subset)
names(pcs.ols$nri) <- c(COORD_COLNAMES, VARS_NEW_NAMES)
pcs.ols$nti <- lapply(c(COORD_COLNAMES, VARS_NEW_NAMES), FUN=fit_ols, "nti", data.cels.subset)
names(pcs.ols$nti) <- c(COORD_COLNAMES, VARS_NEW_NAMES)

pcs.ols.param <- list()
pcs.ols.param$nri <- lapply(pcs.ols$nri, FUN=get_ols_param)
pcs.ols.param$nti <- lapply(pcs.ols$nti, FUN=get_ols_param)

pcs.ols.param <- melt(pcs.ols.param, id.vars=colnames(pcs.ols.param$nri[[1]]))
colnames(pcs.ols.param) <- c(colnames(pcs.ols.param[c(1:36)]), "variable", "pcs_metric")

write.csv(pcs.ols.param, file="Results/Tables/OLS_models.csv")
rm(pcs.ols.param)


# OLS Spatial Autocorrelation tests ---------------------------------------

get_ols_moran <- function(models, nb.weight) {
  
  require (spdep)
  
  mod.1 <- models$m1
  mod.2 <- models$m2
  mod.3 <- models$m3
  
  output <- data.frame(moran_i_1=double(),
                       p.value_1=double(),
                       moran_i_2=double(),
                       p.value_2=double(),
                       moran_i_3=double(),
                       p.value_3=double())
  
  moran.1 <- lm.morantest(mod.1, nb.weight, zero.policy= TRUE)
  output[1,1] <- moran.1$estimate[1]
  output[1,2] <- moran.1$p.value
  
  moran.2 <- lm.morantest(mod.2, nb.weight, zero.policy= TRUE)
  output[1,3] <- moran.2$estimate[1]
  output[1,4] <- moran.2$p.value
  
  moran.3 <- lm.morantest(mod.3, nb.weight, zero.policy= TRUE)
  output[1,5] <- moran.3$estimate[1]
  output[1,6] <- moran.3$p.value
  
  return(output)
}

# test for spatial autoregression in NRI and NTI residuals from OLS
pcs.ols.moran <- list()
pcs.ols.moran$nri <- list()
pcs.ols.moran$nri$"120" <- lapply(pcs.ols$nri, FUN=get_ols_moran, w.d1)
pcs.ols.moran$nri$"360" <- lapply(pcs.ols$nri, FUN=get_ols_moran, w.d2)
pcs.ols.moran$nri$"480" <- lapply(pcs.ols$nri, FUN=get_ols_moran, w.d3)

pcs.ols.moran$nti <- list()
pcs.ols.moran$nti$"120" <- lapply(pcs.ols$nti, FUN=get_ols_moran, w.d1)
pcs.ols.moran$nti$"360" <- lapply(pcs.ols$nti, FUN=get_ols_moran, w.d2)
pcs.ols.moran$nti$"480" <- lapply(pcs.ols$nti, FUN=get_ols_moran, w.d3)

# melt results
pcs.ols.moran <- melt(pcs.ols.moran, id.vars=c("moran_i_1", "p.value_1", "moran_i_2", "p.value_2", "moran_i_3", "p.value_3"))
colnames(pcs.ols.moran) <- c("moran_i_1", "p.value_1", "moran_i_2", "p.value_2", "moran_i_3", "p.value_3", "variable", "distance", "pcs_metric")

write.csv(pcs.ols.moran, "Results/Tables/OLS_moran.csv")



# Fit SAR models ----------------------------------------------------------

fit_sar <- function(var, pcs_metric, data, w, sar.type=c("lag", "err")){
  
  formula.1 <- formula(paste0(pcs_metric, "~", var))
  formula.2 <- formula(paste0(pcs_metric, "~", var, " + year"))
  formula.3 <- formula(paste0(pcs_metric, "~", var, " * year"))
  
  if(sar.type == "lag"){
    sar.1 <- eval(bquote(lagsarlm(.(formula.1), data, listw = w, zero.policy= TRUE, na.action=na.exclude)))
    sar.2 <- eval(bquote(lagsarlm(.(formula.2), data, listw = w, zero.policy= TRUE, na.action=na.exclude)))
    sar.3 <- eval(bquote(lagsarlm(.(formula.3), data, listw = w, zero.policy= TRUE, na.action=na.exclude)))
  }
  if(sar.type == "err"){
    sar.1 <- eval(bquote(errorsarlm(.(formula.1), data, listw = w, zero.policy= TRUE, na.action=na.exclude)))
    sar.2 <- eval(bquote(errorsarlm(.(formula.2), data, listw = w, zero.policy= TRUE, na.action=na.exclude)))
    sar.3 <- eval(bquote(errorsarlm(.(formula.3), data, listw = w, zero.policy= TRUE, na.action=na.exclude)))
  }  
  sar <- list(sar.1, sar.2, sar.3)
  names(sar) <- c("m1", "m2", "m3")
  return(sar)
}


# Fit SAR-lag models
pcs.sar.lag <- list()
pcs.sar.lag$nri <- list()
pcs.sar.lag$nri$"120" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nri", data.cels.subset, w.d1, "lag")
pcs.sar.lag$nri$"360" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nri", data.cels.subset, w.d2, "lag")
pcs.sar.lag$nri$"480" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nri", data.cels.subset, w.d3, "lag")

pcs.sar.lag$nti <- list()
pcs.sar.lag$nti$"120" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nti", data.cels.subset, w.d1, "lag")
pcs.sar.lag$nti$"360" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nti", data.cels.subset, w.d2, "lag")
pcs.sar.lag$nti$"480" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nti", data.cels.subset, w.d3, "lag")

names(pcs.sar.lag$nri$"120") <- names(pcs.sar.lag$nri$"360") <- names(pcs.sar.lag$nri$"480") <- VARS_NEW_NAMES
names(pcs.sar.lag$nti$"120") <- names(pcs.sar.lag$nti$"360") <- names(pcs.sar.lag$nti$"480") <- VARS_NEW_NAMES

# save(pcs.sar.lag, file="Results/RObjects/SAR_lag_models.RData")
# load("Results/RObjects/SAR_lag_models.RData")


# Fit SAR-err models
pcs.sar.err <- list()
pcs.sar.err$nri <- list()
pcs.sar.err$nri$"120" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nri", data.cels.subset, w.d1, "err")
pcs.sar.err$nri$"360" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nri", data.cels.subset, w.d2, "err")
pcs.sar.err$nri$"480" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nri", data.cels.subset, w.d3, "err")

pcs.sar.err$nti <- list()
pcs.sar.err$nti$"120" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nti", data.cels.subset, w.d1, "err")
pcs.sar.err$nti$"360" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nti", data.cels.subset, w.d2, "err")
pcs.sar.err$nti$"480" <- lapply(VARS_NEW_NAMES, FUN=fit_sar, "nti", data.cels.subset, w.d3, "err")

names(pcs.sar.err$nri$"120") <- names(pcs.sar.err$nri$"360") <- names(pcs.sar.err$nri$"480") <- VARS_NEW_NAMES
names(pcs.sar.err$nti$"120") <- names(pcs.sar.err$nti$"360") <- names(pcs.sar.err$nti$"480") <- VARS_NEW_NAMES

# save(pcs.sar.err, file="Results/RObjects/SAR_err_models.RData")
# load("Results/RObjects/SAR_err_models.RData")


# Get parameters from the SAR models
get_sar_param <- function(sar.list){
  require(spdep)
  .get_sar_param <- function(model) {
    require(spdep)
    output <- data.frame(intercept=double(), slope=double(),
                         r.square=double(), p.value=double(),
                         lr.value=double(), lr.p.value=double(),
                         aic=double(), lm.aic=double(), 
                         lm.value=double())
    
    output[1, "intercept"] <- summary(model, Nagelkerke = T)$coefficients[1]
    output[1, "slope"] <- summary(model, Nagelkerke = T)$coefficients[2]
    output[1, "r.square"] <- summary(model, Nagelkerke = T)$NK
    output[1, "p.value"] <- summary(model, Nagelkerke = T)$Coef[2,4]
    output[1, "lr.value"] <- summary(model, Nagelkerke = T)$LR1$statistic[[1]]
    output[1, "lr.p.value"] <- summary(model, Nagelkerke = T)$LR1$p.value[[1]]
    output[1, "aic"] <- AIC(model)
    output[1, "lm.aic"] <- summary(model, Nagelkerke = T)$AIC_lm.model
    if(is.null(summary(model, Nagelkerke=T)$LMtest)){
      output[1, "lm.value"] <- NA
    }else{
      output[1, "lm.value"] <- summary(model, Nagelkerke = T)$LMtest
    }
    return(output)
  }
  
  output.1 <- .get_sar_param(sar.list$m1)
  output.2 <- .get_sar_param(sar.list$m2)
  output.3 <- .get_sar_param(sar.list$m3)
  output.colnames <- colnames(output.1)
  output <- cbind(output.1, output.2, output.3)
  colnames(output) <- c(paste0("m1.", output.colnames), paste0("m2.", output.colnames), paste0("m3.", output.colnames))
  return(output)
}

pcs.sar.lag.param <- list()

pcs.sar.lag.param$nri <- list()
pcs.sar.lag.param$nri$"120" <- lapply(pcs.sar.lag$nri$"120", FUN=get_sar_param)
pcs.sar.lag.param$nri$"360" <- lapply(pcs.sar.lag$nri$"360", FUN=get_sar_param)
pcs.sar.lag.param$nri$"480" <- lapply(pcs.sar.lag$nri$"480", FUN=get_sar_param)

pcs.sar.lag.param$nti <- list()
pcs.sar.lag.param$nti$"120" <- lapply(pcs.sar.lag$nti$"120", FUN=get_sar_param)
pcs.sar.lag.param$nti$"360" <- lapply(pcs.sar.lag$nti$"360", FUN=get_sar_param)
pcs.sar.lag.param$nti$"480" <- lapply(pcs.sar.lag$nti$"480", FUN=get_sar_param)

pcs.sar.lag.param <- melt(pcs.sar.lag.param, id.vars=colnames(pcs.sar.lag.param$nri$"120"[[1]]))
colnames(pcs.sar.lag.param) <- c(colnames(pcs.sar.lag.param)[c(1:27)], "variable", "distance", "pcs_metric")

write.csv(pcs.sar.lag.param, file="Results/Tables/SAR_lag_parameters.csv")
rm(pcs.sar.lag.param)


pcs.sar.err.param <- list()

pcs.sar.err.param$nri <- list()
pcs.sar.err.param$nri$"120" <- lapply(pcs.sar.err$nri$"120", FUN=get_sar_param)
pcs.sar.err.param$nri$"360" <- lapply(pcs.sar.err$nri$"360", FUN=get_sar_param)
pcs.sar.err.param$nri$"480" <- lapply(pcs.sar.err$nri$"480", FUN=get_sar_param)

pcs.sar.err.param$nti <- list()
pcs.sar.err.param$nti$"120" <- lapply(pcs.sar.err$nti$"120", FUN=get_sar_param)
pcs.sar.err.param$nti$"360" <- lapply(pcs.sar.err$nti$"360", FUN=get_sar_param)
pcs.sar.err.param$nti$"480" <- lapply(pcs.sar.err$nti$"480", FUN=get_sar_param)

pcs.sar.err.param <- melt(pcs.sar.err.param, id.vars=colnames(pcs.sar.err.param$nri$"120"[[1]]))
colnames(pcs.sar.err.param) <- c(colnames(pcs.sar.err.param)[c(1:27)], "variable", "distance", "pcs_metric")

write.csv(pcs.sar.err.param, file="Results/Tables/SAR_err_parameters.csv")
rm(pcs.sar.err.param)



# SAR Spatial Autocorrelation tests ---------------------------------------

get_sar_moran <- function(sar.list, nb.weight){
  .get_sar_moran <- function(models, nb.weight) {
    
    require (spdep)
    
    mod.1 <- models[[1]]
    mod.2 <- models[[2]]
    mod.3 <- models[[3]]
    
    output <- data.frame(moran_i_1=double(),
                         p.value_1=double(),
                         moran_i_2=double(),
                         p.value_2=double(),
                         moran_i_3=double(),
                         p.value_3=double())
    
    moran.1 <- moran.test(residuals(mod.1),
                          listw=nb.weight,
                          zero.policy= TRUE,
                          na.action=na.exclude)
    output[1,1] <- moran.1$estimate[1]
    output[1,2] <- moran.1$p.value
    
    moran.2 <- moran.test(residuals(mod.2),
                          listw = nb.weight,
                          zero.policy = TRUE,
                          na.action = na.exclude)
    output[1,3] <- moran.2$estimate[1]
    output[1,4] <- moran.2$p.value
    
    moran.3 <- moran.test(residuals(mod.3),
                          listw = nb.weight,
                          zero.policy = TRUE,
                          na.action = na.exclude)
    output[1,5] <- moran.3$estimate[1]
    output[1,6] <- moran.3$p.value
    
    return(output)
  }
  
  results <- lapply(sar.list, FUN=.get_sar_moran, nb.weight)
  names(results) <- names(sar.list)
  return(results)
}

pcs.sar.lag.moran <- list()
pcs.sar.lag.moran$nri <- list()
pcs.sar.lag.moran$nri$"120" <- lapply(pcs.sar.lag$nri, FUN=get_sar_moran, w.d1)
pcs.sar.lag.moran$nri$"360" <- lapply(pcs.sar.lag$nri, FUN=get_sar_moran, w.d2)
pcs.sar.lag.moran$nri$"480" <- lapply(pcs.sar.lag$nri, FUN=get_sar_moran, w.d3)

pcs.sar.lag.moran$nti <- list()
pcs.sar.lag.moran$nti$"120" <- lapply(pcs.sar.lag$nti, FUN=get_sar_moran, w.d1)
pcs.sar.lag.moran$nti$"360" <- lapply(pcs.sar.lag$nti, FUN=get_sar_moran, w.d2)
pcs.sar.lag.moran$nti$"480" <- lapply(pcs.sar.lag$nti, FUN=get_sar_moran, w.d3)

pcs.sar.err.moran <- list()
pcs.sar.err.moran$nri <- list()
pcs.sar.err.moran$nri$"120" <- lapply(pcs.sar.err$nri, FUN=get_sar_moran, w.d1)
pcs.sar.err.moran$nri$"360" <- lapply(pcs.sar.err$nri, FUN=get_sar_moran, w.d2)
pcs.sar.err.moran$nri$"480" <- lapply(pcs.sar.err$nri, FUN=get_sar_moran, w.d3)

pcs.sar.err.moran$nti <- list()
pcs.sar.err.moran$nti$"120" <- lapply(pcs.sar.err$nti, FUN=get_sar_moran, w.d1)
pcs.sar.err.moran$nti$"360" <- lapply(pcs.sar.err$nti, FUN=get_sar_moran, w.d2)
pcs.sar.err.moran$nti$"480" <- lapply(pcs.sar.err$nti, FUN=get_sar_moran, w.d3)

pcs.sar.lag.moran <- melt(pcs.sar.lag.moran, id.vars=colnames(pcs.sar.lag.moran$nri$"120"$"120"$Tmin))
pcs.sar.err.moran <- melt(pcs.sar.err.moran, id.vars=colnames(pcs.sar.err.moran$nri$"120"$"120"$Tmin))

colnames(pcs.sar.lag.moran) <- colnames(pcs.sar.err.moran) <- c("moran_i_1", "p.value_1", "moran_i_2", "p.value_2", "moran_i_3", "p.value_3", "variable", "test_distance", "model_distance", "pcs_metric")

write.csv(pcs.sar.lag.moran, file = "Results/Tables/SAR_lag_moran.csv")
write.csv(pcs.sar.err.moran, file = "Results/Tables/SAR_err_moran.csv")



# OLS - ANOVA -------------------------------------------------------------

get_ols_anova <- function(model_list){
  aov <- anova(model_list$m1, model_list$m2, model_list$m3)
  
  output <- data.frame(m1.Res.Df=double(),
                       m2.Res.Df=double(),
                       m3.Res.Df=double(),
                       m1.RSS=double(),
                       m2.RSS=double(),
                       m3.RSS=double(),
                       m1.Df=double(),
                       m2.Df=double(),
                       m3.Df=double(),
                       m1.Sum_of_sq=double(),
                       m2.Sum_of_sq=double(),
                       m3.Sum_of_sq=double(), 
                       m1.F=double(), 
                       m2.F=double(),
                       m3.F=double(),
                       m1.p.value=double(), 
                       m2.p.value=double(),
                       m3.p.value=double())
  output[1, 1:3] <- aov$Res.Df
  output[1, 4:6] <- aov$RSS
  output[1, 7:9] <- aov$Df
  output[1, 10:12] <- aov$`Sum of Sq`
  output[1, 13:15] <- aov$F
  output[1, 16:18] <- aov$`Pr(>F)`
  return(output)
}

pcs.ols.aov <- list()
pcs.ols.aov$nri <- lapply(pcs.ols$nri, FUN=get_ols_anova)
pcs.ols.aov$nti <- lapply(pcs.ols$nti, FUN=get_ols_anova)

names(pcs.ols.aov) <- PCS_METRICS

pcs.ols.aov <- melt(pcs.ols.aov, id.vars=colnames(pcs.ols.aov$nri$Tmin))
colnames(pcs.ols.aov) <- c(colnames(pcs.ols.aov)[1:18], "variable", "pcs_metric")

write.csv(pcs.ols.aov, file="Results/Tables/OLS_anova.csv")



# SAR - ANOVA ---------------------------------------------------------

get_sar_anova <- function(model_list, variable, metric){

  aov <- anova(model_list$m1, model_list$m2, model_list$m3)
  
  output <- data.frame(m1.df=integer(), m2.df=integer(), m3.df=integer(),
                       m1.AIC=double(), m2.AIC=double(), m3.AIC=double(),
                       m1.logLik=double(), m2.logLik=double(), m3.logLik=double(),
                       m1.Test=character(), m2.Test=character(), m3.Test=character(), 
                       m1.L.Ratio=double(), m2.L.Ratio=double(), m3.L.Ratio=double(),
                       m1.p.value=double(), m2.p.value=double(), m3.p.value=double(),
                       stringsAsFactors=FALSE)
  
  output[1, paste0("m", 1:3, ".df")] <- aov$df
  output[1, paste0("m", 1:3, ".AIC")] <- aov$AIC
  output[1, paste0("m", 1:3, ".logLik")] <- aov$logLik
  output[1, paste0("m", 1:3, ".Test")] <- aov$Test
  output[1, paste0("m", 1:3, ".L.Ratio")] <- aov$L.Ratio
  output[1, paste0("m", 1:3, ".p.value")] <- aov$`p-value`

  return(output)  
}


# SAR lag

pcs.sar.lag.aov <- list()

pcs.sar.lag.aov$nri <- list()
pcs.sar.lag.aov$nri$"120" <- mapply(FUN = get_sar_anova, pcs.sar.lag$nri$"120", names(pcs.sar.lag$nri$"120"), MoreArgs = list(metric="nri"), SIMPLIFY=FALSE)
pcs.sar.lag.aov$nri$"360" <- mapply(FUN = get_sar_anova, pcs.sar.lag$nri$"360", names(pcs.sar.lag$nri$"360"), MoreArgs = list(metric="nri"), SIMPLIFY=FALSE)
pcs.sar.lag.aov$nri$"480" <- mapply(FUN = get_sar_anova, pcs.sar.lag$nri$"480", names(pcs.sar.lag$nri$"480"), MoreArgs = list(metric="nri"), SIMPLIFY=FALSE)

pcs.sar.lag.aov$nti <- list()
pcs.sar.lag.aov$nti$"120" <- mapply(FUN = get_sar_anova, pcs.sar.lag$nti$"120", names(pcs.sar.lag$nti)[[1]], MoreArgs = list(names(pcs.sar.lag)[[2]]), SIMPLIFY=FALSE)
pcs.sar.lag.aov$nti$"360" <- mapply(FUN = get_sar_anova, pcs.sar.lag$nti$"360", names(pcs.sar.lag$nti)[[2]], MoreArgs = list(names(pcs.sar.lag)[[2]]), SIMPLIFY=FALSE)
pcs.sar.lag.aov$nti$"480" <- mapply(FUN = get_sar_anova, pcs.sar.lag$nti$"480", names(pcs.sar.lag$nti)[[3]], MoreArgs = list(names(pcs.sar.lag)[[2]]), SIMPLIFY=FALSE)


pcs.sar.lag.aov <- melt(pcs.sar.lag.aov, id.vars=colnames(pcs.sar.lag.aov$nri$"120"[[1]]))
colnames(pcs.sar.lag.aov) <- c(colnames(pcs.sar.lag.aov)[1:18], "variable", "distance", "pcs_metric")

write.csv(pcs.sar.lag.aov, file="Results/Tables/SAR_lag_anova.csv")
rm(pcs.sar.lag.aov)


# SAR err

pcs.sar.err.aov <- list()

pcs.sar.err.aov$nri <- list()
pcs.sar.err.aov$nri$"120" <- mapply(FUN = get_sar_anova, pcs.sar.err$nri$"120", names(pcs.sar.err$nri)[[1]], MoreArgs = list(names(pcs.sar.err)[[1]]), SIMPLIFY=FALSE)
pcs.sar.err.aov$nri$"360" <- mapply(FUN = get_sar_anova, pcs.sar.err$nri$"360", names(pcs.sar.err$nri)[[2]], MoreArgs = list(names(pcs.sar.err)[[1]]), SIMPLIFY=FALSE)
pcs.sar.err.aov$nri$"480" <- mapply(FUN = get_sar_anova, pcs.sar.err$nri$"480", names(pcs.sar.err$nri)[[3]], MoreArgs = list(names(pcs.sar.err)[[1]]), SIMPLIFY=FALSE)

pcs.sar.err.aov$nti <- list()
pcs.sar.err.aov$nti$"120" <- mapply(FUN = get_sar_anova, pcs.sar.err$nti$"120", names(pcs.sar.err$nti)[[1]], MoreArgs = list(names(pcs.sar.err)[[2]]), SIMPLIFY=FALSE)
pcs.sar.err.aov$nti$"360" <- mapply(FUN = get_sar_anova, pcs.sar.err$nti$"360", names(pcs.sar.err$nti)[[2]], MoreArgs = list(names(pcs.sar.err)[[2]]), SIMPLIFY=FALSE)
pcs.sar.err.aov$nti$"480" <- mapply(FUN = get_sar_anova, pcs.sar.err$nti$"480", names(pcs.sar.err$nti)[[3]], MoreArgs = list(names(pcs.sar.err)[[2]]), SIMPLIFY=FALSE)


pcs.sar.err.aov <- melt(pcs.sar.err.aov, id.vars=colnames(pcs.sar.err.aov$nri$"120"[[1]]))
colnames(pcs.sar.err.aov) <- c(colnames(pcs.sar.err.aov)[1:18], "variable", "distance", "pcs_metric")

write.csv(pcs.sar.err.aov, file="Results/Tables/SAR_err_anova.csv")
rm(pcs.sar.err.aov)

