library(reshape2)

change_p_values <- function(p_values){
  output <- NULL
  output[which(p_values > 0.05)] <- "ns"
  output[which(p_values <= 0.05)] <- "*"
  output[which(p_values <= 0.01)] <- "**"
  output[which(p_values <= 0.001)] <- "***"
  return(output)
}

format_values <- function(values, p_values, dec = 3){
  
  values <- round(values, dec)
  
  p_values <- change_p_values(p_values)
  
  output <- paste0(values, p_values)
  return(output)
}



# Table 1
ols.columns <- c("m1.intercept", "m1.slope",
                      "m1.r.squared", "m1.p.value",
                      "m2.slope", "m2.r.squared", 
                      "m2.p.value", "m3.r.squared",
                      "m3.p.value")
ols.models <- read.csv("Results/Tables/OLS_models.csv", row.names=1)
ols.models <- ols.models[,c("pcs_metric", "variable", ols.columns)]
ols.models[,ols.columns] <- apply(ols.models[,ols.columns], FUN=round, MARGIN=2, digits=3)
ols.models[,c("m1.p.value", "m2.p.value", "m3.p.value")] <- apply(ols.models[,c("m1.p.value", "m2.p.value", "m3.p.value")],
                                                                  FUN = change_p_values, MARGIN = 2)
ols.models$distance <- NA
ols.models$model_type <- "OLS"

sar.columns <- c("m1.intercept", "m1.slope",
                 "m1.r.square", "m1.p.value",
                 "m2.slope", "m2.r.square", 
                 "m2.p.value", "m3.r.square",
                 "m3.p.value")
sar.models <- read.csv("Results/Tables/SAR_err_parameters.csv", row.names = 1)
sar.models <- sar.models[,c("pcs_metric", "distance", "variable", sar.columns)]
sar.models[, sar.columns] <- apply(sar.models[, sar.columns], FUN = round, MARGIN = 2, digits = 3)
sar.models[, c("m1.p.value", "m2.p.value", "m3.p.value")] <- apply(sar.models[, c("m1.p.value", "m2.p.value", "m3.p.value")],
                                                                   FUN = change_p_values, MARGIN = 2)
sar.models <- subset(sar.models, (pcs_metric == "nri" & distance == 480) | (pcs_metric == "nti" & distance == 360))
colnames(sar.models) <- c("pcs_metric", "distance", "variable", ols.columns)
sar.models$model_type <- "SARerr"

all.models <- rbind(ols.models, sar.models)
all.models$model <- paste0(all.models$model_type, all.models$distance)
all.models <- all.models[,c("pcs_metric", "model", "variable", ols.columns)]
all.models$variable <- factor(all.models$variable, levels=c("x", "y", "Tmin", "Tmax",
                                                               "Pmin", "Pmax", "AET", "ETR",
                                                               "WDI", "Deglac."))
all.models$pcs_metric <- factor(all.models$pcs_metric, levels = c("nri", "nti"))
write.csv2(all.models, file = "Results/Manuscript_figures_and_tables/Table_1.csv")

rm(ols.models, ols.columns, sar.models, sar.columns, all.models)

# Table 2
ols.columns <- c("m1.Res.Df","m1.F","m2.Res.Df","m2.F","m3.Res.Df","m3.F","m2.p.value", "m3.p.value")
ols.aov <- read.csv("Results/Tables/OLS_anova.csv", row.names = 1)
ols.aov$model_type <- "OLS"
ols.aov <- subset(ols.aov, select=c("pcs_metric", "model_type", "variable", ols.columns))
ols.aov$m2.F <- round(ols.aov$m2.F, digits = 3)
ols.aov$m3.F <- round(ols.aov$m3.F, digits = 3)
ols.aov$m2.p.value <- change_p_values(ols.aov$m2.p.value)
ols.aov$m3.p.value <- change_p_values(ols.aov$m3.p.value)
ols.aov$selected_mod <- paste(paste0("2", ols.aov$m2.p.value), paste0("3", ols.aov$m3.p.value), sep=", ")

sar.columns <- c("m1.df", "m1.L.Ratio", "m2.df", "m2.L.Ratio", "m3.df", "m3.L.Ratio", "m2.p.value", "m3.p.value")
sar.aov <- read.csv("Results/Tables/SAR_err_anova.csv", row.names = 1)
sar.aov <- subset(sar.aov, (pcs_metric == "nri" & distance == 480) | (pcs_metric == "nti" & distance == 360))
sar.aov <- subset(sar.aov, select=c("pcs_metric", "distance", "variable", sar.columns))
sar.aov$m2.L.Ratio <- round(sar.aov$m2.L.Ratio, digits = 3)
sar.aov$m3.L.Ratio <- round(sar.aov$m3.L.Ratio, digits = 3)
sar.aov$m2.p.value <- change_p_values(sar.aov$m2.p.value)
sar.aov$m3.p.value <- change_p_values(sar.aov$m3.p.value)
sar.aov$model_type <- paste0("SARerr", sar.aov$distance)
sar.aov$selected_mod <- paste(paste0("2", sar.aov$m2.p.value), paste0("3", sar.aov$m3.p.value), sep=", ")
sar.aov <- sar.aov[, c("pcs_metric", "model_type", "variable", sar.columns, "selected_mod")]

write.csv2(ols.aov, file = "Results/Manuscript_figures_and_tables/Table2_1.csv")
write.csv2(sar.aov, file = "Results/Manuscript_figures_and_tables/Table2_2.csv")

# Supplement A2.1
raw.moran <- read.csv("Results/Tables/RAW_moran.csv", row.names=1)
raw.moran$format_values <- format_values(raw.moran$moran.i, raw.moran$p.value)
raw.moran <- raw.moran[,c("pcs_metrics", "distance", "format_values")]
colnames(raw.moran) <- c("PCS metric", "Distance", "Moran's I")
write.csv2(raw.moran, file = "Results/Manuscript_figures_and_tables/Supplement_A2_1.csv")
rm(raw.moran)

# Supplement A2.2
ols.moran <- read.csv("Results/Tables/OLS_moran.csv", row.names=1)
ols.moran$format_val_1 <- format_values(ols.moran$moran_i_1, ols.moran$p.value_1)
ols.moran$format_val_2 <- format_values(ols.moran$moran_i_2, ols.moran$p.value_2)
ols.moran$format_val_3 <- format_values(ols.moran$moran_i_3, ols.moran$p.value_3)
ols.moran <- ols.moran[,c("pcs_metric", "distance", "variable", "format_val_1", "format_val_2", "format_val_3")]
colnames(ols.moran) <- c("PCS metric", "Distance", "Var.", "Model 1", "Model 2", "Model 3")
write.csv2(ols.moran, file = "Results/Manuscript_figures_and_tables/Supplement_A2_2.csv")
rm(ols.moran)

# Supplement A2.3
sar.moran <- read.csv("Results/Tables/SAR_err_moran.csv", row.names=1)
sar.moran$format_val_1 <- format_values(sar.moran$moran_i_1, sar.moran$p.value_1)
sar.moran$format_val_2 <- format_values(sar.moran$moran_i_2, sar.moran$p.value_2)
sar.moran$format_val_3 <- format_values(sar.moran$moran_i_3, sar.moran$p.value_3)
sar.moran <- sar.moran[, c("pcs_metric", "test_distance", "model_distance", "variable", "format_val_1", "format_val_2", "format_val_3")]
colnames(sar.moran) <- c("PCS metric", "Moran dist.", "SAR dist.", "Var.", "Model 1", "Model 2", "Model 3")
write.csv2(sar.moran, file = "Results/Manuscript_figures_and_tables/Supplement_A2_3.csv")
rm(sar.moran)

# Supplement A2.4
ols.models <- read.csv("Results/Tables/OLS_models.csv", row.names=1)
ols.models <- ols.models[,c("pcs_metric", "variable", "m1.AIC", "m2.AIC", "m3.AIC")]
ols.models$m1.AIC <- round(ols.models$m1.AIC, 1)
ols.models$m2.AIC <- round(ols.models$m2.AIC, 1)
ols.models$m3.AIC <- round(ols.models$m3.AIC, 1)
ols.models <- melt(ols.models, id.vars=c("pcs_metric", "variable"))
colnames(ols.models) <- c("pcs_metric", "variable", "model", "aic")
ols.models$model_type <- "OLS"
ols.models$distance <- NA

sar.models <- read.csv("Results/Tables/SAR_err_parameters.csv", row.names=1)
sar.models <- sar.models[,c("pcs_metric", "distance", "variable", "m1.aic", "m2.aic", "m3.aic")]
sar.models$m1.AIC <- round(sar.models$m1.aic, 1)
sar.models$m2.AIC <- round(sar.models$m2.aic, 1)
sar.models$m3.AIC <- round(sar.models$m3.aic, 1)
sar.models <- sar.models[,c("pcs_metric", "distance", "variable", "m1.AIC", "m2.AIC", "m3.AIC")]
sar.models <- melt(sar.models, id.vars=c("pcs_metric", "distance", "variable"))
colnames(sar.models) <- c("pcs_metric", "distance", "variable", "model", "aic")
sar.models$model_type <- "SAR"

aic.models <- rbind(ols.models, sar.models)
aic.models <- acast(aic.models, model + variable ~ pcs_metric + model_type + distance, value.var="aic")
write.csv2(aic.models, file = "Results/Manuscript_figures_and_tables/Supplement_A2_4.csv")
rm(ols.models, sar.models, aic.models)
