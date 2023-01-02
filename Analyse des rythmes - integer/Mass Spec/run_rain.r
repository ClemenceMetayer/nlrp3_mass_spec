library(jsonlite)
require(rain)
library(readr)

nrep <- 4
period <- 24
deltat <- 4 # time interval between two CTs

data <- read_csv("~/Library/CloudStorage/OneDrive-INSTITUTCURIE/Doctorat/Analyse des rythmes - int/Mass Spec/data/ms_data_rain_ctrl.csv")

names <- read_csv("~/Library/CloudStorage/OneDrive-INSTITUTCURIE/Doctorat/Analyse des rythmes - int/Mass Spec/data/names_wcond_ctrl.csv")
names$`0`[names$`0`=="KPNA3__uid5876013_ctrl"]="KPNA3_ctrl"
names$`0`[names$`0`=="KPNA3__uid5876013_nlrp3"]="KPNA3_nlrp3"

results <- rain(t(data), deltat=deltat, period=period, nr.series=nrep,
                 peak.border=c(0.1, 0.9), verbose=TRUE, na.rm=TRUE) # last argument is to allow handling NaN values

results$names <- names$`0`
results_pval <- results[results$pVal <= 0.05,]
write_csv(results_pval, "~/Library/CloudStorage/OneDrive-INSTITUTCURIE/Doctorat/Analyse des rythmes - int/Mass Spec/Results/rain_results.csv")



