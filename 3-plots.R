# Process Picarro data for DWP lab experiment
# Ben Bond-Lamberty April 2015

source("0-functions.R")

SCRIPTNAME  	<- "3-plots.R"
SUMMARYDATA      <- paste0(OUTPUT_DIR, "summarydata.csv")  # output from script 2


# ==============================================================================
# Main 

sink(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), split=T) # open log

printlog("Welcome to", SCRIPTNAME)

printlog("Reading in summary data...")
summarydata <- read_csv(SUMMARYDATA)
summarydata$DATETIME <- ymd_hms(summarydata$DATETIME)
print_dims(summarydata)
print(summary(summarydata))

p1 <- ggplot(summarydata, aes(DATETIME, CO2_dry, color=DWP)) + geom_point()
p1 <- p1 + facet_grid(Source~Rep, scales="free")
p1 <- p1 + theme(axis.text.x = element_text(angle = 90))
print(p1)
save_plot("CO2_source_rep")
p2 <- ggplot(summarydata[summarydata$CH4_dry < 15,], aes(DATETIME, CH4_dry, color=DWP)) + geom_point()
p2 <- p2 + facet_grid(Source~Rep, scales="free")
p2 <- p2 + theme(axis.text.x = element_text(angle = 90))
print(p2)
save_plot("CH4_source_rep")

for(dwp_core in unique(summarydata$DWP)) {
  printlog("Plotting", dwp_core)
  d <- subset(summarydata, DWP==dwp_core) %>% melt(measure.vars=c("CO2_dry", "CH4_dry"))
  p <- ggplot(d, aes(ELAPSED_TIME/60/60, value, color=Rep)) + geom_point()
  p <- p + facet_grid(variable~Trt, scales="free")
  p <- p + ggtitle(paste("DWP core", dwp_core)) + xlab("Elapsed time (hours)")
  print(p)
  save_plot(paste0("DWP_", dwp_core))
}


printlog("All done with", SCRIPTNAME)
print(sessionInfo())
sink() # close log
