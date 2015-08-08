# Process Picarro data for DWP lab experiment
# Ben Bond-Lamberty April 2015

source("0-functions.R")

SCRIPTNAME  	<- "3-fluxes.R"
SUMMARYDATA      <- file.path(OUTPUT_DIR, "summarydata.csv")  # output from script 2

# ==============================================================================
# Main 

sink(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), split=T) # open log

printlog("Welcome to", SCRIPTNAME)

printlog("Reading in summary data...")
summarydata <- read_csv(SUMMARYDATA)

library(lubridate)
summarydata$DATETIME <- ymd_hms(summarydata$DATETIME)  # UTC
summarydata$STARTDATETIME <- mdy_hm(summarydata$Start, tz="America/Los_Angeles")
print_dims(summarydata)
print(str(summarydata))

# Compute ELAPSED_TIME based on 'Start' field (originally in the valve map data)
printlog("Computing elapsed times...")
summarydata$ELAPSED_TIME <- with(summarydata, as.numeric(difftime(DATETIME, STARTDATETIME, units="secs")))

# Observations with negative elapsed times occurred before the injection
printlog("QC plot showing elapsed time and treatment...")
p <- ggplot(summarydata, aes(ELAPSED_TIME/60/60, max_CO2, color=Trt))
p <- p + geom_point() + geom_vline()
p <- p + facet_wrap(~DWP_core) + xlim(c(-10, 10))
print(p)
save_plot("QC_elapsed_time")
# TODO: note there appears to be a lot of mis-categorized data here

printlog("Reading and merging headspace data...")
hs <- read_csv("data/DWP2013 headspace_cm.csv")
summarydata <- merge(summarydata, hs, all.x=TRUE)
print_dims(summarydata)

printlog("Filtering data...")
fluxdata <- summarydata %>%
  select(Trt, DWP_core, SamplePoint, Injection, Rep, samplenum, ELAPSED_TIME, 
         min_CO2, min_CO2_time, min_CH4, min_CH4_time, 
         max_CO2, max_CO2_time, max_CH4, max_CH4_time,
         Site, Depth_cm, MinDepth_cm, CoreMassPostInjection_g, headspace_in_core_cm)

# ----------------------------------------------------------------------

printlog("Computing flux rates...")
# Flux rates are ppm/s (CO2) or ppb/s (CH4) in `summarydata`
# We want to convert these to mg C/g soil/s
# A = dC/dt * V/M * Pa/RT (cf. Steduto et al. 2002), where
# 	A is the flux (µmol/g/s)
#	  dC/dt is raw respiration as above (mole fraction/s)
# 	V is total chamber volume (cm3)
#	  M is [dry] soil mass (g)
#	  Pa is atmospheric pressure (kPa)
#	  R is universal gas constant (8.3 x 10^3 cm3 kPa mol-1 K-1)
#	  T is air temperature (K)

# Each core had 6 ml CH4 injected = 0.006 L / 22.413 L/mol = 0.0002677017802 mol x 12 gC/mol = 3.21 mg C

# The instrument tubing is 455 cm long by ID 1/16"
V_tubing <- (1/16 * 2.54 / 2 ) ^ 2 * pi * 455
# Headspace on the core is 7.3 cm diameter by 4 cm height.
V_headspace <- (7.3 / 2) ^ 2 * pi * fluxdata$headspace_in_core_cm
# Internal volume of Picarro?
V_picarro <- 0
V_cm3 <- V_tubing + V_headspace + V_picarro

Pa 			<- 101						# kPa				(Richland is ~120 m asl)
R 			<- 8.3145e+3			# cm3 kPa K−1 mol−1
Tair    <- 273.1 + 20     # unknown

m_CO2 <- with(fluxdata, (max_CO2 - min_CO2) / (max_CO2_time - min_CO2_time))  # ppm/s
m_CH4 <- with(fluxdata, (max_CH4 - min_CH4) / (max_CH4_time - min_CH4_time))  # ppb/s
fluxdata$V_cm3 <- V_cm3

# Calculate mass-corrected respiration, µmol/g soil/s
fluxdata$CO2_flux_umol_g_s <- m_CO2 / 1 * # from ppm/s to µmol/s
  V_cm3 / fluxdata$CoreMassPostInjection_g * Pa / (R * Tair)
fluxdata$CH4_flux_umol_g_s <- m_CH4 / 1e3 * # from ppb/s to µmol/s
  V_cm3 / fluxdata$CoreMassPostInjection_g * Pa / (R * Tair)

# Convert from µmol/g soil/s to mgC/s
fluxdata$CO2_flux_mgC_s <- fluxdata$CO2_flux_umol_g_s * fluxdata$CoreMassPostInjection_g / # get rid of /g soil
  1e6 * # to mol 
  12 *  # to g C
  1000  # to mg C
fluxdata$CH4_flux_mgC_s <- fluxdata$CH4_flux_umol_g_s * fluxdata$CoreMassPostInjection_g / # get rid of /g soil
  1e6 * # to mol 
  16 *  # to g C
  1000  # to mg C

print_dims(fluxdata)
printlog("Removing incomplete cases...")
fluxdata <- fluxdata[complete.cases(fluxdata),]
print_dims(fluxdata)

# ----------------------------------------------------------------------

printlog("Computing pre-injection rates...")
fd_preinjection <- fluxdata %>%
  filter(ELAPSED_TIME <= 0) %>%
  group_by(Rep, DWP_core) %>%
  summarise(CO2_flux_umol_g_s_pre = mean(CO2_flux_umol_g_s, na.rm=TRUE),
            CO2_flux_umol_g_s_presd = sd(CO2_flux_umol_g_s, na.rm=TRUE),
            CH4_flux_umol_g_s_pre = mean(CH4_flux_umol_g_s, na.rm=TRUE),
            CH4_flux_umol_g_s_presd = sd(CH4_flux_umol_g_s, na.rm=TRUE))

# Do a bunch of QC plots check pre- versus postinjection fluxes
p <- qplot(as.numeric(DWP_core), CO2_flux_umol_g_s_pre, data=fd_preinjection)
p <- p + geom_text(aes(label=DWP_core), size=4, vjust=-.5, hjust=-.5)
p <- p + geom_errorbar(aes(ymin=CO2_flux_umol_g_s_pre - CO2_flux_umol_g_s_presd,
                           ymax=CO2_flux_umol_g_s_pre + CO2_flux_umol_g_s_presd))
print(p)
save_plot("QC_CO2_preinjection")
p <- qplot(as.numeric(DWP_core), CH4_flux_umol_g_s_pre, data=fd_preinjection)
p <- p + geom_text(aes(label=DWP_core), size=4, vjust=-.5, hjust=-.5)
p <- p + geom_errorbar(aes(ymin=CH4_flux_umol_g_s_pre - CH4_flux_umol_g_s_presd,
                           ymax=CH4_flux_umol_g_s_pre + CH4_flux_umol_g_s_presd))
print(p)
save_plot("QC_CH4_preinjection")

for(dwp in unique(fluxdata$DWP_core)) {
  printlog("QC preinjection for core", dwp)
  d <- subset(fluxdata, DWP_core==dwp)
  p1 <- ggplot(d, aes(ELAPSED_TIME/60/60, CO2_flux_umol_g_s, color=Trt)) + geom_point()
  p1 <- p1 + scale_color_manual(values=c("red", "blue"))
  p1 <- p1 + geom_vline(linetype=2) + ggtitle(paste("DWP core", dwp))
  p1 <- p1 + geom_hline(yintercept=mean(d[d$ELAPSED_TIME <= 0, "CO2_flux_umol_g_s"], na.rm=TRUE), color="red", linetype=2) 
  p1 <- p1 + geom_hline(yintercept=mean(d[d$ELAPSED_TIME > 0, "CO2_flux_umol_g_s"], na.rm=TRUE), color="blue", linetype=2) 
  p2 <- ggplot(fluxdata, aes(ELAPSED_TIME/60/60, CO2_flux_umol_g_s, group=DWP_core)) 
  p2 <- p2 + geom_line(alpha=I(.5)) + geom_line(data=d, color="red")
  pdf(file.path(outputdir(), paste0("QC_core_", dwp, "_CO2.pdf")))
  multiplot(p1, p2)
  dev.off()
  # CH4
  p1 <- ggplot(d, aes(ELAPSED_TIME/60/60, CH4_flux_umol_g_s, color=Trt)) + geom_point()
  p1 <- p1 + scale_color_manual(values=c("red", "blue")) + scale_y_log10()
  p1 <- p1 + geom_vline(linetype=2) + ggtitle(paste("DWP core", dwp))
  p1 <- p1 + geom_hline(yintercept=mean(d[d$ELAPSED_TIME <= 0, "CH4_flux_umol_g_s"], na.rm=TRUE), color="red", linetype=2) 
  p1 <- p1 + geom_hline(yintercept=mean(d[d$ELAPSED_TIME > 0, "CH4_flux_umol_g_s"], na.rm=TRUE), color="blue", linetype=2) 
  p2 <- ggplot(fluxdata, aes(ELAPSED_TIME/60/60, CH4_flux_umol_g_s, group=DWP_core)) 
  p2 <- p2 + geom_line(alpha=I(.5)) + geom_line(data=d, color="red") + scale_y_log10()
  pdf(file.path(outputdir(), paste0("QC_core_", dwp, "_CH4.pdf")))
  multiplot(p1, p2)
  dev.off()
  
}

# merge back into main data
fluxdata <- left_join(fluxdata, fd_preinjection, by=c("Rep", "DWP_core"))
print_dims(fluxdata)

# ----------------------------------------------------------------------

printlog("Computing cumulative C emission...")
fluxdata <- fluxdata %>%
  filter(ELAPSED_TIME > 0) %>% # cumulative emissions only after injection
  select(Rep, DWP_core, ELAPSED_TIME, CO2_flux_mgC_s, CH4_flux_mgC_s) %>%
  group_by(Rep, DWP_core) %>%
  arrange(ELAPSED_TIME) %>%
  mutate(CO2_flux_mgC = CO2_flux_mgC_s * (ELAPSED_TIME - lag(ELAPSED_TIME)),
         cumCO2_flux_mgC = c(0, cumsum(CO2_flux_mgC[!is.na(CO2_flux_mgC)])),
         CH4_flux_mgC = CH4_flux_mgC_s * (ELAPSED_TIME - lag(ELAPSED_TIME)),
         cumCH4_flux_mgC = c(0, cumsum(CH4_flux_mgC[!is.na(CH4_flux_mgC)]))) %>%
  select(-CO2_flux_mgC_s, -CH4_flux_mgC_s, -CO2_flux_mgC, -CH4_flux_mgC) %>%
  # merge back into main data
  left_join(fluxdata, ., by=c("Rep", "DWP_core", "ELAPSED_TIME"))

print_dims(fluxdata)

# ...and change NAs (pre-injection) to zero
fluxdata$cumCO2_flux_mgC[is.na(fluxdata$cumCO2_flux_mgC)] <- 0.0
fluxdata$cumCH4_flux_mgC[is.na(fluxdata$cumCH4_flux_mgC)] <- 0.0

# ----------------------------------------------------------------------

# We ran a subsequent check using cores, 2, 4, and 7, monitoring them 
# continuously to make sure we didn't miss any methane or CO2 'burps'. 
# Split off those data separately.
printlog("Splitting data by injection...")
fluxdata_247check <- filter(fluxdata, Injection == 2)
fluxdata <- filter(fluxdata, Injection != 2)

save_data(fluxdata, scriptfolder=FALSE)
save_data(fluxdata_247check, scriptfolder=FALSE)

printlog("All done with", SCRIPTNAME)
print(sessionInfo())
sink() # close log
