# Process Picarro data for DWP lab experiment
# Ben Bond-Lamberty April 2015

source("0-functions.R")

SCRIPTNAME  	<- "3-fluxes.R"
SUMMARYDATA      <- paste0(OUTPUT_DIR, "summarydata.csv")  # output from script 2


# ==============================================================================
# Main 

sink(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), split=T) # open log

printlog("Welcome to", SCRIPTNAME)

printlog("Reading in summary data...")
summarydata <- read_csv(SUMMARYDATA)

library(lubridate)
summarydata$DATETIME <- ymd_hms(summarydata$DATETIME)
print_dims(summarydata)
print(str(summarydata))

printlog("Reading and merging headspace data...")
hs <- read_csv("data/DWP2013 headspace_cm.csv")
summarydata <- merge(summarydata, hs)
print_dims(summarydata)

printlog("Filtering data...")
fluxdata <- summarydata %>%
  filter(Trt == "injection data", Source == "Core") %>%
  select(DWP_core, Rep, samplenum, m_CO2, m_CH4, ELAPSED_TIME, 
         Site, MinDepth_cm, CoreMassPostInjection_g, headspace_in_core_cm)

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

# The instrument tubing is 455 cm long by ID 1/16"
V_tubing <- (1/16 * 2.54 / 2 ) ^ 2 * pi * 455
# Headspace on the core is 7.3 cm diameter by 4 cm height.
V_headspace <- (7.3 / 2) ^ 2 * pi * fluxdata$headspace_in_core_cm
# Internal volume of Picarro?
V_picarro <- 0
fluxdata$V <- V_tubing + V_headspace + V_picarro

Pa 			<- 101						# kPa				(Richland is ~120 m asl)
R 			<- 8.3145e+3			# cm3 kPa K−1 mol−1
Tair    <- 273.1 + 20     # unknown

# Calculate mass-corrected respiration, µmol/g soil/s
fluxdata$CO2_flux_umol_g_s <- fluxdata$m_CO2 / 1e6 * # from ppm/s to mole fraction/s
  fluxdata$V / fluxdata$CoreMassPostInjection_g * Pa / (R * Tair)
fluxdata$CH4_flux_umol_g_s <- fluxdata$m_CH4 / 1e6 * # from ppm/s to mole fraction/s
  fluxdata$V / fluxdata$CoreMassPostInjection_g * Pa / (R * Tair)

# Convert from µmol/g soil/s to mgC/s
fluxdata$CO2_flux_mgC_s <- fluxdata$CO2_flux_umol_g_s * fluxdata$CoreMassPostInjection_g / # get rid of /g soil
  1e6 * # to mol 
  12 *  # to g C
  1000  # to mg C
fluxdata$CH4_flux_mgC_s <- fluxdata$CH4_flux_umol_g_s * fluxdata$CoreMassPostInjection_g / # get rid of /g soil
  1e6 * # to mol 
  16 *  # to g C
  1000  # to mg C

# Right now looks like I'm off by x1000 (too small)? 

fluxdata <- fluxdata[complete.cases(fluxdata),]

# Each core had 6 ml CH4 injected = 0.006 L / 22.413 L/mol = 0.0002677017802 mol x 12 gC/mol = 3.21 mg C


stop('ok')

meandata$DWP_core <- as.numeric(meandata$DWP_core)

md_long <- melt(meandata, id.vars=c("Trt", "Source", "DWP_core", "Depth_cm", "CoreMassPostInjection_g"))


md_long$value <- md_long$value / md_long$CoreMassPostInjection_g
md_long <- subset(md_long, !is.na(value))

md_long$lbl <- NA
highCO2 <- md_long$variable == "CO2_dry" & md_long$value > .3
highCH4 <- md_long$variable == "CH4_dry" & md_long$value > 0.002
md_long$lbl[highCO2] <- md_long$DWP_core[highCO2]
md_long$lbl[highCH4] <- md_long$DWP_core[highCH4]

p1 <- ggplot(md_long, aes(DWP_core, value, fill=Depth_cm)) + geom_bar(stat="identity")
p1 <- p1 + geom_text(aes(label=lbl), size=3)
p1 <- p1 + facet_grid(variable~., scales="free")
p1 <- p1 + theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(p1)
save_plot("highfluxes")

save_data(md_long)

printlog("All done with", SCRIPTNAME)
print(sessionInfo())
sink() # close log
