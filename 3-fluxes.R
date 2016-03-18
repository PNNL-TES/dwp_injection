# Process Picarro data for DWP lab experiment
# Ben Bond-Lamberty April 2015

source("0-functions.R")

SCRIPTNAME  	<- "3-fluxes.R"
SUMMARYDATA      <- file.path(OUTPUT_DIR, "summarydata.csv")  # output from script 2

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE) # open log

printlog("Welcome to", SCRIPTNAME)

printlog("Reading in summary data...")
summarydata <- read_csv(SUMMARYDATA)

library(lubridate)
summarydata$DATETIME <- ymd_hms(summarydata$DATETIME)  # UTC
summarydata$STARTDATETIME <- ymd_hms(summarydata$STARTDATETIME)  # UTC
print_dims(summarydata)

PREINJECTION_CUTOFF_HRS <- 24   # only use data within 24 hours pre-injection
POSTINJECTION_CUTOFF_HRS <- 72   # only use data within 72 hours post-injection
printlog("PREINJECTION_CUTOFF_HRS =", PREINJECTION_CUTOFF_HRS)
printlog("POSTINJECTION_CUTOFF_HRS =", POSTINJECTION_CUTOFF_HRS)
printlog("Applying cutoff and labeling data...")
summarydata %>%
  filter(ELAPSED_TIME_s > -PREINJECTION_CUTOFF_HRS * 60 * 60) %>%
  filter(ELAPSED_TIME_s < POSTINJECTION_CUTOFF_HRS * 60 * 60) %>%
  mutate(Trt = "Pre-injection") ->
  summarydata
summarydata$Trt[summarydata$ELAPSED_TIME_s > 0] <- "Post-injection"
print_dims(summarydata)

# Observations with negative elapsed times occurred before the injection
printlog("QC plot showing elapsed time and treatment...")
summarydata$RepCore <- paste(summarydata$Rep, summarydata$DWP_core)
p <- ggplot(summarydata, aes(ELAPSED_TIME_s/60/60, max_CO2, color = Trt))
p <- p + geom_point() + geom_vline(xintercept = 0) + xlab("Time since injection (hr)")
p <- p + facet_wrap(~RepCore) + coord_cartesian(xlim = c(-10, 10))
print(p)
save_plot("QC_ELAPSED_TIME_s")
# TODO: note there appears to be a lot of mis-categorized data here

printlog("Reading and merging headspace data...")
hs <- read_csv("data/DWP2013 headspace_cm.csv")
summarydata <- merge(summarydata, hs, all.x = TRUE)
print_dims(summarydata)

printlog("Selecting data...")
fluxdata <- summarydata %>%
  select(Trt, DWP_core, SamplePoint, Injection, Rep, samplenum, ELAPSED_TIME_s, 
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

printlog("Computing pre-injection rates for Injection 1...")
fluxdata %>%
  filter(ELAPSED_TIME_s <= 0 & Injection == 1) %>%
  group_by(Rep, DWP_core) %>%
  summarise(CO2_flux_mgC_s_pre = mean(CO2_flux_mgC_s, na.rm = TRUE),
            CO2_flux_mgC_s_presd = sd(CO2_flux_mgC_s, na.rm = TRUE),
            CH4_flux_mgC_s_pre = mean(CH4_flux_mgC_s, na.rm = TRUE),
            CH4_flux_mgC_s_presd = sd(CH4_flux_mgC_s, na.rm = TRUE)) ->
  fd_preinjection

# Do a bunch of QC plots check pre- versus postinjection fluxes
p <- qplot(as.numeric(DWP_core), CO2_flux_mgC_s_pre, data = fd_preinjection)
p <- p + geom_text(aes(label = DWP_core), size = 4, vjust = -0.5, hjust = -0.5)
p <- p + ggtitle("Pre-injection CO2 flux rates by core")
p <- p + geom_errorbar(aes(ymin = CO2_flux_mgC_s_pre - CO2_flux_mgC_s_presd,
                           ymax = CO2_flux_mgC_s_pre + CO2_flux_mgC_s_presd))
print(p)
save_plot("QC_CO2_preinjection")
p <- qplot(as.numeric(DWP_core), CH4_flux_mgC_s_pre, data = fd_preinjection)
p <- p + geom_text(aes(label = DWP_core), size = 4, vjust = -0.5, hjust = -0.5)
p <- p + ggtitle("Pre-injection CH4 flux rates by core")
p <- p + geom_errorbar(aes(ymin = CH4_flux_mgC_s_pre - CH4_flux_mgC_s_presd,
                           ymax = CH4_flux_mgC_s_pre + CH4_flux_mgC_s_presd))
print(p)
save_plot("QC_CH4_preinjection")

# Plot everything together, by rep
p <- qplot(ELAPSED_TIME_s/60/60, CH4_flux_umol_g_s, data = fluxdata, color = Trt, geom = "line", group = DWP_core) 
p <- p + facet_wrap(~Rep, scales="free") + geom_vline(xintercept=c(-25, 0), linetype=2)
print(p)
save_plot("QC_CH4_rep")
p <- qplot(ELAPSED_TIME_s/60/60, CO2_flux_umol_g_s, data = fluxdata, color = Trt, geom = "line", group = DWP_core) 
p <- p + facet_wrap(~Rep, scales="free") + geom_vline(xintercept=c(-25, 0), linetype=2)
print(p)
save_plot("QC_CO2_rep")


# Individual core graphs, both by themselves and in comparison to others

graph_core <- function(d, fluxdata) {
 
  doplot <- function(gas, d, fluxdata) {
    d <- as.data.frame(d)  # weird
    gasvar <- paste0(gas, "_flux_umol_g_s")
 
    injection <- d$Injection[1]
    rep <- d$Rep[1]
    site <- d$Site[1]
    core <- d$DWP_core[1]
    depth <- d$Depth_cm[1]
    d_samedepth <- subset(fluxdata, Depth_cm == depth)
    d_samedepth$core_rep <- paste(d_samedepth$DWP_core, d_samedepth$Rep)
    printlog("QC preinjection for core", core, "inj", injection, "rep", rep)
    
    p1 <- ggplot(d, aes_string("ELAPSED_TIME_s/60/60", gasvar, color = "Trt")) + geom_point()
    p1 <- p1 + xlab("Time since injection (hr)")
    p1 <- p1 + scale_color_manual(values = c("red", "blue"))
    p1 <- p1 + geom_vline(xintercept = 0, linetype = 2) + ggtitle(paste(gas, "- DWP core", core, rep, site, depth))
    meanpre <- mean(d[d$ELAPSED_TIME_s <= 0, gasvar], na.rm = TRUE)
    meanpost <- mean(d[d$ELAPSED_TIME_s > 0, gasvar], na.rm = TRUE)
    p1 <- p1 + annotate("segment", x = -Inf, xend = 0, y = meanpre, yend = meanpre, color = "blue", linetype = 2)
    p1 <- p1 + annotate("segment", x = 0, xend = Inf, y = meanpost, yend = meanpost, color = "red", linetype = 2)
    p2 <- ggplot(d_samedepth, aes_string("ELAPSED_TIME_s/60/60", gasvar, group = "core_rep")) 
    p2 <- p2 + geom_line(alpha = I(.5)) + geom_line(data = d, color = "red", group = 1)
    p2 <- p2 + xlab("Time since injection (hr)") + ggtitle(paste("...compared to all other", depth, "cores"))
    pdf(file.path(outputdir(), paste0("QC_core_", core, "_", injection, "_", rep, "_", gas, ".pdf")))
    multiplot(p1, p2)
    dev.off()
  }
  
  doplot("CH4", d, fluxdata)
  doplot("CO2", d, fluxdata)
}

fluxdata %>%
  group_by(Injection, DWP_core, Rep) %>%
  do(n = graph_core(., fluxdata))


# ----------------------------------------------------------------------
# December 4, 2015: Some cores are exhibiting wacky pre-injection fluxes.
# Remove them for now (until we decide what to do; see issue #12) TODO
fd_preinjection[fd_preinjection$DWP_core == 16, "CH4_flux_mgC_s_pre"] <- NA
fd_preinjection[fd_preinjection$DWP_core == 35, "CH4_flux_mgC_s_pre"] <- NA
fd_preinjection[fd_preinjection$DWP_core == 9, "CH4_flux_mgC_s_pre"] <- NA
fd_preinjection[fd_preinjection$DWP_core == 4, "CO2_flux_mgC_s_pre"] <- NA
fd_preinjection[fd_preinjection$DWP_core == 7, "CO2_flux_mgC_s_pre"] <- NA
fd_preinjection[fd_preinjection$DWP_core == 9, "CO2_flux_mgC_s_pre"] <- NA

# Merge the preinjection means back into main data
fluxdata <- left_join(fluxdata, fd_preinjection, by=c("Rep", "DWP_core"))
print_dims(fluxdata)

# Just remove 4, 7, 9 entirely for now; see issue #12 (TODO)
fluxdata <- filter(fluxdata, !(DWP_core %in% c(2, 4, 7, 9, 51)))


# ----------------------------------------------------------------------

printlog("Computing cumulative C emission...")
fluxdata %>%
  filter(ELAPSED_TIME_s > 0) %>% # cumulative emissions only after injection
  select(Injection, Rep, DWP_core, ELAPSED_TIME_s, CO2_flux_mgC_s, CH4_flux_mgC_s, 
         CO2_flux_mgC_s_pre, CH4_flux_mgC_s_pre) %>%
  group_by(Injection, Rep, DWP_core) %>%
  arrange(ELAPSED_TIME_s) %>%
  mutate(deltaTime = ELAPSED_TIME_s - lag(ELAPSED_TIME_s),
         # extrapolate pre-injection rates to cumulative fluxes
         CO2_flux_mgC_pre = CO2_flux_mgC_s_pre * deltaTime,
         cumCO2_flux_mgC_pre = c(0, cumsum(CO2_flux_mgC_pre[!is.na(CO2_flux_mgC_pre)])),
         CH4_flux_mgC_pre = CH4_flux_mgC_s_pre * deltaTime,
         cumCH4_flux_mgC_pre = c(0, cumsum(CH4_flux_mgC_pre[!is.na(CH4_flux_mgC_pre)])),
         
         # calculate post-injection cumulative fluxes
         CO2_flux_mgC = CO2_flux_mgC_s * deltaTime,
         cumCO2_flux_mgC = c(0, cumsum(CO2_flux_mgC[!is.na(CO2_flux_mgC)])),
         CH4_flux_mgC = CH4_flux_mgC_s * deltaTime,
         cumCH4_flux_mgC = c(0, cumsum(CH4_flux_mgC[!is.na(CH4_flux_mgC)])),
         
         # calculate net fluxes
         net_cum_CO2_mgC = cumCO2_flux_mgC - cumCO2_flux_mgC_pre,
         net_cum_CH4_mgC = cumCH4_flux_mgC - cumCH4_flux_mgC_pre) %>%
  # remove unneeded fields
  select(-CO2_flux_mgC_s, -CH4_flux_mgC_s, -CO2_flux_mgC, -CH4_flux_mgC,
         -CO2_flux_mgC_s_pre, -CH4_flux_mgC_s_pre, -CO2_flux_mgC_pre, -CH4_flux_mgC_pre) %>%
  # ...and merge back into main data
  left_join(fluxdata, ., by=c("Injection", "Rep", "DWP_core", "ELAPSED_TIME_s")) ->
  fluxdata

print_dims(fluxdata)

# ...and change NAs (pre-injection) to zero
fluxdata$cumCO2_flux_mgC[is.na(fluxdata$cumCO2_flux_mgC)] <- 0.0
fluxdata$cumCH4_flux_mgC[is.na(fluxdata$cumCH4_flux_mgC)] <- 0.0
fluxdata$net_cum_CO2_mgC[is.na(fluxdata$net_cum_CO2_mgC)] <- 0.0
fluxdata$net_cum_CH4_mgC[is.na(fluxdata$net_cum_CH4_mgC)] <- 0.0

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
closelog() # close log
