# Process Picarro data for DWP lab experiment
# Ben Bond-Lamberty April 2015

source("0-functions.R")

SCRIPTNAME  	<- "2-summarize.R"
RAWDATA      <- paste0(OUTPUT_DIR, "rawdata.csv.gz")  # output from script 1

# ==============================================================================
# Main 

sink(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), split=T) # open log

printlog("Welcome to", SCRIPTNAME)

printlog("Reading in raw data...")
rawdata <- gzfile(RAWDATA) %>% readr::read_csv()
print_dims(rawdata)
print(summary(rawdata))

rawdata %>% group_by(rep,trt) %>% summarise(n()) %>% print()

# Fractional solenoid values mean that the analyzer was shifting
# between two samples. Discard these.
printlog( "Removing fractional MPVPosition" )
rawdata <- subset(rawdata, MPVPosition == trunc(MPVPosition))

# Make true dates
printlog( "Converting date/time info to POSIXct..." )
rawdata$DATETIME <- ymd_hms(paste(rawdata$DATE, rawdata$TIME), tz = "UTC") # tz="America/Los_Angeles")
printlog("First timestamp:")
print(min(rawdata$DATETIME))
printlog("Last timestamp:")
print(max(rawdata$DATETIME))

printlog( "Sorting by date..." )
rawdata <- arrange(rawdata, DATETIME)

# Assign a different sample number to each sample group 
# (we know we're on a new sample when MPVPosition changes)
printlog("Assigning sample numbers...")
oldsampleflag <- with(rawdata, c(FALSE, 
                                  MPVPosition[-length(MPVPosition)] == MPVPosition[-1] &
                                    trt[-length(trt)] == trt[-1]))
rawdata$samplenum <- cumsum(!oldsampleflag)

printlog("Sample number summary:")
rawdata %>% 
  group_by(trt, rep) %>% 
  summarise(unique_samplenums = length(unique(samplenum)),
            unique_MPVPositions = length(unique(MPVPosition))) %>%
  print()

printlog("Computing elapsed seconds...")
rawdata <- rawdata %>%
  group_by(samplenum) %>%
#  filter(trt == "injection data") %>%
  mutate(elapsed_seconds = (FRAC_HRS_SINCE_JAN1 - min(FRAC_HRS_SINCE_JAN1)) * 60 * 60)

printlog( "Computing summary statistics for each sample..." )
# The window in which we look for min and max concentrations
MAX_MINCONC_TIME <- 10  # the minimum concentration has to occur in first 10 s
MAX_MAXCONC_TIME <- 45  # the maximum concentration can't occur in first 2 s

# We want to apply different criteria here, so three different pipelines
# to compute the min and max gas concentrations
summarydata_min <- rawdata %>%
  filter(elapsed_seconds <= MAX_MINCONC_TIME) %>%
  group_by(samplenum) %>%
  summarise(
    min_CO2 = min(CO2_dry),
    min_CO2_time = nth(elapsed_seconds, which.min(CO2_dry)),
    min_CH4 = min(CH4_dry),
    min_CH4_time = nth(elapsed_seconds, which.min(CH4_dry))
  )

# Now we want to look for the max concentration AFTER the minimum
rawdata_temp <- rawdata %>%
  left_join(summarydata_min, by="samplenum") 
summarydata_maxCO2 <- rawdata_temp %>%
  filter(elapsed_seconds > min_CO2_time & elapsed_seconds < MAX_MAXCONC_TIME) %>%
  summarise(
    max_CO2 = max(CO2_dry),
    max_CO2_time = nth(elapsed_seconds, which.max(CO2_dry))
  )
summarydata_maxCH4 <- rawdata_temp %>%
  filter(elapsed_seconds > min_CH4_time & elapsed_seconds < MAX_MAXCONC_TIME) %>%
  summarise(
    max_CH4 = max(CH4_dry),
    max_CH4_time = nth(elapsed_seconds, which.max(CH4_dry))
  )

summarydata_other <- rawdata %>%
  group_by(samplenum) %>%
  summarise(
    Injection = mean(injection),
    Rep = paste(unique(rep), collapse=","),
    Trt = paste(unique(trt), collapse=","),
    
    DATETIME = mean(DATETIME),
    N = n(),
    Valve	= mean(MPVPosition),
    h2o_reported = mean(h2o_reported)
  )

printlog("summarydata_other summary:")
summarydata_other %>% group_by(Rep, Trt) %>% summarise(n()) %>% print()

# Merge pieces together to form final summary data set
summarydata_raw <- summarydata_other %>%
  left_join(summarydata_min, by="samplenum") %>%
  left_join(summarydata_maxCO2, by="samplenum") %>% 
  left_join(summarydata_maxCH4, by="samplenum")

printlog("summarydata_raw summary:")
summarydata_raw %>% group_by(Rep, Trt) %>% summarise(n()) %>% print()

# Load MPVPosition map
printlog("Loading valve map data and merging...")
valvemap <- read_csv("data/DWP_valve assignment 3 March2105.csv")
summarydata <- merge(summarydata_raw, valvemap, by=c("Injection", "Rep", "Valve"))

printlog("Post valvemap merge summarydata summary:")
summarydata %>% group_by(Rep, Trt) %>% summarise(n()) %>% print()

# Injection 2 was different than injection 1: it consisted of three cores in sequence,
# always on valve 12, monitored continuously (alternating with an ambient on valve 10).
# Valve 11 appears in the raw data but is an artifact. From Sarah Fansler's email:
# Core #4 injected with 6 ml methane at 3 pm on June 26.
# Core #2 injected with 6 ml methane at 10:22 am on June 30.
# Core #7 injected at 10:56 am on July 3.
# Data pulled at 3:25 pm on 7 July.
printlog("Splitting off injection 2 data...")
inj2 <- with(summarydata, Injection == 2 & Valve == 12)
summarydata$DWP_core[inj2] <- "7"
# Note that the Picarro is on UTC already
core2time <- with_tz(ymd_hm("2015-07-03 10:56", tz = "America/Los_Angeles"), tz = "UTC")
core4time <- with_tz(ymd_hm("2015-06-30 10:22", tz = "America/Los_Angeles"), tz = "UTC")
summarydata$DWP_core[inj2 & summarydata$DATETIME < core2time] <- "2"
summarydata$DWP_core[inj2 & summarydata$DATETIME < core4time] <- "4"

printlog("Post injection 2 removal summarydata summary:")
summarydata %>% group_by(Rep, Trt) %>% summarise(n()) %>% print()

# Assign 'Source' field
summarydata$Source <- "Core"
summarydata$Source[summarydata$DWP_core == "ambient"] <- "Ambient"
summarydata$Source[summarydata$DWP_core == "CH4 blank"] <- "Blank"

# Load field data to get depth information
printlog("Loading field data and merging...")
fielddata <- read_csv("data/DWP Core field data.csv")
fielddata$Notes <- fielddata$SampleDate <- NULL
summarydata <- merge(summarydata, fielddata, all.x = TRUE)
summarydata$Depth_cm[summarydata$Source == "Ambient"] <- "Ambient"
summarydata$Depth_cm[summarydata$Source == "Blank"] <- "Blank"

printlog("Post field data merge summarydata summary:")
summarydata %>% group_by(Rep, Trt) %>% summarise(n()) %>% print()

# Load mass data
printlog("Loading mass data and merging...")
massdata <- read_csv("Core inj log_mass 9March2015.csv", datadir = "data/")
summarydata <- merge(summarydata, massdata, by = c("Injection", "Rep", "DWP_core"), all.x = TRUE)

printlog("Post mass data merge summarydata summary:")
summarydata %>% group_by(Rep, Trt) %>% summarise(n()) %>% print()

printlog("Computing min depth...")
summarydata$MinDepth_cm <- as.numeric(str_extract(summarydata$Depth_cm, "^[0-9]*"))


# Done!

save_data(summarydata, scriptfolder=FALSE)

printlog("All done with", SCRIPTNAME)
print(sessionInfo())
sink() # close log
