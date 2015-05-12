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

# Fractional solenoid values mean that the analyzer was shifting
# between two samples. Discard these.
printlog( "Removing fractional MPVPosition" )
rawdata <- subset(rawdata, MPVPosition==trunc(MPVPosition))

# Make true dates
printlog( "Converting date/time info to POSIXct..." )
rawdata$DATETIME <- ymd_hms(paste(rawdata$DATE, rawdata$TIME))

# Assign a different sample number to each sample group (we know we're on a new group when MPVPosition changes)
printlog("Assigning sample numbers...")
rawdata$samplenum <- with(rawdata, cumsum(c(TRUE, !MPVPosition[-length(MPVPosition)] == MPVPosition[-1])))

printlog("Computing CO2 and CH4 models...")
mods <- rawdata %>%
  group_by(samplenum) %>%
  do(co2mod = lm(CO2_dry ~ FRAC_HRS_SINCE_JAN1, data = .),
     ch4mod = lm(CH4_dry ~ FRAC_HRS_SINCE_JAN1, data = .))

modsummary <- mods %>% 
  summarise(R2_CO2=summary(co2mod)$r.squared,
            R2_CH4=summary(ch4mod)$r.squared,
            m_CO2=coef(co2mod)[2],
            m_CH4=coef(co2mod)[2])
modsummary$samplenum <- unique(rawdata$samplenum)

printlog( "Computing summary statistics for each sample..." )
summarydata <- rawdata %>%
  group_by(samplenum) %>%
  summarise(
    EPOCH_TIME = mean(EPOCH_TIME),
    DATETIME = mean(DATETIME),
    ALARM_STATUS = mean(ALARM_STATUS),
    INST_STATUS	= mean(INST_STATUS),
    N = n(),
    Valve	= mean(MPVPosition),
    CH4_dry = mean(CH4_dry),
    CO2_dry = mean(CO2_dry),
    h2o_reported = mean(h2o_reported),
    Injection = mean(injection),
    Rep = paste(unique(rep), collapse=","),
    Trt = paste(unique(trt), collapse=",")
  )

printlog("Merging summary and model data...")
summarydata <- merge(summarydata, modsummary)

# Load MPVPosition map
printlog("Loading valve map data and merging...")
valvemap <- read_csv("data/DWP_valve assignment 3 March2105.csv")
summarydata <- merge(summarydata, valvemap)
summarydata$Source <- "Core"
summarydata$Source[summarydata$DWP == "ambient"] <- "Ambient"
summarydata$Source[summarydata$DWP == "CH4 blank"] <- "Blank"

printlog("Computing elapsed times...")
summarydata <- summarydata %>% 
  group_by(Rep, Valve, Trt) %>% 
  mutate(ELAPSED_TIME=as.numeric(DATETIME - min(DATETIME)))

# Done!
save_data(summarydata, scriptfolder=FALSE)

printlog("All done with", SCRIPTNAME)
print(sessionInfo())
sink() # close log
