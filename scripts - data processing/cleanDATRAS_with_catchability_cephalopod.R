#############################################################
#### Code to download and clean survey data from DATRAS
#### Coding: originally started by Aurore Maureaud, March 2021 

#### updated: Daniel van Denderen July 2021
#### Coding updated for wingspread / catch per unit area information
#### Coding updated with gear efficiency estimates based on Walker et al. 2017
#### see https://doi.org/10.1093/icesjms/fsw250

#### updated Daniel van Denderen and Daniel Ottmann November 2022
#### Coding updated to include the class cephalopoda
#############################################################
Sys.setenv(LANG = "en")
rm(list=ls())


##########################################################################################
#### LOAD LIBRARIES
##########################################################################################
library(data.table)
library(dplyr)
library(worms)
library(worrms)
library(crul)
library(urltools) # check if you have Rcpp installed, no need to load

# ------------------------------------------------------------------------------
# LOAD FILES DATRAS -- done once -- 2021 July 17 -- continue next step
# ------------------------------------------------------------------------------
library(icesDatras)
last.year <- 2019

# Haul info from Datras
hh.ns <- getDATRAS(record='HH', survey='NS-IBTS', years=c(1967:last.year), quarters=c(1,3))
hh.baltic <- getDATRAS(record='HH', survey='BITS', years=c(1991:last.year), quarters=c(1,4))
hh.evhoe <- getDATRAS(record='HH', survey='EVHOE', years=c(1997:last.year), quarters=4)
hh.cgfs <- getDATRAS(record='HH', survey='FR-CGFS', years=c(1998:last.year), quarters=4)
hh.igfs <- getDATRAS(record='HH', survey='IE-IGFS', years=c(2003:last.year), quarters=4)
hh.nigfs <- getDATRAS(record='HH', survey='NIGFS', years=c(2005:last.year), quarters=c(1:4))
hh.pt <- getDATRAS(record='HH', survey='PT-IBTS', years=c(2002:last.year), quarters=c(3:4))
hh.rock <- getDATRAS(record='HH', survey='ROCKALL', years=c(1999:2009), quarters=3)
hh.scorock <- getDATRAS(record='HH', survey='SCOROC', years=c(2011:last.year), quarters=3)
hh.swc <- getDATRAS(record='HH', survey='SWC-IBTS', years=c(1985:2010), quarters=c(1:4))
hh.scowcgfs <- getDATRAS(record='HH', survey='SCOWCGFS', years=c(2011:last.year), quarters=c(1:4))
hh.canmar <- getDATRAS(record='HH', survey='Can-Mar', years=c(1970:last.year), quarters=c(1:4))

hh <- rbind(hh.ns, hh.baltic, hh.evhoe, hh.cgfs, hh.igfs, hh.nigfs, hh.pt, hh.rock, hh.scorock, 
            hh.swc, hh.scowcgfs, hh.canmar)

# save (save outside github as file is large)
# setwd("")
# save(hh, file='haulinfoDatras_210717.RData')

# Length info from DATRAS
hl.ns <- getDATRAS(record='HL', survey='NS-IBTS', years=c(1967:last.year), quarters=c(1,3))
hl.baltic <- getDATRAS(record='HL', survey='BITS', years=c(1991:last.year), quarters=c(1,4))
hl.evhoe <- getDATRAS(record='HL', survey='EVHOE', years=c(1997:last.year), quarters=4)
hl.cgfs <- getDATRAS(record='HL', survey='FR-CGFS', years=c(1998:last.year), quarters=4)
hl.igfs <- getDATRAS(record='HL', survey='IE-IGFS', years=c(2003:last.year), quarters=4)
hl.nigfs <- getDATRAS(record='HL', survey='NIGFS', years=c(2005:last.year), quarters=c(1:4))
hl.pt <- getDATRAS(record='HL', survey='PT-IBTS', years=c(2002:last.year), quarters=c(3:4))
hl.rock <- getDATRAS(record='HL', survey='ROCKALL', years=c(1999:2009), quarters=3)
hl.scorock <- getDATRAS(record='HL', survey='SCOROC', years=c(2011:last.year), quarters=3)
hl.swc <- getDATRAS(record='HL', survey='SWC-IBTS', years=c(1985:2010), quarters=c(1:4))
hl.scowcgfs <- getDATRAS(record='HL', survey='SCOWCGFS', years=c(2011:last.year), quarters=c(1:4))
hl.canmar <- getDATRAS(record='HL', survey='Can-Mar', years=c(1970:last.year), quarters=c(1:4))

hl <- rbind(hl.ns, hl.baltic, hl.evhoe, hl.cgfs, hl.igfs, hl.nigfs, hl.pt, hl.rock, hl.scorock,
            hl.swc, hl.scowcgfs, hl.canmar)

# save (save outside github as file is large)
# setwd("")
#save(hl, file='lengthinfoDatras_210717.RData')

rm(hl.ns, hl.baltic, hl.evhoe, hl.cgfs, hl.igfs, hl.nigfs, hl.pt, hl.rock, hl.scorock, hl.swc, hl.scowcgfs, hl.canmar, 
   hh.ns, hh.baltic, hh.evhoe, hh.cgfs, hh.igfs, hh.nigfs, hh.pt, hh.rock, hh.scorock, hh.swc, hh.scowcgfs, hh.canmar)


# ------------------------------------------------------------------------------
# CREATE A UNIQUE HAUL ID
# ------------------------------------------------------------------------------
# load hl and hh (this is done outside github as files are large - Dani to decide)
load('data/data datras/haulinfoDatras_210717.RData')
load('data/data datras/lengthinfoDatras_210717.RData')

hl$HaulID <- paste(hl$Survey, hl$Year,hl$Quarter, hl$Country, hl$Ship, hl$Gear, hl$StNo, hl$HaulNo)
hh$HaulID <- paste(hh$Survey, hh$Year,hh$Quarter, hh$Country, hh$Ship, hh$Gear, hh$StNo, hh$HaulNo)

# Is the HaulID unique?
hhn <- unique(hh$HaulID)
length(hhn)==nrow(hh)

# check which one is not
# pb <- c()
# for (i in 1:length(hhn)){
#  j <- which(hh$HaulID==hhn[i])
#  if(length(j)>1){pb <- hhn[i]}
# }

# > pb
# [1] "NS-IBTS 1995 1 NA AA36 GOV 999 999"
hh <- hh %>% filter(HaulID!="NS-IBTS 1995 1 NA AA36 GOV 999 999") # remove the non-unique HaulID in hh and hl
hl <- hl %>% filter(HaulID!="NS-IBTS 1995 1 NA AA36 GOV 999 999")

# Only keep hauls where there is the length composition. 70273 hauls in hh and 70273 in hl
hh <- subset(hh, hh$HaulID %in% hl$HaulID)
hl <- subset(hl, hl$HaulID %in% hh$HaulID)

# ------------------------------------------------------------------------------
# MERGE HH and HL FILES
# ------------------------------------------------------------------------------

haulidhl <- sort(unique(hl$HaulID))
haulidhh <- sort(unique(hh$HaulID))
identical(haulidhh, haulidhl)

# remove some columns in hl
hl$SweepLngt <- hl$SpecCodeType <- hl$SpecCode <- hl$Sex <- hl$DateofCalculation <- hl$RecordType <- hl$GearEx <- NULL

# remove some columns in hh
hh$DateofCalculation <- hh$ThClineDepth <- hh$ThermoCline <- hh$SwellHeight <- hh$SwellDir <- hh$WindSpeed <- hh$WindDir <- hh$BotCurSpeed <- NULL
hh$BotCurDir <- hh$SurCurSpeed <- hh$SurCurDir <- hh$SpeedWater <- hh$TowDir <- hh$WgtGroundRope <- hh$KiteDim <- hh$Buoyancy <- NULL
hh$DoorWgt <- hh$DoorSurface <- hh$WarpDen <- hh$Warpdia <- hh$Warplngt <- hh$Tickler <- hh$Rigging  <- NULL
hh$HydroStNo <- hh$HaulLat <-  hh$HaulLong <- hh$DayNight <- hh$Stratum <- hh$TimeShot <- hh$Day <- hh$RecordType <- hh$GearExp <- hh$DoorType <- NULL


#survey <- merge(hh, hl, by='HaulID', all.x=FALSE, all.y=TRUE)
survey <- right_join(hh, hl, by=c('HaulID','Survey','Quarter','Country','Ship','Gear','StNo','HaulNo','Year'))
nrow(survey)==nrow(hl)

survey <- survey %>% 
  dplyr::rename(SBT = BotTemp,
                SST = SurTemp,
                Speed = GroundSpeed,
                AphiaID = Valid_Aphia)

### Check if the HaulID is unique
### Not the case for the baltic sea, a lot of duplicates!!!
#ids <- unique(hh$HaulID)
#pb <- vector()
# for(i in 1:length(ids)){
#   x <- which(hh$HaulID==ids[i])
#   if(length(x)>1){pb[length(pb)+1] <- ids[i]}
# }
# print(pb) # dim 0 ok!


# ------------------------------------------------------------------------------
# REMOVE INVALID DATA
# ------------------------------------------------------------------------------
survey <- survey %>% 
  filter(HaulVal %in% 'V', #Remove invalid hauls
         !is.na(AphiaID), # Remove invalid species records
         SpecVal %in% c(1,10,4,7),
         DataType %in% c('S','R','C'))


# ------------------------------------------------------------------------------
# RESCALE DATA INTO ABUNDANCE FOR THE HAUL DURATION AND ABUNDANCE AT LENGTH
# ------------------------------------------------------------------------------
# If Data Type=='C', abundance at length already re-adjusted with time so get back the abundance for the actual duration of the haul.
# If data type=='R', abundance at length is multiplied by sub-factor and adjusted to time
survey$CatCatchWgt = as.numeric(survey$CatCatchWgt)

survey <- survey %>% 
  mutate(HLNoAtLngt = case_when(DataType=='C' ~ HLNoAtLngt*SubFactor*HaulDur/60,
                                DataType %in% c('S','R') ~ HLNoAtLngt*SubFactor),
         TotalNo = case_when(DataType=='C' ~ TotalNo*HaulDur/60, 
                             DataType %in% c('S','R') ~ TotalNo),
         CatCatchWgt = case_when(DataType=='C' ~ CatCatchWgt*HaulDur/60,
                                 DataType %in% c('S','R') ~ CatCatchWgt)) %>% 
  select(-HaulVal, -DataType, -StdSpecRecCode, -SpecVal, -SubWgt, -SubFactor) %>% 
  mutate(Survey = if_else(Survey=='SCOWCGFS', 'SWC-IBTS', Survey)) %>% 
  mutate(Survey = if_else(Survey=='SCOROC','ROCKALL',Survey)) %>% 
  filter(!(Survey=="NS-IBTS" & BySpecRecCode %in% c(0,2,3,4,5)), # remove hauls where not all species are recorded
         !(Survey=="BITS" & BySpecRecCode==0))
 
# ------------------------------------------------------------------------------
# GET THE SWEPT AREA in km2
# ------------------------------------------------------------------------------
survey <- survey %>% 
  mutate(WingSpread = replace(WingSpread, WingSpread==-9, NA),
         DoorSpread = replace(DoorSpread, DoorSpread==-9, NA),
         Speed = replace(Speed, Speed==-9, NA),
         Distance = replace(Distance, Distance==-9, NA),
         Depth = replace(Depth, Depth==-9, NA))

survey <- survey %>% 
  mutate(WingSpread = replace(WingSpread, WingSpread== 0, NA),
         DoorSpread = replace(DoorSpread, DoorSpread== 0, NA),
         Distance = replace(Distance, Distance == 0, NA))

# select only certain gears 
# summary of gears per survey
gears <- data.frame(survey) %>% 
  dplyr::group_by(Survey, Gear) %>% 
  dplyr::summarise(hauls = length(unique(HaulID)), years = length(unique(Year))) %>% 
  select(Survey, Gear, hauls, years)

# only select certain gears per survey (GOV and/or most dominant in cases without GOV)
survey <- survey %>% 
  filter(!(Survey=="NS-IBTS" & Gear %in% c('ABD', 'BOT', 'DHT', 'FOT', 'GRT', 'H18', 'HOB', 'HT', 'KAB', 'VIN')),
         !(Survey=="BITS" & Gear %in% c('CAM', 'CHP', 'DT', 'EGY', 'ESB', 'EXP', 'FOT', 'GRT', 'H20', 'HAK', 'LBT','SON')),
         !(Survey=="PT-IBTS" & Gear=='CAR'),
         !(Survey=="Can-Mar" & Gear=='Y36'))

source('scripts - data processing/source_DATRAS_wing_doorspread.r')

# ------------------------------------------------------------------------------
# GET CPUEs AND RIGHT COLUMNS NAMES
# ------------------------------------------------------------------------------

# Remove data without length composition or negative values
xx <- subset(survey, HLNoAtLngt<0 | is.na(LngtClass))
no_length_hauls <- sort(unique(xx$HaulID))

# Only keep abundances/weight
survey <- survey %>%
  filter(!(HaulID %in% no_length_hauls)) %>% # remove hauls without length data
  mutate(numcpue = TotalNo/Area.swept, # abundance/km2
         wtcpue = CatCatchWgt/(Area.swept*1000), #weight in kg/km2
         numh = (TotalNo*60)/HaulDur, # abundance/hour
         wgth = CatCatchWgt*60/(HaulDur*1000), #weight in kg/h
         num = TotalNo, #raw number of individuals
         wgt = CatCatchWgt/1000, # raw weight in kg         
         numlencpue = HLNoAtLngt/Area.swept, #abundance/km2 per length class
         numlenh = HLNoAtLngt*60/HaulDur, #abundance/h per length class
         Season = 'NA',
         Depth = replace(Depth, Depth<0, NA),
         SBT = replace(SBT, SBT<0, NA),
         SST = replace(SST, SST<0, NA),
         LngtClass = ifelse(LngtCode %in% c('.','0'), LngtClass*0.1, LngtClass)) %>% # fix unit of length class
  dplyr::rename(Length = LngtClass) %>% 
  # group_by(Survey, HaulID, StatRec, Year, Month, Quarter, Season, ShootLat, ShootLong, HaulDur, Area.swept, Gear, Depth, SBT, SST, AphiaID, Length) %>%
  # summarize_at(.vars=c('numcpue', 'wtcpue', 'numh', 'wgth', 'num', 'wgt', 'numlencpue','numlenh'), .funs=function(x) sum(x, na.rm=T)) %>%
  select(Survey, HaulID, StatRec, Year, Month, Quarter, Season, ShootLat, ShootLong, HaulDur, Area.swept, Area.doors, Gear, Depth, SBT, SST,
         AphiaID, CatIdentifier, numcpue, wtcpue, numh, wgth, num, wgt, Length, numlencpue, numlenh)
survey <- data.frame(survey)

# ------------------------------------------------------------------------------
# Clean species names
# ------------------------------------------------------------------------------

survey$Species <- NA
survey$AphiaID <- as.numeric(survey$AphiaID)
dat.ices <- survey
aphia_list <- unique(dat.ices$AphiaID)
aphia_list <- aphia_list[!duplicated(aphia_list)]

# creating taxonomy tables for each species
spt <- c(seq(1,length(aphia_list),by=50),length(aphia_list))
my_sp_taxo <- c()
for (j in 1:(length(spt)-1)){
spsub <- wm_record(aphia_list[spt[j]:spt[j+1]])
my_sp_taxo <- rbind(my_sp_taxo,spsub)
}

# row binds all the results and pass to data frame. 
df_test <- as.data.frame(my_sp_taxo)
df_test$url <- df_test$lsid <- df_test$citation <- NULL
df_test$isExtinct <- df_test$modified <- df_test$valid_authority <- df_test$unacceptreason <- NULL
df_test$authority <- df_test$status <- df_test$taxonRankID <- df_test$isBrackish <- df_test$isFreshwater <- df_test$isTerrestrial <- df_test$match_type <- NULL
#check if it identified everything
dim(subset(df_test, is.na(df_test$phylum))) # ok

# In the class column, we only keep the 5 groups we want. 
df_test <- subset(df_test, class %in% c("Elasmobranchii","Actinopteri","Holocephali",
                                        "Myxini","Petromyzonti","Cephalopoda")) 
# get all cephalopoda
ceph <- subset(df_test,df_test$class == "Cephalopoda")
yy <- subset(xx, AphiaID %in% survey$AphiaID ) # Get cephalopoda without length measurements


keep_sp <- data.frame(df_test) # subsetting
keep_sp <- data.frame(unlist(keep_sp$valid_name)) #unlisting
names(keep_sp) <- 'ScientificName'
keep_ap <- data.frame(df_test) # subsetting
keep_ap <- data.frame(unlist(keep_ap$AphiaID))
names(keep_ap) <- 'AphiaID'
keep_gen <- data.frame(df_test) # subsetting
keep_gen <- data.frame(unlist(keep_gen$genus))
names(keep_gen) <- 'Genus'
keep_fa <- data.frame(df_test) # subsetting
keep_fa <- data.frame(unlist(keep_fa$family))
names(keep_fa) <- 'Family'
keep <- cbind(keep_ap, keep_sp, keep_gen, keep_fa)

dat.ices <- subset(dat.ices, dat.ices$AphiaID %in% keep_ap$AphiaID)
dat.ices <- left_join(dat.ices, keep, by='AphiaID')
dat.ices$Species <- dat.ices$ScientificName
dat.ices$ScientificName <- NULL
survey <- dat.ices

survey <- survey %>%
  dplyr::select(Survey, HaulID, StatRec, Year, Month, Quarter, Season, ShootLat, ShootLong, HaulDur, Area.swept, Area.doors, 
         Gear, Depth, SBT, SST,AphiaID, Family, Genus, Species, CatIdentifier, numcpue, wtcpue, numh, wgth, num, wgt, Length, numlencpue, numlenh)

### Code to integrate from Anna Rindorf on species bycatch corrections
survey <- data.frame(survey)
survey <- survey %>%
  mutate(Species = recode(Species,'Dipturus batis'='Dipturus','Dipturus flossada'='Dipturus',
                          'Dipturus batis-complex'='Dipturus','Dipturus intermedia'='Dipturus',
                          'Dipturus'='Dipturus','Liparis montagui'='Liparis',
                          'Liparis liparis'='Liparis','Liparis liparis liparis'='Liparis',
                          'Chelon aurata'='Chelon','Chelon ramada'='Chelon',
                          'Mustelus mustelus/asterias'='Mustelus','Mustelus'='Mustelus',
                          'Mustelus mustelus'='Mustelus','Mustelus asterias'='Mustelus',
                          'Alosa'='Alosa','Alosa alosa'='Alosa','Alosa fallax'='Alosa',
                          'Argentina'='Argentina','Argentinidae'='Argentina',
                          'Argentina silus'='Argentina','Argentina sphyraena'='Argentina',
                          'Callionymus reticulatus'='Callionymus','Callionymus maculatus'='Callionymus',
                          'Ciliata mustela'='Ciliata','Ciliata septentrionalis'='Ciliata',
                          'Gaidropsarus'='Gaidropsarus','Gaidropsaurus macrophthalmus'='Gaidropsarus',
                          'Gaidropsaurus mediterraneus'='Gaidropsarus','Gaidropsaurus vulgaris'='Gaidropsarus',
                          'Sebastes'='Sebastes','Sebastes norvegicus'='Sebastes','Sebastes mentella'='Sebastes',
                          'Sebastes marinus'='Sebastes','Syngnathus'='Syngnatus',
                          'Syngnathus rostellatus'='Syngnatus','Syngnathus acus'='Syngnatus',
                          'Syngnathus typhle'='Syngnatus','Nerophis ophidion'='Syngnatus',
                          'Pomatoschistus'='Pomatoschistus','Pomatoschistus microps'='Pomatoschistus',
                          'Pomatoschistus minutus'='Pomatoschistus','Pomatoschistus pictus'='Pomatoschistus',
                          'Lesueurigobius'='Gobius','Gobius cobitis'='Gobius','Gobius niger'='Gobius',
                          'Leusueurigobius friesii'='Gobius','Neogobius melanostomus'='Gobius',
                          'Neogobius'='Gobius'))


# ------------------------------------------------------------------------------
# CHECK % OF CEPHALOPODA ENTRIES LACKING LENGTH MEASUREMENT 
# ------------------------------------------------------------------------------
survey_ceph <- subset(survey, AphiaID %in% ceph$AphiaID)
100 * nrow(yy) / (nrow(yy) + nrow(survey_ceph))

# ------------------------------------------------------------------------------
# ADD GEAR EFFICIENCY CORRECTIONS 
# ------------------------------------------------------------------------------
#source('scripts - data processing/Merge_names_walker-catchability_with_datras.R')
load("traits and species/Names_DATRAS_Walker_match.Rdata")

# load all aphiaID of cephalopods --> Dani this can be updated if there is better information
ceph_q <- data.frame(AphiaID = ceph$AphiaID,Species=NA,Genus=NA,Family=NA,q_group = "GRP4") 
q_names <- rbind(q_names,ceph_q)

conversions <- read.csv("data/Walkeretal_2017_supp/EfficiencyTab.csv",header=T,sep=",")
conversions <- subset(conversions,conversions$Gear =="GOV")

# combine q_group/q_species code with survey 
survey <- cbind(survey,q_names[match(survey$AphiaID,q_names$AphiaID),c("q_group")])
colnames(survey)[ncol(survey)] <- "Code"

# get unique length (q's are length based) and code
uni <- survey %>%
  distinct(Code, Length)

uni <- left_join(uni, conversions, by = "Code") %>%
          mutate(Diff = abs(Length.x + 0.02 -Length.y)) %>%   # added 0.02 to avoid ending in the middle
          group_by(Length.x,Code) %>%
          filter(Diff == min(Diff)) 
uni$uni <- paste(uni$Code,uni$Length.x)

survey$uni <- paste(survey$Code,survey$Length)

# and combine to get efficiency
survey <- cbind(survey,uni[match(survey$uni,uni$uni),c("Efficiency")])
colnames(survey)[ncol(survey)] <- "q_eff"  
survey$q_eff <- ifelse(survey$q_eff<0.01, 0.01,survey$q_eff) # to avoid too high corrections at the boundaries of the GAM model (Walker) -- check sensitivity!!

survey <- survey %>%
  mutate(numlencpue_q    = numlencpue/q_eff,
         numlenh_q       = numlenh/q_eff) %>%
  select(Survey, HaulID, StatRec, Year, Month, Quarter, Season, ShootLat, ShootLong, HaulDur, Area.swept, Area.doors, 
         Gear, Depth, SBT, SST, AphiaID, Family, Genus, Species, CatIdentifier, numcpue, wtcpue, numh, wgth, num, wgt, Length, numlencpue, numlenh,
         numlencpue_q,numlenh_q)

# ------------------------------------------------------------------------------
# RE-CALCULATE WEIGHTS
# ------------------------------------------------------------------------------
detach(package:worms)
detach(package:plyr)

# 1. associate an LME to each haul and make final list of species
### Prepare list for estimating length-weight parameters
list.taxa <- survey %>% 
  select(HaulID, Survey, ShootLat, ShootLong, Family, Genus, Species) %>% 
  distinct()

# get LME
library(rgdal)
shape1 <- readOGR(dsn = "data/LME shapefile",layer="lme66")
coords <- list.taxa %>%
  dplyr::select(ShootLat, ShootLong, Survey) %>%
  distinct()
str(coords)

coordinates(coords) <- ~ ShootLong + ShootLat
proj4string(coords) <- proj4string(shape1)
lme <- over(coords, shape1)

coords <- list.taxa %>%
  dplyr::select(ShootLat, ShootLong, Survey) %>%
  distinct()
coords <- cbind(coords, lme$LME_NUMBER)
setnames(coords, old='lme$LME_NUMBER', new='lme')

# plot
# ggplot(coords,aes(ShootLong,ShootLat))+
#   borders('world', xlim=c(-90,50), ylim=c(25,85), fill='black', colour='black') +
#   coord_quickmap(xlim=c(-90,50), ylim=c(25,85))+
#   geom_polygon(data=shape1, aes(y=lat, x=long, group=group), fill='lightgrey', colour='black')+
#   theme_bw()+xlab('')+ylab('')+
#   geom_point(cex = 0.2, colour='blue')

coords$lme <- as.factor(coords$lme)
#Select from each LME 50 long and lat
ind <- c()
for (i in 1:nlevels(coords$lme)){
  ind <- c(ind, sample(which(coords$lme==levels(coords$lme)[i]), 50, replace = FALSE))
}
long50 <- coords$ShootLong[ind]
lat50 <- coords$ShootLat[ind]
lme50 <- rep(levels(coords$lme), each=50)
#For each haul without LME find a close LME that has an LME number already
nlme <- subset(coords, is.na(lme)) # many hauls without LME 710
nlme$ShootLat <- as.numeric(as.vector(nlme$ShootLat))
nlme$ShootLong <- as.numeric(as.vector(nlme$ShootLong))
long50 <- as.numeric(as.vector(long50))
lat50 <- as.numeric(as.vector(lat50))
dilme <- c()
for (i in 1:length(lme50)){
  dilme <- cbind(dilme, (nlme$ShootLat-lat50[i])**2 + (nlme$ShootLong-long50[i])**2)
}
mindi <- apply(dilme, 1, which.min) 
coords$lme[is.na(coords$lme)] <- lme50[mindi] # assign the closest LME number to each haul without LME

#Check
coords$ShootLat <- as.numeric(as.vector(coords$ShootLat))
coords$ShootLong <- as.numeric(as.vector(coords$ShootLong))

# rockall not assigned to Faroe plateau but to celtic sea LME
coords$lme <- as.character(coords$lme)
coords <- coords %>%
  mutate(lme = replace(lme, Survey  =='ROCKALL', '60')) %>%
  as.data.frame()
#plot(coords$ShootLong, coords$ShootLat, col=rainbow(length(unique(coords$lme)))[as.factor(coords$lme)], pch=".")

survey <- left_join(survey, coords, by=c('ShootLat', 'ShootLong','Survey')) 
survey <- survey %>% filter(Species!='Gobioidei')

# this was run before cephalopods - it is used to get a LxW relationship for 
# each LME; this information has not been updated for the cephalopod code
# a fixed lenght x weight is used for all cephalopods
# library(tidyverse)
# list.taxa <- survey %>% 
#  select(Family, Genus, Species, lme) %>% 
# distinct() %>%
#  mutate(fao = 27,
#         Subspecies = str_split(Species, pattern = " ", simplify=T)[,3],
#         Species = str_split(Species, pattern = " ", simplify=T)[,2], 
#        Species = if_else(Subspecies!="", paste(Species, Subspecies, sep=" "), Species))
#write.csv(data.frame(list.taxa), file="traits and species/taxa.DATRAS.FB.tofill5.csv", row.names=FALSE)
#save(survey, file = "data/DATRAS_before_lw_xxxxx.RData") # could save intermediate stage

# 2. re-calculate weights with length-weight relationships
datalw <- read.csv('traits and species/taxa.DATRAS.FB_filled5.csv') %>% 
  mutate(Taxon = case_when(level=='family' ~ family,
                           level=='genus' ~ genus,
                           level=='species' ~ paste(genus, species, sep=" ")),
         lme = as.factor(lme)) %>% 
  select(-fao,-family,-genus,-species)

ceph_q <- cbind(ceph_q,survey[match(ceph_q$AphiaID,survey$AphiaID),c("Species")])
colnames(ceph_q)[ncol(ceph_q)] <- "Taxon"
datalw_ceph <- data.frame(lme=1000,level="species",FB_E_Code=NA,source=NA,type.length=NA,
                          taxo=NA,b=3,a=0.01,Taxon=ceph_q$Taxon)
datalw_ceph$lme <- as.factor(datalw_ceph$lme)
datalw <- rbind(as.data.frame(datalw),datalw_ceph)

survey <- survey %>% 
  mutate(Taxon = case_when(is.na(Species) & is.na(Genus) ~ Family,
                           Species=="" & is.na(Genus) ~ Family,
                           is.na(Species) & !is.na(Genus) ~ Genus,
                           Species=="" & !is.na(Genus) ~ Genus,
                           !is.na(Species) ~ Species))

survey[survey$Species=='Syngnatus',]$Taxon <- 'Syngnathus'
survey[survey$Species=='Syngnatus',]$Species <- 'Syngnathus'

# now change lme for each cephalopoda 
survey$lme[survey$AphiaID %in% ceph_q$AphiaID] <- 1000

# summarize abundance/weight at the haul level
survey.num <- left_join(survey, datalw, by=c('Taxon','lme')) %>% 
  select(Survey,HaulID,StatRec,Year,Month,Quarter,Season,ShootLat,ShootLong,HaulDur,Area.swept,Area.doors,Gear,Depth,SBT,SST,Family,Genus,Species,Taxon,
         CatIdentifier,numcpue,numh,num) %>% 
  distinct() %>% 
  group_by(Survey,HaulID,StatRec,Year,Month,Quarter,Season,ShootLat,ShootLong,HaulDur,Area.swept,Area.doors,Gear,Depth,SBT,SST,Family,Genus,Species,Taxon) %>%
  summarize_at(.vars=c('numcpue', 'numh', 'num'), .funs = function(x) sum(x)) %>% 
  ungroup()

survey.wgt <- left_join(survey, datalw, by=c('Taxon','lme')) %>% 
  select(Survey,HaulID,StatRec,Year,Month,Quarter,Season,ShootLat,ShootLong,HaulDur,Area.swept,Area.doors,Gear,Depth,SBT,SST,Family,Genus,Species,Taxon,
         CatIdentifier,wtcpue,wgth,wgt) %>% 
  distinct() %>% 
  group_by(Survey,HaulID,StatRec,Year,Month,Quarter,Season,ShootLat,ShootLong,HaulDur,Area.swept,Area.doors,Gear,Depth,SBT,SST,Family,Genus,Species,Taxon) %>%
  summarize_at(.vars=c('wtcpue', 'wgth', 'wgt'), .funs = function(x) sum(x)) %>% 
  ungroup()

survey1 <- full_join(survey.num, survey.wgt, 
                     by=c('Survey','HaulID','StatRec','Year','Month','Quarter',
                          'Season','ShootLat','ShootLong','HaulDur','Area.swept','Area.doors',
                          'Gear','Depth','SBT','SST','Family','Genus','Species','Taxon'))

# summarize abundance/weight from length data
survey2 <- left_join(survey, datalw, by=c('Taxon','lme')) %>% 
  mutate(wgtlencpue = numlencpue*a*Length^b/1000, # divide by 1000 to get kg/km2
         wgtlenh = numlenh*a*Length^b/1000, # divide by 1000 to get kg/h
         wgtlencpue_q = numlencpue_q*a*Length^b/1000, # divide by 1000 to get kg/km2
         wgtlenh_q = numlenh_q*a*Length^b/1000) %>% # divide by 1000 to get kg/h
  group_by(Survey,HaulID,StatRec,Year,Month,Quarter,Season,ShootLat,ShootLong,HaulDur,Area.swept,Area.doors,Gear,Depth,SBT,SST,Family,Genus,Species,Taxon, a, b) %>% 
  summarize_at(.vars=c('numlencpue','numlenh','wgtlencpue','wgtlenh','numlencpue_q','numlenh_q','wgtlencpue_q','wgtlenh_q'), .funs=function(x) sum(x)) %>% 
  ungroup()

# merge both and compare
nrow(survey1)==nrow(survey2)
survey3 <- full_join(survey1, survey2, by=c('Survey','HaulID','StatRec','Year','Month','Quarter','Season','ShootLat','ShootLong','HaulDur',
                                             'Area.swept','Area.doors','Gear','Depth','SBT','SST','Family','Genus','Species','Taxon'))

library(ggplot2)

# correlation between abundances to check calculations are right
cor(x = survey3$numh, y = survey3$numlenh, method = 'pearson')
xx <- subset(survey3, !is.na(numcpue))
cor(x = xx$numcpue, y = xx$numlencpue, method = 'pearson')

# check weights
xx <- subset(survey3, wtcpue   >0 & wgtlencpue>0)
cor(x = xx$wtcpue , y = xx$wgtlencpue, method = 'pearson')

xx <- subset(survey3, wgth>0 & wgtlenh>0)
cor(x = xx$wgth, y = xx$wgtlenh, method = 'pearson')

### cor = 0.88 and 0.89 -> check per survey.

# ------------------------------------------------------------------------------
# CHECK PER SURVEY
# ------------------------------------------------------------------------------

# no zeros
xx <- subset(survey3, wgth>0 & wgtlenh>0)

# rockall looks OK
ggplot(subset(xx, Survey=='ROCKALL'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10() 

# IE-IGFS looks OK
ggplot(subset(xx, Survey=='IE-IGFS'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10() 

# NIGFS looks OK
ggplot(subset(xx, Survey=='NIGFS'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10() 

# PT-IBTS looks OK
ggplot(subset(xx, Survey=='PT-IBTS'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10() 

# FR-CGFS looks OK
ggplot(subset(xx, Survey=='FR-CGFS'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10() 

# Can-MARS looks OK
ggplot(subset(xx, Survey=='Can-Mar'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10() 

# SWC-IBTS issue
ggplot(subset(xx, Survey=='SWC-IBTS'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10()

comp <- subset(xx, Survey=='SWC-IBTS') %>% 
  select(HaulID,wgtlenh,wgth) %>% 
  distinct() %>% 
  group_by(HaulID) %>%
  summarize_at(.vars=c('wgtlenh', 'wgth'), .funs = function(x) sum(x)) %>% 
  ungroup() %>%
  as.data.frame()

ggplot(comp, aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10()

comp$factor <-   comp$wgtlenh / comp$wgth 
plot(comp$factor)
resc <- comp$HaulID[comp$factor > 40]

# after check with original haul length data (HL) for some resc haulid, weight is clearly wrong factor 100  
survey3 <- survey3 %>%
  mutate(wtcpue = if_else(HaulID %in% resc, wtcpue*100,wtcpue),
         wgth = if_else(HaulID %in% resc , wgth*100,wgth),
         wgt = if_else(HaulID %in% resc , wgt*100,wgt))

# BITS issue
ggplot(subset(xx, Survey=='BITS'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10() 

comp <- subset(xx, Survey=='BITS') %>% 
  select(HaulID,wgtlenh,wgth) %>% 
  distinct() %>% 
  group_by(HaulID) %>%
  summarize_at(.vars=c('wgtlenh', 'wgth'), .funs = function(x) sum(x)) %>% 
  ungroup() %>%
  as.data.frame()

ggplot(comp, aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10()

comp$factor <-   comp$wgtlenh / comp$wgth
plot(comp$factor)
resc <- comp$HaulID[comp$factor > 40]

# after check with original haul length data (HL) for some resc haulid, weight is clearly wrong factor 100  
survey3 <- survey3 %>%
  mutate(wtcpue = if_else(HaulID %in% resc, wtcpue*100,wtcpue),
         wgth = if_else(HaulID %in% resc , wgth*100,wgth),
         wgt = if_else(HaulID %in% resc , wgt*100,wgt))

# EVHOE may have an issue, no changes as not very clear
ggplot(subset(xx, Survey=='EVHOE'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10()

comp <- subset(xx, Survey=='EVHOE') %>% 
  select(HaulID,wgtlenh,wgth) %>% 
  distinct() %>% 
  group_by(HaulID) %>%
  summarize_at(.vars=c('wgtlenh', 'wgth'), .funs = function(x) sum(x)) %>% 
  ungroup() %>%
  as.data.frame()

ggplot(comp, aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10()

comp$factor <-   comp$wgtlenh / comp$wgth 
plot(comp$factor)

# NS - IBTS issue
ggplot(subset(xx, Survey=='NS-IBTS'), aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10() 

comp <- subset(xx, Survey=='NS-IBTS') %>% 
  select(HaulID,wgtlenh,wgth) %>% 
  distinct() %>% 
  group_by(HaulID) %>%
  summarize_at(.vars=c('wgtlenh', 'wgth'), .funs = function(x) sum(x)) %>% 
  ungroup() %>%
  as.data.frame()

ggplot(comp, aes(x=wgth, y=wgtlenh)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10()

comp$factor <-   comp$wgtlenh / comp$wgth
comp$uni <- c(1:nrow(comp))
plot(comp$factor~comp$uni,ylim=c(0,120))
points(comp$factor[comp$factor > 20]~comp$uni[comp$factor > 20],col="red")
points(comp$factor[comp$factor > 8 & comp$factor <20]~comp$uni[comp$factor  > 8 & comp$factor <20],col="blue")

# two issues - one estimate 100 times higher based on length, the other 10 times
resc <- comp$HaulID[comp$factor > 20] 
resc2 <- comp$HaulID[comp$factor > 8 & comp$factor <20]

# after check with original haul length data (HL) for some resc haulid, weight is clearly wrong factor 100 
# and also a cluster of factor 10
survey3 <- survey3 %>%
  mutate(wtcpue = if_else(HaulID %in% resc, wtcpue*100,wtcpue),
         wgth = if_else(HaulID %in% resc , wgth*100,wgth),
         wgt = if_else(HaulID %in% resc , wgt*100,wgt))

survey3 <- survey3 %>%
  mutate(wtcpue = if_else(HaulID %in% resc2, wtcpue*10,wtcpue),
         wgth = if_else(HaulID %in% resc2 , wgth*10,wgth),
         wgt = if_else(HaulID %in% resc2 , wgt*10,wgt))

# check again correlations
xx <- subset(survey3, wtcpue> 0 & wgtlencpue>0)
cor(x = xx$wtcpue , y = xx$wgtlencpue, method = 'pearson') # is (a bit) better

xx <- subset(survey3, wgth>0 & wgtlenh>0)
cor(x = xx$wgth, y = xx$wgtlenh, method = 'pearson') # is (a bit) better

# now check per haul without zeros, NAs
xx <- subset(survey3, wtcpue>0 & wgtlencpue>0)

comp <- xx %>% 
  select(HaulID,wgtlencpue,wtcpue) %>% 
  distinct() %>% 
  group_by(HaulID) %>%
  summarize_at(.vars=c('wgtlencpue', 'wtcpue'), .funs = function(x) sum(x)) %>% 
  ungroup() %>%
  as.data.frame()

ggplot(comp, aes(x=wtcpue, y=wgtlencpue)) + geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", size=0.5) + scale_x_log10() + scale_y_log10()

cor(x = xx$wtcpue , y = xx$wgtlencpue, method = 'pearson')
# [1] 0.9346

# ------------------------------------------------------------------------------
# SAVE DATA
# ------------------------------------------------------------------------------
survey3 <- survey3 %>% 
  select(-num, -wgt) %>%
  as.data.frame()

survey3 <- cbind(survey3,survey[match(survey3$Taxon,survey$Taxon),c("AphiaID")])
colnames(survey3)[ncol(survey3)] <- "AphiaID"

survey3 <- cbind(survey3,my_sp_taxo[match(survey3$AphiaID,my_sp_taxo$AphiaID),c("class")])
colnames(survey3)[ncol(survey3)] <- "class"

save(survey3, file='data/ICESsurveys30Nov_withq_ceph.RData')
