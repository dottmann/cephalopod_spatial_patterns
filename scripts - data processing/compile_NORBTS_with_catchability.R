#############################################################
#### Code to download and clean survey data from IMR survey
#### Coding: Aurore Maureaud, March 2021

#### 
#### Coding updated by Daniel van Denderen:
#### added wingspread / catch per unit area information
#### added gear efficiency estimates based on Walker et al. 2017
#### see https://doi.org/10.1093/icesjms/fsw250
#############################################################

rm(list=ls())

### Libraries
library(data.table)
library(dplyr)

### Load files
# First setting the working directory to the main folder with the unzipped csv files
data_dir <- paste(getwd(),"/data/data Norway/",sep="")

# we create a vector list with the filenames that match with a .csv ending
files = list.files(data_dir,pattern="*.csv")

# then we call a lapply function that takes x (every csv) and calls it back to a rbind. Check the seperator to see if it's correct
norw_dat = do.call(rbind, lapply(files, function(x) read.csv(paste(data_dir,x,sep=''), stringsAsFactors = FALSE, header = TRUE, sep = ";")))
rm(files)

# change colnames from Norwegian to new names in English
setnames(norw_dat, old = c("aar","mnd","lengde","bredde","redskap","starttid","stopptid","taueT","bunndyp",
                           "opening","dist","tilstand","kvalitet","delnr","akode","art","latin","maal_fangst","fangstKvant",
                           "fangstAnt","maal_lprov","lengdemaal","lengdeProveKv","lengdeProveAnt","interv","kjonn"), 
         new = c("Year","Month","ShootLong","ShootLat","Gear","ShootTimeB","ShootTimeE","HaulDur",
                 "Depth","Netopening","Distance","quality_gear","quality_haul","SubSampleNr",
                 "SpecCode","AkodeName","ScientificName","MeasureType","Weight","NoMeas","MeasureType2","LengthMethod",
                 "WeightSubSample","AbundanceSubSample","Interv","Sex"))


##########################################################################################
#### CHANGE SPECIES LATIN NAMES
##########################################################################################

# make sure our species names start with a capital letter, and then only lower case
#str(norw_dat)

# First we steal a function from the interwebs
capitalize <- function(x){
  first <- toupper(substr(x, start=1, stop=1)) ## capitalize first letter
  rest <- tolower(substr(x, start=2, stop=nchar(x)))   ## everything else lowercase
  paste0(first, rest)
}

# It requires the column to be factors
norw_dat$ScientificName <- as.factor(norw_dat$ScientificName)

# and then we run the function on it
levels(norw_dat$ScientificName) <- capitalize(levels(norw_dat$ScientificName))
rm(capitalize)


##########################################################################################
#### CREATE HAULD ID
##########################################################################################

# Give survey name
norw_dat$Survey <- rep("NorBTS",each=length(unique(rownames(norw_dat))))

# Haulid
norw_dat$HaulID <- paste(norw_dat$Survey, norw_dat$Year,norw_dat$Month,norw_dat$Gear,norw_dat$ShootLong,norw_dat$ShootLat, norw_dat$Depth, norw_dat$ShootTimeB)

# Recalculate the haul duration because the column has weird values
# start time: ShootTimeB in XYZW where XY are hours from 0 to 24 and ZW are minutes from 0 to 59
# end time: ShootTimeE
norw_dat[norw_dat$ShootTimeE==-1,]$ShootTimeE <- 'NA'
norw_dat[norw_dat$ShootTimeB==-1,]$ShootTimeB <- 'NA'
norw_dat$ShootTimeB <- as.numeric(as.vector(norw_dat$ShootTimeB))
norw_dat$ShootTimeE <- as.numeric(as.vector(norw_dat$ShootTimeE))

times <- data.frame(cbind(norw_dat$HaulID, norw_dat$ShootTimeB, norw_dat$ShootTimeE))
names(times) <- c('HaulID','ShootTimeB','ShootTimeE')
times <- subset(times, !is.na(times$ShootTimeB))
times <- subset(times, !is.na(times$ShootTimeE))
for(i in 1:ncol(times)){times[,i] <- as.character(times[,i])}
# add 0 as characters to have length 4 of times
times[nchar(times$ShootTimeB)==2,]$ShootTimeB <- paste('00',times[nchar(times$ShootTimeB)==2,]$ShootTimeB, sep='')
times[nchar(times$ShootTimeB)==3,]$ShootTimeB <- paste('0',times[nchar(times$ShootTimeB)==3,]$ShootTimeB, sep='')
times[nchar(times$ShootTimeB)==1,]$ShootTimeB <- paste('000',times[nchar(times$ShootTimeB)==1,]$ShootTimeB, sep='')
times[nchar(times$ShootTimeE)==2,]$ShootTimeE <- paste('00',times[nchar(times$ShootTimeE)==2,]$ShootTimeE, sep='')
times[nchar(times$ShootTimeE)==3,]$ShootTimeE <- paste('0',times[nchar(times$ShootTimeE)==3,]$ShootTimeE, sep='')
times[nchar(times$ShootTimeE)==1,]$ShootTimeE <- paste('000',times[nchar(times$ShootTimeE)==1,]$ShootTimeE, sep='')

# count minutes and hours for begining and end
times$minB <- as.numeric(as.vector(substr(times$ShootTimeB, start=1, stop=2)))*60+as.numeric(as.vector(substr(times$ShootTimeB, start=3, stop=4))) 
times$minE <- as.numeric(as.vector(substr(times$ShootTimeE, start=1, stop=2)))*60+as.numeric(as.vector(substr(times$ShootTimeE, start=3, stop=4))) 
times$duration <- times$minE-times$minB
times[times$minB>1320 & times$minE<120,]$duration <- times[times$minB>1320 & times$minE<120,]$minE-times[times$minB>1320 & times$minE<120,]$minB+1440
times[times$minB>1080 & times$minE<420,]$duration <- times[times$minB>1080 & times$minE<420,]$minE-times[times$minB>1080 & times$minE<420,]$minB+1440
# all remaining times are too long or start before begining time -> to be removed
times <- subset(times, times$duration>0)
# let's check the very high times: higher than 8h?
times <- unique(times)
times$ShootTimeB <- times$ShootTimeE <- times$minB <- times$minE <- NULL
setnames(times, old='duration', new='HaulDur2')

# join back with norw_dat
norw_dat0 <- left_join(norw_dat, times, by='HaulID')
nrow(norw_dat)==nrow(norw_dat0)
norw_dat <- norw_dat0


##########################################################################################
#### SELECT GEAR TYPES
##########################################################################################

# Remove all hauls done with "wrong" gear types. Keeping "torske", "reke", "benthos
removed_gear <- c("3113","3114","3115","3118", # konsumtål
                  "3119", # sestad trål
                  "3130","3173","3174","3175","3176","3177", # industritrål
                  "3131","3171","3172", # tobistrål
                  "3132","3133","3134", # single, dobbelt and trippel
                  "3136", #bomtrål
                  "3279", #fisketrål
                  "3400", #trål
                  "3401", #IKMT
                  "3410", "3411", "3412", #semipelagisk
                  "3415", # partrål
                  "3420","3421","3422", "3423", #krabbetrål
                  "3430", #plankton
                  "3440" #bomtrål
)
norw_dat <- norw_dat[!norw_dat$Gear %in% removed_gear,]
rm(removed_gear)


##########################################################################################
#### REMOVE BAD QUALITY HAULS
##########################################################################################
# Remove bad quality hauls and gears
norw_dat <- subset(norw_dat, norw_dat$quality_gear %in% c(1,2))
norw_dat <- subset(norw_dat, norw_dat$quality_haul %in% c(1,2))

# Is there still empty species names and abundances?
check.sp <- subset(norw_dat, norw_dat$ScientificName=='') # all hauls from 1981 and 1982 with no ab/weight/spp specified
norw_dat <- subset(norw_dat, norw_dat$ScientificName!='') # remove rows with empty rows

check.ab <- subset(norw_dat, is.na(norw_dat$NoMeas)) # ok
check.sub.ab <- subset(norw_dat, is.na(norw_dat$AbundanceSubSample))
check.sum <- subset(norw_dat, is.na(norw_dat$Sum))
check.sub.w <- subset(norw_dat, is.na(norw_dat$WeightSubSample)) # same as abundance


##########################################################################################
#### STANDARDIZE UNITS AND REMOVE NEGATIVE VALUES
##########################################################################################

# HaulDuration: if the range 1-60m then minutes. If 0-1, in hours
# ICES data in minutes, convert all in minutes 1h <-> 60min
# -1, data unavailable, so insert NA

norw_dat[norw_dat$HaulDur<=1,]$HaulDur <- norw_dat[norw_dat$HaulDur<=1,]$HaulDur*60
norw_dat[norw_dat$HaulDur<10 & norw_dat$Distance>2,]$HaulDur <- norw_dat[norw_dat$HaulDur<10 & norw_dat$Distance>2,]$HaulDur*60
norw_dat[norw_dat$HaulDur<0,]$HaulDur <- NA

# Transform distance nautical miles to km
# 1nm <-> 1.852km

norw_dat$Distance <- norw_dat$Distance*1.852/1
norw_dat[norw_dat$Distance<0,]$Distance <- NA

# remove strange distances based on speed and duration
norw_dat <- norw_dat %>%
 mutate(fact = Distance / (HaulDur/60*(3 * 1.825 )), # speed is around 3 knots (see Jakobsen et al. 1998 - ICES C.M. 1997N: 17)
        Distance = ifelse(fact <0.5 | fact >2, NA, Distance)) %>%
          select(-fact)
        
# Change net opening to DoorSpread
setnames(norw_dat, old = "Netopening", new="DoorSpread")
norw_dat[norw_dat$DoorSpread<0,]$DoorSpread <- NA
norw_dat$DoorSpread <- norw_dat$DoorSpread/1000 # transform m into km
norw_dat$WingSpread <- norw_dat$DoorSpread * 0.3

# Transform abundance and weight into the same units, transform weight measures all in kg
# for column Weight, use MeasureType
# for column NoMeas, use MeasureType as is category 6, *1000 individuals
# No document for conversion factors from L weight measurements!!!
# No liters measurements after 2001, so ok if we only select from 2005
# Two rows are in MeasureType or MeasureType==6, but in 1993 and 1995, so will be removed
norw_dat[norw_dat$MeasureType==5,]$Weight <- norw_dat[norw_dat$MeasureType==5,]$Weight*1 # changed from x 1000 correct??
norw_dat[norw_dat$MeasureType==6,]$Weight <- norw_dat[norw_dat$MeasureType==6,]$Weight*1000*1000
norw_dat[norw_dat$MeasureType==6,]$NoMeas <- norw_dat[norw_dat$MeasureType==6,]$NoMeas*1000
norw_dat[norw_dat$MeasureType==7,]$Weight <- norw_dat[norw_dat$MeasureType==7,]$Weight*1000
norw_dat[norw_dat$MeasureType==8,]$Weight <- norw_dat[norw_dat$MeasureType==8,]$Weight*1000
norw_dat[norw_dat$MeasureType==9,]$Weight <- norw_dat[norw_dat$MeasureType==9,]$Weight/1000

# Correction factors for gutted/without head and L transfo. might exist, but cannot find it

# Transform units from the sub-samples not possible because of NAs
norw_dat[is.na(norw_dat$WeightSubSample),]$WeightSubSample <- -1
norw_dat[is.na(norw_dat$AbundanceSubSample),]$AbundanceSubSample <- -1

norw_dat[norw_dat$MeasureType2==5,]$WeightSubSample <- norw_dat[norw_dat$MeasureType2==5,]$WeightSubSample*1000
norw_dat[norw_dat$MeasureType2==6,]$WeightSubSample <- norw_dat[norw_dat$MeasureType2==6,]$WeightSubSample*1000*1000
norw_dat[norw_dat$MeasureType2==6,]$AbundanceSubSample <- norw_dat[norw_dat$MeasureType2==6,]$AbundanceSubSample*1000
norw_dat[norw_dat$MeasureType2==7,]$WeightSubSample <- norw_dat[norw_dat$MeasureType2==7,]$WeightSubSample*1000
norw_dat[norw_dat$MeasureType2==8,]$WeightSubSample <- norw_dat[norw_dat$MeasureType2==8,]$WeightSubSample*1000
norw_dat[norw_dat$MeasureType2==9,]$WeightSubSample <- norw_dat[norw_dat$MeasureType2==9,]$WeightSubSample/1000

# Replace all -1 by NAs
norw_dat[norw_dat$WeightSubSample==(-1),]$WeightSubSample <- NA
norw_dat[norw_dat$AbundanceSubSample==(-1),]$AbundanceSubSample <- NA
norw_dat[norw_dat$Weight==(-1),]$Weight <- NA
norw_dat[norw_dat$NoMeas==(-1),]$NoMeas <- NA


##########################################################################################
#### CHANGE FORMAT AND AGGREGATE
##########################################################################################

library(reshape2)
library(dplyr)
# rehape format with length measurements and delete 0
norw_dat <- melt(norw_dat, c(names(norw_dat)[c(1:28)], names(norw_dat)[70:74]), c(29:69), variable.name='Length', value.name='NumLen')

sum.pos <- norw_dat %>% #data with abundance at length, 354061 unique speciesxHauldIDs, it works!
  filter(Sum>0) %>%
  mutate(#HaulSpe = paste(HaulID, ScientificName, sep=' '),
    Length=replace(Length, is.na(NumLen), NA)) %>%
  filter(!is.na(Length))

sum.na <- norw_dat %>%
  filter(is.na(Sum)) %>%
  mutate(Length=replace(Length, is.na(NumLen), NA)) %>% # give NA to Length when there is no abundance at length
  distinct() # remove duplicates without length compisition data

norw_dat <- rbind(sum.pos, sum.na)

# Estimate missing swept areas
norw_dat <- norw_dat %>%
  mutate(Area.doors = DoorSpread*Distance,
         Area.swept = WingSpread*Distance)

# set NA for one very high haul 
norw_dat <- norw_dat %>%
  mutate(Area.swept = ifelse(Area.swept > 2, NA, Area.swept),
         Area.doors = ifelse(Area.doors > 2/0.3, NA, Area.doors))

nor <- norw_dat %>%
  select(HaulID, Year, Area.swept, Area.doors, HaulDur, Gear, Depth, Distance) %>%
  filter(Year>1989,
         !is.na(HaulDur)) %>%
  distinct()

par(mfrow=c(1,2))
plot(Area.swept ~ HaulDur, data=nor)
plot(Area.swept ~ Depth, data=nor)

nor$Dur2 <- (nor$HaulDur-mean(nor$HaulDur))^2
lm0 <- lm(Area.swept ~ HaulDur + Dur2, data=nor)

pred0 <- predict(lm0, newdata=nor, interval='confidence', level=0.95)
nor <- cbind(nor,pred0)
nor[is.na(nor$Area.swept),]$Area.swept <- nor[is.na(nor$Area.swept),]$fit
nor[is.na(nor$Area.doors),]$Area.doors <- nor[is.na(nor$Area.doors),]$fit/0.3

nor <- nor %>%
  select(HaulID, Area.swept,Area.doors) %>%
  dplyr::rename(Area2=Area.swept, Door2 = Area.doors) %>%
  filter(Area2>=0)

nor2 <- left_join(norw_dat, nor, by='HaulID')
nor2 <- nor2 %>%
  mutate(Area.swept = coalesce(Area.swept,Area2),
         Area.doors = coalesce(Area.doors,Door2))
norw_dat <- nor2  

# Continue cleaning
norw_dat <- norw_dat %>%
  mutate(Quarter = ceiling(as.numeric(Month)/3),
         numcpue = NoMeas/Area.swept, # nbr / km2
         wtcpue = Weight/Area.swept, # kg / km2
         numh = NoMeas*60/HaulDur2, # nbr / hour
         wgth = Weight*60/HaulDur2, # kg / h
         SubFactor.Ab = NoMeas/Sum,
         numlencpue = NumLen*SubFactor.Ab/Area.swept, # raise abundance at length to the whole sample / swept area of haul Dur
         numlenh = NumLen*SubFactor.Ab*60/HaulDur2,
         Survey = 'NorBTS',
         Season = 'NA',
         SBT=NA, 
         SST=NA,
         HaulDur = HaulDur2,
         Species = ScientificName) %>%
  select(Survey, HaulID, Year, Month, Quarter, Season, ShootLat, ShootLong, HaulDur, Area.swept, Area.doors, Gear, Depth, SBT, SST, Species,
         numcpue, wtcpue, numh, wgth, Length, numlencpue, numlenh)


# Change lengths
norw_dat$Length <- as.factor(norw_dat$Length)
kk <- levels(norw_dat$Length)[1:40]
Min <- as.numeric(substr(kk, start=2, stop=4)) ## get lower interval value
Max <- as.numeric(substr(kk, start=6, stop=8))   ## get higher interval value
kkk <- cbind(Min,Max)
kkkk <- apply(kkk, 1, FUN=mean)
levels(norw_dat$Length) <- c(kkkk, 5500)

# To keep for th merge script
#norw_dat <- subset(norw_dat, !norw_dat$Distance==0)
# Remove haul duration lower than 20' and higher than 2h
norw_dat <- subset(norw_dat, norw_dat$HaulDur<120 & norw_dat$HaulDur>20)



##########################################################################################
#### CLEAN SPECIES NAMES
##########################################################################################

library(worms)
library(worrms)
library(crul)
library(urltools)

# First we create a species list with unique species names
sp_list <- as.data.frame(norw_dat$Species)
colnames(sp_list)[1] <- c("species") # changing the column name
sp_list <- data.frame(sp_list[!duplicated(sp_list[c("species")]),])
colnames(sp_list)[1] <- c("species") # changing the column name
sp_list <- subset(sp_list, !is.na(sp_list$species))

# Clean species names
sp_list$species <- as.character(sp_list$species, stringsAsFactors = FALSE)

cleanspl <- function(name){
  temp <- paste0(toupper(substr(name, 1, 1)), tolower(substr(name, 2, nchar(name))))
  temp <- gsub("\\((.+?)\\)$", "", temp)
  temp <- gsub("\\.$", "", temp) #remove last .
  temp <- gsub(" $", "", temp) #remove empty space
  temp <- gsub("            $", "", temp) #remove empty spaces
  temp <- gsub("^ ", "", temp)
  temp <- gsub(" .$", "", temp) #remove last single character
  temp <- gsub(" .$$", "", temp) #remove last double characters
  temp <- gsub(" unid$", "", temp)
  temp <- gsub(" unident$", "", temp)
  temp <- gsub(" spp$", "", temp)
  temp <- gsub(" spp.$", "", temp)
  temp <- gsub(" sp$", "", temp)
  temp <- gsub(" sp.$", "", temp)
  temp <- gsub(" s.c$", "", temp)
  temp <- gsub(" s.p$", "", temp)
  temp <- gsub(" s.f$", "", temp)
  temp <- gsub(" s.o$", "", temp)
  temp <- gsub(" so$", "", temp)
  temp <- gsub(" yoy$", "", temp)
  temp <- gsub(" YOY$", "", temp)
  temp <- gsub(" n. gen$", "", temp)
  temp <- gsub(" aff.$", "", temp)
  temp <- gsub(" kroyer", "", temp)
  temp <- gsub("-occidentale", "", temp)
  temp <- gsub("_rubricing", "", temp)
  temp <- gsub("-dilectus", "", temp)
  temp <- gsub(" cf tubicola", "", temp)
  temp <- gsub(",$", "", temp)
  temp <- gsub(",", " ", temp)
  return(temp)
}
sp_list$species.cleaned <- cleanspl(sp_list$species)
sp_list$species.cleaned <- cleanspl(sp_list$species.cleaned)
sp_list$speciesurl <- gsub(" ", "+",sp_list$species)

# Creating a for loop that check each species names with WoRMS and returns aphiaID to a new column in the species list
sp_list$aphiaID <- NA
for(i in 1:length(unique(sp_list$species.cleaned))){
  print(i)
  y <- sp_list[i,3]
  sp_list$aphiaID[i] <- tryCatch(as.data.frame(wm_name2id(name = y)), error=function(err) NA)
}
sp_list <- sp_list[,-3]

sp_list$aphiaID <- as.character(sp_list$aphiaID)
#setwd('C:/Users/auma/Documents/PhD DTU Aqua/(iv) Clean surveys/2. Clean taxonomy')
#write.csv(sp_list, file='sp_list_europe_2004_2017.csv', row.names=FALSE)
# save(sp_list, file = "traits and species/spp_norway.Rdata")
# load(file = "traits and species/spp_norway.Rdata")


# -------------------------DANI--------------------------------------------------
# the norway names file doesn't need to be saved each run
# not all cephalopods got identified due to missing name information
# we can manually update in the file: check.Norway.names.WORMS_squid.csv
# afterwards all species with an aphiaid are included 
# I made a first check - Dani you might want to see if I missed anything
# -------------------------DANI--------------------------------------------------

clean.names <- read.delim('traits and species/check.Norway.names.WORMS_ceph.csv', sep = ";", stringsAsFactors = T)
clean.names$aphiaID <- NULL
setnames(sp_list, old='species.cleaned', new='ScientificName')
# species = to match back with dataset
# scientific name = merge with worms correct names
# ScientificName_worms= correct names for some species (the ones that did not have the aphiaID)
sp_list <- left_join(sp_list, clean.names, by='ScientificName')
sp_list_ok <- subset(sp_list, sp_list$aphiaID!='NA')
sp_list_change <- subset(sp_list, sp_list$aphiaID=='NA')
sp_list_change$aphiaID <- sp_list_change$AphiaID_worms
sp_list_change$species.cleaned <- NULL
sp_list_change$ScientificName <- sp_list_change$ScientificName_worms

sp_list <- rbind(sp_list_ok, sp_list_change)
sp_list$AphiaID_worms <- sp_list$ScientificName_worms <- NULL
sp_list <- subset(sp_list, sp_list$ScientificName!='')
sp_list$aphiaID <- as.numeric(sp_list$aphiaID)

# The Aphia ID loop returns NAs and negative ID. Check, nrows has to be 0
sp.pb<-subset(sp_list, is.na(aphiaID) | aphiaID<0) 
dim(sp.pb) # good

# list with only aphia IDs
aphia_list <- c(sp_list$aphiaID) 

# remove duplicates
aphia_list <- aphia_list[!duplicated(aphia_list)]

# creating taxonomy tables for each species
# creating taxonomy tables for each species
spt <- c(seq(1,length(aphia_list),by=50),length(aphia_list))
my_sp_taxo <- c()
for (j in 1:(length(spt)-1)){
  spsub <- wm_record(aphia_list[spt[j]:spt[j+1]])
  my_sp_taxo <- rbind(my_sp_taxo,spsub)
}

# Remove duplicates:
my_sp_taxo <- my_sp_taxo[!duplicated(my_sp_taxo), ]

# Save the table:
# save(my_sp_taxo, file = "data/data Norway/ices_spp_taxo.Rdata")

# Load table:
load(file = "data/data Norway/ices_spp_taxo.Rdata")

# row binds all the results and pass to data frame. 
df_test <- data.frame(my_sp_taxo)
df_test$url <- df_test$lsid <- df_test$citation <- NULL
df_test$isExtinct <- df_test$modified <- df_test$valid_authority <- df_test$unacceptreason <- NULL
df_test$authority <- df_test$status <- df_test$taxonRankID <- df_test$isBrackish <- df_test$isFreshwater <- df_test$isTerrestrial <- df_test$match_type <- NULL
#check if it identified everything
dim(subset(df_test, is.na(df_test$phylum))) # ok

# In the class column, we only keep the 5 class we want, corresponding to fish species
df_test <- subset(df_test, class %in% c("Elasmobranchii","Actinopteri","Holocephali","Myxini","Petromyzonti","Cephalopoda", "Teleostei")) 

# ------------------ DANI ------------------------------------------------------
# get all cephalopoda - not used just so you can see what is included
ceph <- subset(df_test,df_test$class == "Cephalopoda")
write.table(ceph, "C:/Users/danot/My Drive/_postdoc/projects/cephalopods_WP1/traits data base/norway_cephalopoda_traits.txt", append = F, quote = F, sep = ";", row.names = F, col.names = T)
# ------------------ DANI -----------------------------------------------------

# List of names to keep
keep_sp <- data.frame(df_test) # subsetting
keep_sp <- data.frame(unlist(keep_sp$scientificname)) #unlisting
names(keep_sp) <- 'rightname'
sp_list <- subset(sp_list, ScientificName %in% keep_sp$rightname)
sp_list <- sp_list[!duplicated(sp_list$aphiaID),]
keep_sp$rightname <- as.character(keep_sp$rightname)
#identical(sort(sp_list$ScientificName), sort(keep_sp$rightname))

norw_dat <- subset(norw_dat, norw_dat$Species %in% keep_sp$rightname)
setnames(sp_list, old='species', new='Species')
norw_dat <- left_join(norw_dat, sp_list, 'Species')
norw_dat$AphiaID <- norw_dat$aphiaID
norw_dat$aphiaID <- NULL
norw_dat <- norw_dat %>% 
  mutate(Species = ScientificName)

# add class, family and genus
keep_sp <- data.frame(df_test) # subsetting
keep_sp <- data.frame(unlist(keep_sp$valid_name)) #unlisting
names(keep_sp) <- 'ScientificName'
keep_ap <- data.frame(df_test) # subsetting
keep_ap <- data.frame(unlist(keep_ap$AphiaID))
names(keep_ap) <- 'AphiaID'
keep_class <- data.frame(df_test) # subsetting
keep_class <- data.frame(unlist(keep_class$class))
names(keep_class) <- 'Class'
keep_gen <- data.frame(df_test) # subsetting
keep_gen <- data.frame(unlist(keep_gen$genus))
names(keep_gen) <- 'Genus'
keep_fa <- data.frame(df_test) # subsetting
keep_fa <- data.frame(unlist(keep_fa$family))
names(keep_fa) <- 'Family'
keep <- cbind(keep_ap, keep_sp, keep_gen, keep_fa,keep_class)

norw_dat <- cbind(norw_dat,keep[match(norw_dat$AphiaID,keep$AphiaID),c("Genus","Family","Class")])

norw_dat <- norw_dat %>%
  mutate(StatRec=NA) %>% 
  select(Survey, HaulID, StatRec, Year, Month, Quarter, Season, ShootLat, ShootLong, HaulDur, Area.swept,Area.doors, Gear, Depth, SBT, SST, Class, Family, 
         Genus, Species, AphiaID, numcpue, wtcpue, numh, wgth, Length, numlencpue, numlenh)

rm(df_test, keep_sp, keep,keep_class, keep_ap,keep_fa,keep_gen, my_sp_taxo, norw_dat0, sp_list, sp_list_change, sp_list_ok, sp.pb, check.ab, check.sp, check.sub.ab)
rm(check.sub.w, check.sum, clean.names, sum.na, sum.pos, times)
rm(i, aphia_list, kk, kkk, Max, Min, y, cleanspl, kkkk,pred0,lm0,nor,nor2,spsub,spt,j)

#save(norw_dat, file='cleaned data/NORBTSJuly2022_intermediate.RData')
##########################################################################################
#### RE-CALCULATE WEIGHTS based on L x W relationship
##########################################################################################
detach(package:worms)
detach(package:plyr)

# 1. associate an LME to each haul and make final list of species
### Prepare list for estimating length-weight parameters
list.taxa <- norw_dat %>% 
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

coords$lme <- as.character(coords$lme)
coords <- coords %>%
  mutate(lme = replace(lme, lme  == '60', '21'),
         lme = replace(lme, lme  == '19', '20')) %>%
  as.data.frame()

# now for each haul without LME select closest LME
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
nlme <- subset(coords, is.na(lme)) # hauls without LME 31
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

norw_dat <- left_join(norw_dat, coords, by=c('ShootLat', 'ShootLong','Survey')) 

library(tidyverse)
list.taxa <- norw_dat %>% 
  select(Family, Genus, Species, lme) %>% 
  distinct() %>%
  mutate(fao = 27,
         Subspecies = str_split(Species, pattern = " ", simplify=T)[,3],
         Species = str_split(Species, pattern = " ", simplify=T)[,2],
         Species = if_else(Subspecies!="", paste(Species, Subspecies, sep=" "), Species))
#write.csv(data.frame(list.taxa), file="traits and species/taxa.NorwDat.FB.tofill.csv", row.names=FALSE)

# 2. re-calculate weights with length-weight relationships
datalw <- read.csv('traits and species/taxa.NorwDat.FB_filled_EU.csv') %>% 
  mutate(Taxon = case_when(level=='family' ~ family,
                           level=='genus' ~ genus,
                           level=='species' ~ paste(genus, species, sep=" ")),
         lme = as.factor(lme)) %>% 
  select(-fao,-family,-genus,-species,-uni,-mergeLME)

norw_dat <- norw_dat %>% 
  mutate(Taxon = Species)

# get length weight fish
norw_dat <- left_join(norw_dat, datalw, by=c('Taxon','lme'))

# -------------------DANI-------------------------------------------------------------
# get length weight cephal - load from separate file
# Add cephalopod traits:
ceph_traits <- read.delim("traits and species/norway_cephalopoda_traits2.txt", sep = "\t", stringsAsFactors = T)
ceph_traits <- ceph_traits %>%
  dplyr::select(Taxon = scientificname,
                a = LW_a,
                b = LW_b)
ceph_taxa <- unique(ceph_traits$Taxon)

# Create "exclude" function:
"%ni%" <- Negate("%in%")

datalw_ceph <- norw_dat %>%
  filter(Taxon %in% ceph_taxa) %>%
  dplyr::select(-a, -b) %>%
  left_join(ceph_traits, by = "Taxon")

norw_dat <- norw_dat %>%
  filter(Taxon %ni% ceph_taxa) %>%
  bind_rows(datalw_ceph)
  
# -------------------DANI-----------------------------------------------------------

norw_dat$Length <- as.numeric(as.character(norw_dat$Length))
norw_dat$wgtlenh <- norw_dat$a*norw_dat$Length^norw_dat$b*norw_dat$numlenh/1000
norw_dat$wgtlencpue <- norw_dat$a*norw_dat$Length^norw_dat$b*norw_dat$numlencpue/1000

# get species catchabilities per length class
source('scripts - data processing/Merge_names_walker-catchability_with_Norway.R')
conversions <- read.csv("data/Walkeretal_2017_supp/EfficiencyTab.csv",header=T,sep=",")
conversions <- subset(conversions,conversions$Gear =="GOV")

# combine q_group/q_species code with norw_dat 
norw_dat <- cbind(norw_dat,datq[match(norw_dat$AphiaID,datq$AphiaID),c("q_group")])
colnames(norw_dat)[ncol(norw_dat)] <- "Code"

# -------------------DANI-------------------------------------------------------------
# add catchability group cephalopods
# --------------------------------------------------------------------------------
norw_dat$Code <- ifelse(norw_dat$Class == "Cephalopoda", "GRP4", norw_dat$Code)
# -------------------DANI-------------------------------------------------------------

# get unique length (q's are length based) and code
uni <- norw_dat %>%
  distinct(Code, Length)

uni <- left_join(uni, conversions, by = "Code") %>%
  mutate(Diff = abs(Length.x + 0.02 -Length.y)) %>%   # added 0.02 to avoid ending in the middle
  group_by(Length.x,Code) %>%
  filter(Diff == min(Diff)) 
uni$uni <- paste(uni$Code,uni$Length.x)

norw_dat$uni <- paste(norw_dat$Code,norw_dat$Length)

# and combine to get efficiency
norw_dat <- cbind(norw_dat,uni[match(norw_dat$uni,uni$uni),c("Efficiency")])
colnames(norw_dat)[ncol(norw_dat)] <- "q_eff"  
norw_dat$q_eff <- ifelse(norw_dat$q_eff<0.01, 0.01,norw_dat$q_eff) # to avoid too high corrections at the boundaries of the GAM model (Walker) -- check sensitivity!!

# estimate with correction of catchability
norw_dat$numlencpue_q  <- norw_dat$numlencpue/ norw_dat$q_eff
norw_dat$numlenh_q     <- norw_dat$numlenh/ norw_dat$q_eff
norw_dat$wgtlenh_q     <- norw_dat$wgtlenh/ norw_dat$q_eff
norw_dat$wgtlencpue_q  <- norw_dat$wgtlencpue/ norw_dat$q_eff

# Calculate mean q_eff of cephaloopods:
norw_dat %>%
  filter(Class == "Cephalopoda") %>%
  summarise(mean_efficiency_cephalopoda = mean(q_eff, na.rm = T),
            mean_efficiency_cpue_cephalopoda = sum(wtcpue * q_eff, na.rm = T) / sum(wtcpue, na.rm = T),
            mean_efficiency_length_cephalopoda = sum(Length * q_eff, na.rm = T) / sum(Length, na.rm = T))


# part of data is per species sub-group, part of data is per length-class
subgroup <- norw_dat[,c(1:25,40)]
subgroup <- subgroup %>%
  distinct() %>%
  group_by(Survey,HaulID,StatRec, Year,Month,Quarter,Season, ShootLat,ShootLong,HaulDur,Area.swept,Area.doors,Gear,Depth,SBT,SST,Class, Family, Genus, Species,Code)  %>% 
summarize_at(.vars=c('numcpue', 'wtcpue', 'numh','wgth'), .funs = function(x) sum(x)) %>% 
  as.data.frame()

lengthcl <- norw_dat[,c(2,20,27,28,36,37,38,39,43:46)]
lengthcl <- lengthcl %>%
  distinct() %>%
  group_by(HaulID,Species, a,b)  %>% 
  summarize_at(.vars=c('numlencpue', 'wgtlencpue', 'numlenh','wgtlenh','numlencpue_q', 'wgtlencpue_q', 'numlenh_q','wgtlenh_q'), .funs = function(x) sum(x)) %>% 
  as.data.frame()

norw_dat <- left_join(subgroup, lengthcl, by=c('HaulID','Species'))

# ----------------------------DANI ------------------------------------------------
# many cephalopod do not have a measured length observation and we can therefore not estimate 
# weight based on the length-weight relationship
# as a hack we can use the weight and multiply with the catchability factor 
# for all with NA
norw_dat$wgtlencpue_q <- ifelse(is.na(norw_dat$wgtlencpue_q) & norw_dat$Class == "Cephalopoda", norw_dat$wtcpue * 0.3, norw_dat$wgtlencpue_q)

# now reasonable spatial coverage of squid
# ----------------------------DANI ------------------------------------------------

# Save the data:
save(norw_dat, file='data/NORBTSdec2022_Ceph.RData')


# Try removing all hauls lacking at least 1 length measurement
lacking <- unique(subset(norw_dat, is.na(wgtlencpue_q))$HaulID)
temp <- norw_dat %>%
  filter(HaulID %ni% lacking)
# -----
# Less than half the observations are left.
# Ask Daniël about this


##########################################################################################
#### CHECKING SPATIAL DISTRIBUTION OF DATA POINTS
##########################################################################################

require(ggplot2)

coords <- norw_dat %>%
  # filter(Month %in% c(8,9)) %>%
  filter(ShootLat < 62, Class == "Cephalopoda") %>%
  select(ShootLat, ShootLong, wgtlencpue_q) %>%
  distinct()

ggplot(coords, aes(ShootLong,ShootLat))+
  #borders(xlim=c(-120,-110),ylim=c(40,41),fill="azure3",colour = "black") +
  borders(xlim=c(-20,50),ylim=c(54,82),fill="lightgrey",colour = "lightgrey") +
  coord_quickmap(xlim=c(-20,50),ylim=c(54,82))+theme_bw()+
  geom_point(aes(color = wgtlencpue_q), size = 2)


norw_dat %>%
  filter(Class == "Cephalopoda") %>%
  group_by(HaulID, ShootLat) %>%
  summarise(wgtlencpue_q = sum(wgtlencpue_q),
            numlenh = sum(numlenh)) %>%
  ggplot() +
  geom_point(aes(x = ShootLat, y = wgtlencpue_q)) +
  scale_y_log10()

coordY <- norw_dat %>%
  filter(Month %in% c(8,9)) %>%
  dplyr::select(HaulID, ShootLat, ShootLong, Survey, Year) %>%
  distinct() %>%
  group_by(Year) %>%
  summarize(number = length(unique(HaulID)))
