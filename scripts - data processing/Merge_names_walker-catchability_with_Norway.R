
# run to get species catchabilities per length class

load("traits and species/Names_DATRAS_Walker_match.Rdata")

datq <- norw_dat %>% 
  select(Family, Genus, Species, AphiaID) %>% 
  distinct() %>% 
  as.data.frame 

datq <- cbind(datq,q_names[match(datq$AphiaID,q_names$AphiaID),c("q_group")])
colnames(datq)[ncol(datq)] <- "q_group"

datq <- cbind(datq,q_names[match(datq$Genus,q_names$Genus),c("q_group")])
colnames(datq)[ncol(datq)] <- "q_group2"

datq$q_group <- ifelse(is.na(datq$q_group),datq$q_group2,datq$q_group)

# add five species mannually
datq$q_group <- ifelse(datq$AphiaID == 126433 , "GRP4" ,datq$q_group)
datq$q_group <- ifelse(datq$AphiaID == 127213 , "GRP7" ,datq$q_group)
datq$q_group <- ifelse(datq$AphiaID == 126432 , "GRP4" ,datq$q_group)
datq$q_group <- ifelse(datq$AphiaID == 127231 , "GRP7" ,datq$q_group)
datq$q_group <- ifelse(datq$AphiaID == 127119 , "GRP2" ,datq$q_group)
datq$q_group <- ifelse(is.na(datq$AphiaID) , NA ,datq$q_group)

