# --------------------------------------
# FUNCTION clean_ebird
# required packages: none
# description: Clean up e-Bird data downloads
#              keep only the following eBird fields to decrease object size
#              1. Unique identifier
#              7. Scientific name - later, list sensitive species and remove 
#                  from analysis. Also, beware of sp endings.
#              10. Exotic code - for removing escapees X)
#              16. Country
#              18. State
#              20. County
#              29. Latitude
#              30. Longitude
#              31. Observation date
#              32. Time observation started
#              35. Sampling event identifier
#              37. Protocol name - used to have "Protocol type"
#              41. Duration in minutes
#              42. Effort in km traveled
#              43. Effort in ha (area) covered
#              44. Number of observers
#              45. Where all species reported (0/1)
#              46. Group identifier
#              48. Approved (0/1) - check if everything that made it through 
#                  Vivi download is 1. If yes, remove this field.
# inputs: data.frame with eBird data, a focal species name, and a max route length
# outputs: list with data.frame of cleaned eBird data and vector with number
#          of observations kept after applying each cleaning criterion
########################################
clean_ebird <- function(ebdata,fsps,maxrl,maxar) {
  
  # declare vector with number of observations at each stage of cleaning
  obsv <- rep(NA,8)
  
  # columns to keep
  colkeep <- c("GLOBAL.UNIQUE.IDENTIFIER",
               "SCIENTIFIC.NAME",
               "EXOTIC.CODE",
               "COUNTRY",
               "STATE",
               "COUNTY",
               "LATITUDE",
               "LONGITUDE",
               "OBSERVATION.DATE",
               "TIME.OBSERVATIONS.STARTED",
               "SAMPLING.EVENT.IDENTIFIER",
               "PROTOCOL.NAME",
               "DURATION.MINUTES",
               "EFFORT.DISTANCE.KM",
               "EFFORT.AREA.HA",
               "NUMBER.OBSERVERS",
               "ALL.SPECIES.REPORTED",
               "GROUP.IDENTIFIER",
               "HAS.MEDIA")
               
  #ebdata <- ebdata[,c(1,7,10,16,18,29:32,35,37,41:48)]
  ebdata <- ebdata[,which(colnames(ebdata) %in% colkeep)] # use column names instead
  
  ## Count number of focal species records prior to data filtering
  # Define focal species
  # fsps <- c("Pharomachrus mocinno")
  #obsv[1] <- sum(ebdata$SCIENTIFIC.NAME==fsps) # initial number of fsps obs # <--- problem with NAs
  obsv[1] <- sum(ebdata$SCIENTIFIC.NAME==fsps & !is.na(ebdata$SCIENTIFIC.NAME))
  
  ## Get rid of records without explicit approval
  if(any(is.na(ebdata$APPROVED) | ebdata$APPROVED==0)) {
    ebdata <- ebdata[-which(is.na(ebdata$APPROVED) | ebdata$APPROVED==0),]
  }
  #obsv[2] <- sum(ebdata$SCIENTIFIC.NAME==fsps) # count obs after unnaproved out # <--- problem with NAs
  obsv[2] <- sum(ebdata$SCIENTIFIC.NAME==fsps & !is.na(ebdata$SCIENTIFIC.NAME))
  
  # get rid of records of escapee birds (exotic.code="X")
  if (any(is.na(ebdata$EXOTIC.CODE) | ebdata$EXOTIC.CODE == "X")) {
    ebdata <- ebdata[-c(which(is.na(ebdata$EXOTIC.CODE) | ebdata$EXOTIC.CODE == "X")),] 
  }

  # another option:
  #ebdata <- ebdata[ebdata$EXOTIC.CODE != "X" | is.na(ebdata$EXOTIC.CODE),]
  
  #obsv[3] <- sum(ebdata$SCIENTIFIC.NAME==fsps) # count obs after escapees out # <--- problem with NAs
  obsv[3] <- sum(ebdata$SCIENTIFIC.NAME==fsps & !is.na(ebdata$SCIENTIFIC.NAME))

  # get rid of all records that are not from complete lists
  ebdata <-ebdata[c(which(ebdata$ALL.SPECIES.REPORTED==1)),]
  obsv[4] <- sum(ebdata$SCIENTIFIC.NAME==fsps) # count obs after incomplete out
  
  # keep only observations from Area, Stationary, and Traveling protocols
  
  keep <- which(ebdata$PROTOCOL.NAME=="Area" | ebdata$PROTOCOL.NAME=="Stationary" | ebdata$PROTOCOL.NAME=="Traveling")
  ebdata <- ebdata[c(keep),]
  obsv[5] <- sum(ebdata$SCIENTIFIC.NAME==fsps) # count obs after protocol filter
  
  # find observations that have a group identifier. For each group keep event
  # with the biggest number of records
  groupids <- sort(unique(ebdata$GROUP.IDENTIFIER)) # get unique group ids
  groupids <- groupids[-which(groupids=="")]      # remove empty ids
  ng <- length(groupids)                          # get number of groups
  keep_event <- rep("NA",ng)                      # declare empty vector of events to keep
  ## get vector of events to keep and list with events per group
  epg <- lapply(split(ebdata$SAMPLING.EVENT.IDENTIFIER, ebdata$GROUP.IDENTIFIER)[-1], function(x) {
    table_x <- sort(table(x), decreasing = TRUE)
    names_x <- names(table_x)
  })
  keep_event <- sapply(epg, function(x) x[1])
  ## flatten list of events per group and subtract from there events to keep,
  ## in order to obtain events to discard (toss)
  groupevts <- unique(unlist(epg))              # get all events in groups
  toss <- groupevts[!groupevts%in%keep_event]   # gather events to toss in vector
  ebdata <- ebdata[which(!ebdata$SAMPLING.EVENT.IDENTIFIER%in%toss),] # toss
  ## Count number of bird records after data filtering
  obsv[6] <- sum(ebdata$SCIENTIFIC.NAME==fsps) # count obs after tossing mirrored lists
  
  ## identify events with more than maxrl EFFORT.DISTANCE.KM
  longevts <- ebdata$SAMPLING.EVENT.IDENTIFIER[which(ebdata$EFFORT.DISTANCE.KM>maxrl)]
  toss <- unique(longevts) # gather events to toss in vector
  ebdata <- ebdata[which(!ebdata$SAMPLING.EVENT.IDENTIFIER%in%toss),] # toss
  ## Count number of records after data filtering
  obsv[7] <- sum(ebdata$SCIENTIFIC.NAME==fsps) # count obs after tossing long routes
  
  ## identify events with more than maxar EFFORT.AREA.HA
  largevts <- ebdata$SAMPLING.EVENT.IDENTIFIER[which(ebdata$EFFORT.AREA.HA>maxar)]
  toss <- unique(largevts) # gather events to toss in vector
  ebdata <- ebdata[which(!ebdata$SAMPLING.EVENT.IDENTIFIER%in%toss),] # toss
  ## Count number of records after data filtering
  obsv[8] <- sum(ebdata$SCIENTIFIC.NAME==fsps) # count obs after tossing large areas
  
  
  # declare and fill output list
  out <- vector(mode="list",length=2)
  out[[1]] <- ebdata
  out[[2]] <- obsv
  
  return(out)
  
}
# end of function clean_ebird
# --------------------------------------
# clean_ebird()
