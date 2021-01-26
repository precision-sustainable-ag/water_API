#' Identify Infiltration Event within Available Water Data
#'
#' This function identifies infiltration events within calculated available water data.
#'
#' @param data List of data divided into chunks.
#' @param time.col Index of the timestamp column.
#' @param design.col List the induces of the optional columns characterizing the experimental design.
#' @param response.col Index of the column indicating the depth at which the sensors measuring volumetric water content were installed.
#' @param butter.n Butterworth filter order.
#' @param butter.w Critical frequency of the butterworth filter. Smaller w values increase the sensitivity of the created algorithm.
#' @param sensitivity Relative sensitivity of the algorithm. Lower values increase the sensitivity of the created algorithm.
#' @param approx.delay Estimated delay between the provided and filtered data. Default is 8 hours. Can be adjusted to improve algorithm performance.
#' @param temp.resolution Temporal resolution of the provided data, numeric value.
#' @param temp.unit Temporal resolution of the provided data, unit. Should be equal to "hour" or "min".
#' @return A data frame with interpolated missing volumetric water content.
#' @export

IdentifyInfiltrationEvents <- function(data = NA, time.col = NA, design.col = NA, response.col = NA,
                                       butter.n = 2, butter.w = 0.075, sensitivity = 0.1,
                                       approx.delay = 8, temp.resolution = 1, temp.unit = "hour"){ # begin function


  peaks.inventory <- c() # placeholder for final dataframe used to store data

  for(position in 1:length(data)){ # begin iteration over elements in data

    select <- data[[position]] # select element in list
    for(item in 1:length(select)) { # begin iteration over items in select

      subset <- select[[item]] # select item in select
      subset <- subset[!is.na(subset[,response.col]),] # remove observations with missing response variable values
      subset$timediff <- 1:NROW(subset) # calculate time difference between observations
      vwc <- subset[,response.col] # define response variable values

      # update data frame used to store results
      if(!is.data.frame(peaks.inventory)){ # if first time going through the loop, create dataframe used to store results
        peaks.inventory <- subset[0,design.col]; peaks.inventory[1,] <- NA # create columns according to data in subset
        peaks.inventory$position <- NA # add a position column, used to keep track of position in initial data list
        peaks.inventory$item <- NA # add an item column, used to keep track of position in initial data list
        peaks.inventory$time.min <- NA # add a time.min column, used to define absolute time at beginning of infiltration
        peaks.inventory$time.max <- NA # add a time.max column, used to define absolute time at end of infiltration
        peaks.inventory$depth.min <- NA # add a depth.min column, used to define depth of water at the beginning of infiltration
        peaks.inventory$depth.max <- NA # add a depth.max column, used to define depth of water at the end of infiltration
        peaks.inventory$timediff.min <- NA # add a timefiff.min column, used to define the relative beginning time of infiltration
        peaks.inventory$timediff.max <- NA # add a timediff.max column, used to define the relative end of infiltration
        peaks.inventory <- peaks.inventory[0,] # remove first row
      }

      # filter selected soil water content data with butterworth filter
      bf <- butter(butter.n, butter.w, type="low") # generate coefficients for lowpass butterworth filter
      vwc.filter <- signal::filter(bf, vwc) # filter the data


      # calculate local minimums and maximums using the second derivative of the filtered water content data
      local.mins <- which(diff(diff(vwc.filter)>0)==T)
      local.maxs <-  which(diff(diff(vwc.filter)<0)==T)

      ### automatically select local minimums and maximums showing changes in vwc caused by rainfall
      start <- 1
      if(length(na.omit(local.maxs[local.maxs<= local.mins[1]]))>0){
        start <- (local.maxs[local.maxs<= local.mins[1]])-4*approx.delay
        start <- max(start,1)
      }

      local.maxs <- local.maxs[local.maxs>local.mins[1]] #remove all maximums not preceded by a minimum

      if(length(na.omit(local.maxs))>0){ # begin if statement checking if there is at least one maximum preceded by a minimum

        #--------- Identify Peaks ---------

        # remove all the consecutive maxs for which there are no minimums in the middle
        both <- rbind(data.frame(value=local.mins, minormax = "min", stringsAsFactors = F),
                      data.frame(value=local.maxs, minormax = "max", stringsAsFactors = F)) # merge local min and max info into one data frame

        both <- both[order(both$value),] # order values within created data frame
        rownames(both) <- 1:NROW(both) # make sure rows within dataframe are named properly

        # initiate counters for while loop
        nrow.ini <- NROW(both)+1; nrow.end <- NROW(both)

        while(nrow.end<nrow.ini & NROW(both)>1) { #begin while loop established to select local mins and maxs showing changes due to rainfall

          nrow.ini <- NROW(both) # update counter value

          rows.delete <- c() # initiate vector used to determine which local min and max values should be removed
          rownames(both) <- 1:NROW(both) # make sure rows within dataframe are named properly
          max.list <- as.numeric(rownames(both[both$minormax == "max",])) # list all maxs within data frame
          min.list <- as.numeric(rownames(both[both$minormax == "min",])) # list all mins within data frame
          min.list <- split(min.list, cumsum(c(1, diff(min.list) != 1))) # list all consecutive mins within data frame
          max.list <- split(max.list, cumsum(c(1, diff(max.list) != 1))) # list all consecutive maxs within data frame

          for(element in 1:length(min.list)){ #  begin iteration over the list of consecutive mins

            rows.list <- unlist(min.list[[element]]) # select elements within the list
            if(length(rows.list)>1) { # if there is more than one element within the list...

              vwc.list <- vwc.filter[both$value[rows.list]] # ... select vwc of these elements
              rows.list <- rows.list[vwc.list != min(vwc.list)]  # ... select all minimums which does not minimize vwc
              rows.delete <- c(rows.delete, rows.list) # ... add them to vector for removal
            }}  # end if statement and iteration

          for(element in 1:length(max.list)){

            rows.list <- unlist(max.list[[element]]) # select elements within the list
            if(length(max.list[[element]])>1) { # if there is more than one element within the list ...

              vwc.list <- vwc.filter[both$value[rows.list]] # ... select vwc of these elements
              rows.list <- rows.list[vwc.list != max(vwc.list)]  # ... select all maximums which does not maximize vwc
              rows.delete <- c(rows.delete, rows.list) # ... add them to vector for removal
            }} # end if statement and iteration

          if(!is.null(rows.delete)){ both <- both[-rows.delete,]; rownames(both) <- 1:NROW(both)}  # remove selected mins and maxs

          if(NROW(both)>1){ # if at least one local maximum

            rows.delete <- c() # reset vector
            rownames(both) <- 1:NROW(both) # make sure rows within dataframe are named properly
            for (row in 2:(NROW(both))){ # begin iteration checking if there is a significant increase in vwc between mins and maxs
              if (both$minormax[row] == "max" & (vwc.filter[both$value[row]] - vwc.filter[both$value[row-1]]) < sensitivity) {
                rows.delete <- c(rows.delete, row) # if Delta vwc < 0.1 cm3/cm3 -> delete max value
              }}

            if(!is.null(rows.delete)){ both <- both[-rows.delete,]; rownames(both) <- 1:NROW(both)} # remove selected maxs

            nrow.end <- NROW(both) # update counter value

          } # end if statement checking if there is only one local maximum
        } # end while loop

        #--------- Complete peaks.inventory -------

        if(NROW(both)>1) { # if there is at least one local maximum

          for (row in seq(2, NROW(both), by=2)){ # begin iteration over maximum values in dataframe

            temp.inventory <- peaks.inventory[0,]
            temp.inventory[1,] <- NA

            for(col in 1:length(design.col)){temp.inventory[,col] <- unique(subset[,design.col[col]])}
            temp.inventory$position <- position
            temp.inventory$item <- item

            # extract vwc information from original (non-filtered) time series

            # vwc at local maximum
            temp.inventory$depth.max <- max(vwc[(max(start,both$value[row-1]-approx.delay)):both$value[row]])
            select.time <- subset[subset$timediff >= (start + both$value[row-1])/2 & subset$timediff <= both$value[row],]
            select.time <- na.omit(select.time[select.time$cumul.filter == temp.inventory$depth.max,]) # remove missing values

            if(NROW(select.time)>0) {temp.inventory$time.max <- select.time[NROW(select.time),time.col]} # identify time at local maximum

            if(!is.na(temp.inventory$time.max)){ # if time at local maximum was identified

              temp.inventory$timediff.max <- subset$timediff[subset[,time.col] == temp.inventory$time.max]  # determine relative time at local maximum

              # calculate local minimum preceding local maximum
              temp.inventory$depth.min <- min(vwc[(max(start,both$value[row-1]-approx.delay)):temp.inventory$timediff.max[1]])

              select.time <- subset[subset$timediff >= max(start,both$value[row-1]-approx.delay) & subset$timediff <= temp.inventory$timediff.max[1],]
              select.time <- na.omit(select.time[select.time$cumul.filter == temp.inventory$depth.min,]) # select timeframe where local minimum is expected
              if(NROW(select.time)>0) {temp.inventory$time.min  <- select.time[NROW(select.time),time.col]}


              if(!is.na(temp.inventory$time.min)){temp.inventory$timediff.min <- subset$timediff[subset[,time.col] == temp.inventory$time.min] } # determine relative time at local minimum


            } # end if statement checking if time at local maximum was identified

            temp.inventory <- within(temp.inventory, {time.min <- as.character(time.min); time.max <- as.character(time.max)})
            peaks.inventory <-  rbind(peaks.inventory, temp.inventory)  # add information to final dataset

          } # end iteration over maximum values
        } # end if statement checking if there is at least one local maximum



      } # end if statement checking if there is at least one  maximum preceded by a minimum
    } # end iteration over items in select
  } # end iteration over elements in data


  peaks.inventory <- peaks.inventory[(peaks.inventory$timediff.max - peaks.inventory$timediff.min) <= 4*24,] # define maximum time between a min and max to 3 days
  peaks.inventory$timediff.min <- peaks.inventory$timediff.min * temp.resolution
  peaks.inventory$timediff.max <- peaks.inventory$timediff.max * temp.resolution
  peaks.inventory$time.max[nchar(peaks.inventory$time.max) == 10] <- paste(peaks.inventory$time.max[nchar(peaks.inventory$time.max) == 10], '00:00:00')
  peaks.inventory$time.min[nchar(peaks.inventory$time.min) == 10] <- paste(peaks.inventory$time.min[nchar(peaks.inventory$time.min) == 10], '00:00:00')


  return(peaks.inventory)

} # end function
