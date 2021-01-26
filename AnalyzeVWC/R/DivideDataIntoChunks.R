#' Divide Cumulative Water Depth Data into Chunks of Continuous Data
#'
#' This function divides the dataframes of available water into chunks of continuous data
#'
#' @param data Data frame. Must contain at least 2 columns with one timestamp column in POSIXct/POSIXlt format and one column with the calculated depths of water. Additional columns can be added to characterize experimental design.
#' @param design.col List the induces of the optional columns characterizing the experimental design.
#' @param time.col Index of the timestamp column.
#' @param response.col Index of the column indicating the depth at which the sensors measuring volumetric water content were installed.
#' @param min.successive.obs Minimum number of observations in a chunk.
#' @param time.interval Temporal resolution.
#' @param time.unit Temporal unit.
#' @return A data frame with the coefficients used to calculate the weighted average of volumetric water content, a data frame with calculate cumulative depth of water.
#' @export

DivideIntoChunks <- function(data = NA, design.col = NA, time.col = NA, response.col = NA, min.successive.obs=12,
                             time.interval = 1, time.unit = "hour"){ # begin function


  # identify experimental units
  if(T %in% !is.na(design.col) & length(design.col)>1){ # if at least one column was specified in design
    for(var in length(design.col):1){data <- data[order(data[,design.col[var]]),]} # order data according to experimental design
    ncol <- NCOL(data) # save number of columns into variable
    data$design <- "design" # create a new column to identify experimental design
    for (var in design.col){data$design <- paste(data$design, data[,var],sep="#")} # update newly created column
  } # end if statement checking if at least one column was specified in design

  if(!(T %in% !is.na(design.col)) & length(design.col)==1){ # if no column were specified in design
    data$design <- "unit" # then all observations belong to the same experimental unit.
  }

  sensor.data.chunks <- list() # create list used to store final data
  for(unit in unique(data$design)) { # begin iteration of experimental units

    select <- data[data$design == unit,] # subset data
    select <- select[order(select[,time.col]),] # order data according to timestamp
    select <- select[!is.na(select[,response.col]),] # remove missing values
    subset <- list() # create temporary list used to store chunk of data associated with selected experimental unit

    if(NROW(select)>0){ # begin if statement checking if there is at least one observation in select

      row.names(select) <- 1:NROW(select) # rename rows

      select$timelag <- lag(as.character(select[,time.col])) # calculate time difference
      select$timediff <- NA; select$timediff[1] <- time.interval # create new column and set first value
      select$timediff[2:NROW(select)] <- difftime(as.POSIXlt(select[2:NROW(select),time.col], tz="GMT"), as.POSIXlt(select$timelag[2:NROW(select)], tz="GMT"), unit=time.unit) # calculate time difference
      gap <- as.numeric(c(1,row.names(select[select$timediff>time.interval,]),NROW(select))) #  list all rows corresponding to a time gap greater than one hour

      for(element in 1:(length(gap)-1)){ # begin iteration over rows corresponding to a time gap greater than one hour
        if(NROW(select[gap[element]:(gap[element+1]-1),1:(NCOL(select)-3)])>min.successive.obs){ # if there is at least the minimum observations in select
          subset <- list.append(subset,select[gap[element]:(gap[element+1]-1),1:(NCOL(select)-3)]) # add data to temporary list
        }}

    } # end if statement checking if there is at least one observation in select

    if(length(subset)>0){sensor.data.chunks <- list.append(sensor.data.chunks, subset)} # add data to final list

  } # end iteration over experimental units

  return(sensor.data.chunks) # return final list

} # end function
