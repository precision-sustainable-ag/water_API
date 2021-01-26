#' Interpolate Missing Volumetric Water Content Data
#'
#' This function interpolates missing volumetric water content data according to the user's specification. Missing volumetric water content values are estimated using a linear interpolation.
#'
#' @param data Data frame. Must contain at least 3 columns with one timestamp column in POSIXct/POSIXlt format, one column to indicate the depth at which the sensors measuring volumetric water content were installed, and one column with the measured volumetric water content data in the numeric format. Additional columns can be added to characterize experimental design.
#' @param temp.resolution Temporal resolution of the provided data, numeric value.
#' @param temp.unit Temporal resolution of the provided data, unit. Should be equal to "hour" or "min".
#' @param time.col Index of the timestamp column.
#' @param response.col Index of the column indicating the depth at which the sensors measuring volumetric water content were installed.
#' @param design.col List the induces of the optional columns characterizing the experimental design.
#' @param depth.col Index of the column containing the provided measured volumetric water content data.
#' @param depths.values List the depth values at which the sensors measuring volumetric water content were installed.
#' @param missing.max.time.by.depth List the maximum amount of missing time allowed for interpolation. Position of the element in the list corresponds to the corresponding depth in depths.values. Values are provided in number of temp.resolution intervals. Set value to NA if data should be interpolated no matter the number of missing observations.
#' @return A data frame with interpolated missing volumetric water content.
#' @export

InterpolateMissingVWCData <- function(data = NA, temp.resolution = 1, temp.unit = "hours", time.col = NA, response.col = NA,
                                   design.col = NA, depth.col = NA, depths.values = c(-15, -45, -80), missing.max.time.by.depth = c(12,NA,NA)){ # begin function

  # add depth column to design.col
  design.col <- c(design.col, depth.col)

  # identify experimental units
  if(T %in% !is.na(design.col) & length(design.col)>1){ # if at least one column was specified in design.col
    for(var in length(design.col):1){ data <- data[order(data[,design.col[var]]),]} # order data according to experimental design.col
    ncol <- NCOL(data) # save number of columns into variable
    data$design.col <- ""  # create a new column to identify experimental design.col
    for (var in design.col){data$design.col <- paste(data$design.col, data[,var],sep="")} # update newly created column
  }

  if(!(T %in% !is.na(design.col)) & length(design.col)==1){ # if no column were specified in design.col
    data$design.col <- "unit" # then all observations belong to the same experimental unit.
  }

  # begin iteration over design.col parameters
  data.all <- data[0,] # create data frame used to store results

  for(unit in unique(data$design.col)) { # begin iteration of experimental units

    select <- data[data$design.col == unit,] # subset data for selected experimental unit
    select <- select[!(duplicated(select[,time.col])),] # remove duplicates

    timestamp.min <- min(select[,time.col]) # identify minimum timestamp
    timestamp.max <- max(select[,time.col])  # identify maximum timestamp
    timestamps <- as.character(seq(timestamp.min, timestamp.max, by = paste(temp.resolution, temp.unit)))  # list all possible timestamps

    select[,time.col] <- as.character(select[,time.col])# create final dataset

    for (time in timestamps[!(timestamps %in% select[,time.col])]){ #  begin iteration in timestamps not in data

      temp <- select[1,]; temp[,1:NCOL(temp)] <- NA # create new observation
      temp[1,time.col] <- time # add timestamps

      # add constant values
      if(T %in% !is.na(design.col) & length(design.col)>1){ # fill in the blanks for constant values, if applicable
        for(col in design.col){
          temp[1,col] <- unique(select[,col]) }}

      select <- rbind(select, temp) # add output to subset
    } # end iteration over timestamps not in data

    # interpolate missing data
    select <- select[order(select[,time.col]),] # order data by timestamp
    ncol <- NCOL(select) # define number of column

    for(col in response.col){
      select[NCOL(select)+1] <- select[,col] # add column to store interpolated values
    }
    colnames(select)[(ncol+1):(ncol+length(response.col))] <- paste(colnames(select)[response.col],".interpolated", sep="") # rename new column

    for(col in 1:length(response.col)){
      select[,ncol+col] <- na_interpolation(select[,ncol+col], option = "linear")
    }

    # remove interpolated surface data if applicable
    if(!is.na(missing.max.time.by.depth[depths.values == unique(select[,depth.col])])){

      for(col in 1:length(response.col)){

        missing.list <- c()

        for(row in 2:(NROW(select)-1)){ #begin iteration over rows

          if (is.na(select[row,response.col[col]])) {missing.list <- c(missing.list, row)}
          if (length(missing.list)>missing.max.time.by.depth[depths.values == unique(select[,depth.col])] & !is.na(select[row+1,response.col[col]])){select[missing.list,ncol+col] <- NA}
          if (!is.na(select[row+1,response.col[col]])) {missing.list <- c()}
        } # end iteration over rows
      }
    }

    # reorganize data
    colnames <- colnames(select)[response.col]
    select <- select[,-response.col]
    colnames(select)[(ncol-length(response.col)+1):ncol] <- colnames

    data.all <- rbind(data.all, select) # add output to final dataset
  } # end iterations over experimental units

  data.all <- data.all[,-(ncol-length(response.col))]

  for(interpolated.col in 1:length(response.col)){
    data.all[,NCOL(data.all)-interpolated.col + 1] <- round(data.all[,NCOL(data.all)-interpolated.col + 1], digits=1)
  }


  # return final data frame
  return(data.all)


} # end function



