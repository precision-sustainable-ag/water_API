#' Calculate Available Water from Individual Sensor Measurements
#'
#' This function calculate cumulative water depth from individual volumetric soil water content measurements
#'
#' @param data Data frame. Must contain at least 3 columns with one timestamp column in POSIXct/POSIXlt format, one column to indicate the depth at which the sensors measuring volumetric water content were installed, and one column with the measured volumetric water content data in the numeric format. Additional columns can be added to characterize experimental design.
#' @param design.col List the induces of the optional columns characterizing the experimental design.
#' @param time.col Index of the timestamp column.
#' @param depth.col Index of the column containing the provided measured volumetric water content data.
#' @param response.col Index of the column indicating the depth at which the sensors measuring volumetric water content were installed.
#' @param profile.min Lower depth of the soil profile.
#' @param profile.max Upper depth of the soil profile.
#' @param increment Depth of the soil profile layers used to calculate cumulative depths of water in the profile.
#' @param sensor.height Thickness of soil affecting the volumetric soil water content measurement.
#' @param divideby100 Equals T if the volumetric water content values are provided in percent rather than decimal values.
#' @return A data frame with the coefficients used to calculate the weighted average of volumetric water content, a data frame with calculate cumulative depth of water.
#' @export

CalculateCumulativeWaterDepth <- function(data = NA, design.col = NA, time.col = NA, depth.col = NA,
                                          response.col = NA, profile.min = NA, profile.max = NA, increment = 5,
                                          sensor.height = 10, divideby100 = T){ # begin function

  # readjust volumetric water content values if needed
  if (divideby100){data[,response.col] <- data[,response.col]/100}

  # define soil layers according to profile.min, profile.max, and sensor.heigh
  buildformula <- data.frame(layer.low = seq(profile.min, profile.max - increment, by = increment),
                             layer.high = seq(profile.min + increment, profile.max, by=increment))

    # add columns to define formula
    for(col in 1:length(unique(data[,depth.col]))){
      buildformula[,NCOL(buildformula)+1] <- 0 # create column and set values to zero
      colnames(buildformula)[NCOL(buildformula)] <- paste0("D",col) # rename new column
    }


    # define coefficients when sensor is in soil layer
    depths <- unique(data[,depth.col])
    depths <- depths[order(depths)] # order values

    for(depth in 1:length(depths)){ # edit buildformula data frame
      buildformula[between(buildformula$layer.low, depths[depth]-sensor.height/2+sensor.height/100, depths[depth]+sensor.height/2-sensor.height/100) |
                     between(buildformula$layer.high, depths[depth]-sensor.height/2+sensor.height/100, depths[depth]+sensor.height/2-sensor.height/100) ,2+depth] <- 1
    }

      # define coefficients when sensor is not in soil layer
      for(depth in 1:(length(depths)-1)){ # begin iteration by depth

        nb.layers <- (depths[depth+1] - depths[depth] - sensor.height) / increment # calculate the number of layers between two sensors

        if(nb.layers>0){ # begin if statement checking if there is at least one layer without sensors
          select <- buildformula[buildformula$layer.low < (depths[depth+1]-sensor.height/2) &
                                   buildformula$layer.high > (depths[depth]-sensor.height/2) &
                                   buildformula[,2+depth] == 0,] # select layers between sensors

          counter1 <- nb.layers # define counter to keep track of coefficient for lower depth
          counter2 <- 1 # define counter to keep track of coefficient for upper  depth
          for(row in 1:NROW(select)){ # begin iteration over layers between sensors
            buildformula[buildformula$layer.low == select$layer.low[row], 2+depth] <- counter1 # define coefficient for lower depth
            buildformula[buildformula$layer.low == select$layer.low[row], 2+depth+1] <- counter2 # define coefficient for upper depths
            counter1 <- counter1 - 1 # update counter for lower depth
            counter2 <- counter2 + 1 # update counter for upper depth
          } # end iteration over layers between sensors
        } # end if statement checking if there is at least one layer between sensors
      } # end iteration by depth

      # calculate weight
      buildformula$weight <- buildformula$D1 # create a new column
      for(depth in 2:length(depths)){buildformula$weight <- buildformula$weight + buildformula[,2+depth]} # calculate weighted average coefficients

      # format data
      data[,time.col] <- as.character(data[,time.col])
      original.colnames <- colnames(data)
      colnames(data)[time.col] <- "time"
      colnames(data)[depth.col] <- "sensor.depth"
      colnames(data)[response.col] <- "vwc"

      # calculate cumulative depth of water
      data.cumul <- c()

      if(T %in% !is.na(design.col) & length(design.col)>1){ # if at least one column was specified in design.col
        for(var in length(design.col):1){ data <- data[order(data[,design.col[var]]),]} # order data according to experimental design.col
        ncol <- NCOL(data) # save number of columns into variable
        data$design <- "design"  # create a new column to identify experimental design.col
        for (var in design.col){data$design <- paste(data$design, data[,var],sep="#")} # update newly created column
      }

      if(!(T %in% !is.na(design.col)) & length(design.col)==1){ # if no column were specified in design.col
        data$design <- "unit" # then all observations belong to the same experimental unit.
      }

      for(unit in unique(data$design)) { # begin iteration of experimental units


         if(length(data.cumul) == 0){
          data.cumul <- as.data.frame(data[data$design==unit,] %>%
            group_by(time,sensor.depth) %>%
            summarize(mean(vwc)) %>%
            spread(sensor.depth, `mean(vwc)`)) %>%
            add_column(design=unit)
         } else {
           temp <- as.data.frame(data[data$design==unit,] %>%
                                   group_by(time,sensor.depth) %>%
                                   summarize(mean(vwc)) %>%
                                   spread(sensor.depth, `mean(vwc)`)) %>%
                                   add_column(design=unit)
           data.cumul <- rbind(data.cumul, temp)
         }
      } # end iteration over experimental units


    # re-define design columns
      for(col in 1:length(design.col)){
        data.cumul[,NCOL(data.cumul)+1] <- NA
        colnames(data.cumul)[NCOL(data.cumul)] <- original.colnames[design.col[col]]
        for(row in 1:NROW(data.cumul)){ data.cumul[row,NCOL(data.cumul)] <- strsplit(data.cumul$design[row],split="#")[[1]][col+1]}
      }

      # reorganize columns
      data.cumul <- data.cumul[,c((NCOL(data.cumul)-length(design.col)+1):NCOL(data.cumul),1:(NCOL(data.cumul)-length(design.col)-1))]

        # calculate cumulative depth of water
        data.cumul$cumul.depth <- 0 # create new column and define default value as 0

        for(row in 1:NROW(buildformula)){ # begin iteration over soil layer

          addition <- 0 # create counter used to calculate cumulative depth of water across the profile
          for(depth in 1:length(depths)){addition <- addition + buildformula[row,2+depth] * data.cumul[,length(design.col)+1+depth]} # update counter
          data.cumul$cumul.depth <- data.cumul$cumul.depth + addition/buildformula$weight[row] # update data frame provided by the user

        } # end iteration over soil layer

      data.cumul$cumul.depth <- round(data.cumul$cumul.depth / NROW(buildformula) * (profile.max - profile.min), digits = 3) # round to 3 digits
      data.cumul$time <- as.POSIXlt(data.cumul$time, tz='GMT')

        return(list(buildformula, data.cumul))# return dataframe with cumulated depth of water across the profile

  } # end function








