#' Filter Data for High Frequency Noise
#'
#' This functions filters the calculated depth of water using one or multiple consecutive Savitsky-Golay filters.
#'
#' @param data List of data divided into chunks.
#' @param response.col Index of the column indicating the depth at which the sensors measuring volumetric water content were installed.
#' @param nb.gs.filters Number of consecutive Savitsky-Golay filters applied to the data.
#' @return A data frame with the coefficients used to calculate the weighted average of volumetric water content, a data frame with calculate cumulative depth of water.
#' @export

FilterWCData <- function(data = NA, response.col = NA, nb.gs.filters = 2){ # begin function

  sensor.data.filter <- list() # create empty list to store filtered data
  for(element in 1:length(data)){ # begin iteration over elements in data

    select <- data[[element]] # select data in list
    select.filter <- list() # create list to store final results

    for(item in 1:length(select)) { # begin iteration over items in select

      subset <- select[[item]] # subset items in select
      subset$cumul.filter <- subset[,response.col] # create a new column to store filtered values, and set default value to original data

      for(nb in 1:nb.gs.filters){ # begin iteration over consecutive number of filters to be applied to the data
        vwc <- subset$cumul.filter # list variables to be filtered
        subset$cumul.filter <-  (-2*lag(vwc,3)+3*lag(vwc,2)+6*lag(vwc,1) +7*vwc + 6*lead(vwc,1) +3*lead(vwc,2)-2*lead(vwc,3))/21 # compute first filter
      } # end iteration over consecutive  number of filters to be applied to the data.

      select.filter <- list.append(select.filter, subset) # add data to sub-list

    } # end iteration over items in select

    sensor.data.filter <- list.append(sensor.data.filter, select.filter) # add data to final list

  } # end iteration over elements in data

  return(sensor.data.filter) # return final list

} # end function
