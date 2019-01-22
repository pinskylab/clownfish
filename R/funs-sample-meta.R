# functions to help while maintaining the database
# one source for all of the helper functions to work with Leyte database and clownfish data


#' read_db is a helper function to connect to db with read only intent
#' @export
#' @importFrom dplyr src_mysql
#' @name read_db
#' @author Michelle Stuart
#'
#' @param x = which db?
#'
#' @examples
#' db <- read_Db("Leyte")

read_db <- function(x){
  library(dplyr)
  db <- src_mysql(dbname = x, default.file = path.expand("~/myconfig.cnf"), port = 3306, create = F, host = NULL, username = NULL, password = NULL)
  return(db)
}


#' write_db establishes a connection with the mysql db with intent to change it
#'
#' @param db_name
#'
#' @return database connection
#' @export
#'
#' @examples
#' leyte <- write_db("Leyte")

write_db <- function(db_name){
  db <- dbConnect(MySQL(), dbname = db_name, default.file = path.expand("~/myconfig.cnf"), port = 3306, create = F, host = NULL, user = NULL, password = NULL)
  return(db)
}

#' site_recap finds recaptured fish by site
#'
#' @param site_name
#'
#' @return a table of fish that have been recaptured at a given site
#' @export
#'
#' @examples
#' recaps <- site_recap("Haina")

site_recap <- function(site_name){

  fish <- fish_anem_dive() %>%
  # select all dives at a given site
    filter(site == site_name) %>%
    select(anem_table_id, anem_obs, anem_id, old_anem_id, sample_id,
           gen_id, anem_table_id, recap, tag_id) %>%
  # select all fish that have a tag_id or capid (because it is a recap, or are genetically recap)
    filter(recap == "Y" | !is.na(gen_id))

  return(fish)
}


#' for a given date range, list all of the dives
#' @export
#' @import dplyr
#' @name daterange_dive
#' @author Michelle Stuart
#' @param begin_date the beginning of the range
#' @param end_date the end of the range
#' @examples
#' dat <- daterange_dive("2016-01-01", "2017-12-30")

daterange_dive <- function(begin_date, end_date){
    dive <- get_dive() %>%
    filter(date > begin_date & date < end_date)

  return(dive)
}


#' sample_latlon
#'
#' @param sample_ids 
#'
#' @return a table of sample_ids with lat lons
#' @export
#'
#' @examples
#' where <- sample_latlon(c("APCL12_043", "APCL17_267"))
#' 
sample_latlon <- function(sample_ids){
  # find the anem_table_id for the sample
    fish <- fish_anem_dive() %>%
    select(sample_id,anem_table_id, anem_obs_time, dive_table_id, date, gps) %>%
    filter(sample_id %in% sample_ids) %>%
    # identify time zone as Asia
    mutate(anem_obs_time = lubridate::force_tz(ymd_hms(str_c(date, anem_obs_time, sep = " ")), tzone = "Asia/Manila"),
           # convert to UTC
           anem_obs_time = lubridate::with_tz(anem_obs_time, tzone = "UTC"),
           gpx_date = date(.data$anem_obs_time),
           gpx_hour = hour(anem_obs_time),
           minute = minute(anem_obs_time))


  # find the lat lon
  gpx <- leyte %>%
    tbl("GPX") %>%
    select(lat, lon, time, unit) %>%
    collect() %>%
    separate(time, into = c("gpx_date", "gps_time"), sep = " ") %>%
    mutate(gpx_date = date(.data$gpx_date)) %>%
    filter(gpx_date %in% fish$gpx_date) %>%
    separate(gps_time, into = c("gpx_hour", "minute", "second"), sep = ":") %>%
    filter(as.numeric(gpx_hour) %in% fish$gpx_hour & as.numeric(minute) %in% fish$minute) %>%
    mutate(gpx_hour = as.numeric(gpx_hour),
           minute = as.numeric(minute))

  # find matches for times to assign lat long - there are more than one set of seconds (sec.y) that match
  fish <- left_join(fish, gpx, by = c("gps" = "unit",  "gpx_date","gpx_hour", "minute")) %>%
    mutate(lat = as.numeric(lat),
           lon = as.numeric(lon)) # need to make decimal 5 digits - why? because that is all the gps can hold

  # calculate a mean lat lon for each anem observation
  coord <- fish %>%
    group_by(sample_id) %>% # id should be referring to one row of the data
    summarise(lat = mean(lat, na.rm = TRUE),
              lon = mean(lon, na.rm = T))

  return(coord)

}



#' for a given anemone , returns site, date, divetype
#' @export
#' @name anem_dive
#' @author Michelle Stuart
#' @param anem_ids anem(s) for which to get data
#' @examples
#' info <- anem_dive(2183)

anem_dive <- function(anems_ids){

  anem <- get_anem() %>%
    filter(anem_id %in% anem_ids)

  dive <- get_dive() %>%
    filter(dive_table_id %in% anem$dive_table_id)

  anem <- left_join(anem, dive, by = "dive_table_id")

  return(anem)
}




#' once rows have been changed remove them from table and add them back in
#' @export
#' @name change_rows
#' @author Michelle Stuart
#' @param x = whole table
#' @param y = the changed rows
#' @param z = identifying co
#' @examples
#' deer <- change_rows(deer, change)

change_rows <- function(table, change, identifier){
  table <- anti_join(table, change, by = identifier)
  table <- rbind(table, change)
  return(table)
}


#' Get anemone data
#'
#' @return the anemone db table
#' @export
#' @import dplyr
#'
#' @examples
#' anem_17 <- get_anem() %>% filter(anem_id > 2000)

get_anem <- function(){
  leyte <- read_db("Leyte")
  anem <- leyte %>%
    tbl("anemones") %>%
    collect()

  return(anem)
}

#' Get dive data
#'
#' @return the diveinfo db table
#' @export
#' @import dplyr
#'
#' @examples
#' dive <- get_dive()
get_dive <- function(){
  leyte <- read_db("Leyte")
  dive <- leyte %>%
    tbl("diveinfo") %>%
    collect()

  return(dive)
}

#' Get Clownfish data
#'
#' @return the clownfish db table
#' @export
#' @import dplyr
#'
#' @examples
#' fish <- get_fish()
get_fish <- function(){
  leyte <- read_db("Leyte")
  fish <- leyte %>%
    tbl("clownfish") %>%
    collect()

  return(fish)
}

#' fish_anem_dive grabs this data for fish from the db
#'
#'
#' @return a table of all fish, associated anemones, and dive info
#' @export
#'
#' @examples
#' fish_meta <- fish_anem_dive()

fish_anem_dive <- function(){
  fish <- get_fish()

  anem <- get_anem() %>%
    filter(anem_table_id %in% fish$anem_table_id)

  fish <- left_join(fish, anem, by = "anem_table_id")

  dive <- get_dive() %>%
    filter(dive_table_id %in% fish$dive_table_id)

  fish <- left_join(fish, dive, by = "dive_table_id")

  return(fish)

}
