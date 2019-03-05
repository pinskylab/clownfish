# lab helpers - helper functions for lab work

#' plate_from_db recreates the platemap based on well locations in the db
#'
#' @param table_name a table imported from the database that includes id and well columns
#' @param ... id to be displayed on the plate map
#'
#' @return a platemap representing where samples went in a plate
#' @export
#' @import dplyr
#' @importFrom reshape2 acast
#'
#' @examples
#' platemap <- plate_from_db(extractions, sample_id)
plate_from_db <- function(table_name, ...){
  # split the well out into row and column 
  table_name$row <- substr(table_name$well, 1, 1)
  table_name$col <- as.numeric(substr(table_name$well, 2, 3))
  
  # select columns for plate 
    table_name <- table_name %>% 
      select(row, col, ...) %>% #keep row & col, identifier
      arrange(row, col)

  table_name <- as.data.frame(table_name)
  
  # make map
  platemap <<- as.matrix(reshape2::acast(table_name,table_name[,1] ~ table_name[,2]))
  return(table_name)
}


#' Title
#'
#' @param list_of_ids a list of ids
#' @param id_type character string of id column name
#'
#' @return a table of ids and well locations
#' @export
#'
#' @examples
#' plate <- make_plate_with_negs(samples, "extraction_id")


make_plate_with_negs <- function(list_of_ids, id_type){
  # make a dataframe of the list_of_ids
  ids <- data.frame(list_of_ids, stringsAsFactors = F)
 
  # how many rows are in the table (how many samples)?
  y <- nrow(ids)

  # how many plates would these make, 94 samples plus 2 blanks per plate
  (nplates <- floor(y/94)) # extra parenthesis are to print

  # define wells
  well <- 1:(96*nplates)

  # insert the negative controls and set up the plate
  plate <- data.frame() # blank data frame to build upon
  for (i in 1:nplates){
    c <- 96*i-95 # well 1 on a plate
    d <- 96*i-85 # 11
    e <- 96*i-84 # 12 negative control well
    f <- 96*i-83 # 13
    g <- 96*i-36 # 60
    h <- 96*i-35 # 61 negative control well
    j <- 96*i-34 # 62
    k <- 96*i-2  # 94
    l <- 96*i - 37 # 59
    m <- 96*i #96
    str1 <- as.data.frame(cbind(well[c:d], ids[c:d,])) # 1:11
    names(str1) <- c("well", "id_type")
    str2 <- as.data.frame(cbind(well[e], "XXXX")) # because the first blank is in the 12th position
    names(str2) <- c("well", "id_type")
    str3 <- as.data.frame(cbind(well[f:g], ids[e:l,])) #13:60 in plate, 12:59 in list
    names(str3) <- c("well", "id_type")
    str4 <- as.data.frame(cbind(well[h], "XXXX")) # because the 2nd blank is in the 61st position
    names(str4) <- c("well", "id_type")
    str5 <- as.data.frame(cbind(well[j:k], ids[g:k,]))# 62:96 in plate, 60:94 in list
    names(str5) <- c("well", "id_type")
    
    # and stick all of the rows together
    temp <- data.frame(rbind(str1, str2, str3, str4, str5))
    temp$row <- rep(LETTERS[1:8], 12)
    temp$col <- unlist(lapply(1:12, rep, 8))
    temp$plate <- paste("plate", i, sep = "")
    plate <- rbind(plate, temp)
    
  }
  
  # put the samples in order of id (with negative controls inserted)
  plate <- arrange(plate, plate, col, row)
  
  return(plate)
}




#' Make Platemap
#'
#' @param tble a table of samples with well locations
#'
#' @return a platemap
#' @export
#' @import dplyr
#' @importFrom reshape2 acast
#'
#' @examples
#' platmap <- make_platemap(plate)


make_platemap <- function(tble){
  plate <- tble %>% 
    select(contains("id"), well) %>% 
    mutate(row = substr(well, 1, 1), 
      col = as.numeric(substr(well, 2, 3)))
  plate <- plate %>% 
    select(row, col, contains("id"))
  platemap <<- as.matrix(reshape2::acast(plate, plate$row ~ plate$col))
}





#' Remove rows
#'
#' @param table_name table of samples
#' @param how_many_wells how many wells in the plate
#'
#' @return a sliced table
#' @export
#' @import dplyr
#'
#' @examples
#' samp <- remove_rows(samp, 48)
#' 
remove_rows <- function(table_name, how_many_wells){
  x <- nrow(table_name) %% how_many_wells # get the remainder after dividing by 48
  table_name <- table_name %>% 
    select(1) %>% 
    arrange() %>% 
    slice(1:(nrow(table_name) - x))
  
  return(table_name)
  
}


#' lig_ng figure out how many ng to use in making pools for ligations
#'
#' @param dig a table of digests that need to be ligated
#' @param regeno a table of samples to be regenotyped
#' @return a table of ligations with volumes of sample and water to combine
#' @export
#' @import dplyr
#' @importFrom tibble tibble
#' @name lig_ng
#' @author Michelle Stuart
#' @examples 
#' ligs <- lig_ng(dig)
lig_ng <- function(dig, regeno) {
  out <- tibble() # make a blank data frame to write to
  y <- nrow(dig) # get the number of beginning rows
  for(i in c(50, 75, 100, 150, 200)){
    if (nrow(out) < y){ # if we haven't placed all of the candidates yet
      # define the samples that can be ligated at the current ng
      ng <- dig %>%
        mutate(uL_in = round(i/quant, 1)) %>% # round to 1 decimal point
        filter(uL_in < 22.2 & uL_in > 0.5) %>%
        mutate(water = round(22.2-uL_in, 1), 
          DNA = i)
      # define regenos that can be licated at the current ng
      reg <- regeno %>%
        mutate(uL_in = round(i/quant, 1)) %>% # round to 1 decimal point
        filter(uL_in < 22.2 & uL_in > 0.5) %>%
        mutate(water = round(22.2-uL_in, 1), 
          DNA = i)
      # pull out  pools
      while (nrow(ng) >= 47){
        while(nrow(reg) >= 1){
        pool <- ng %>% 
          slice(1:47)
        re <- reg %>%
          slice(1)
        ng <- anti_join(ng, pool, by ="digest_id")
        reg <- anti_join(reg, re, by = "digest_id")
        out <- rbind(out, pool, re)
        dig <- anti_join(dig, ng, by = "digest_id")
        regeno <- anti_join(regeno, ng, by = "digest_id")
      }
      }
    }
  }
  return(out)
}


#' make a plate from a list of sample_ids, extraction_ids, etc.
#'
#' @param list_of_ids a list of ids
#'
#' @export
#' @return a table of plate locations for samples
#' @name make_plate
#' @author Michelle Stuart
#' @examples 
#' test <- make_plate(lig_ids)

make_plate <- function(list_of_ids){
  # make a dataframe of the list_of_ids
  ids <- as.data.frame(list_of_ids)
  
  # how many rows are in the table (how many samples)?
  y <- nrow(ids)
  
  if (y >= 96){
    
  # how many plates would these make
  (nplates <- floor(y/96)) # extra parenthesis are to print
  
  # remove those rows that don't fit into plates
  ids <- remove_rows(ids, 96)
  
  # define wells
  well <- 1:(96*nplates)
  
  # set up the plate
  plate <- data_frame()
  for (i in 1:nplates){
    a <- 96*i-95 # position 1
    b <- 96*i     # position 96
    temp <- cbind(well[a:b], as.data.frame(ids[a:b, ]))
    temp$row <- rep(LETTERS[1:8], 12)
    temp$col = unlist(lapply(1:12, rep, 8))
    temp$plate = paste("plate", i, sep = "")
    plate <- rbind(plate, temp)
  }
  
  # put plate in order
  plate <- arrange(plate, plate, col, row)
  
  }else{
    plate <- data.frame( Row = rep(LETTERS[1:8], 12), Col = unlist(lapply(1:12, rep, 8)))
    plate <- plate[1:y,]
    plate <- cbind(plate, ids)
    plate$plate <- "shortplate1"
    }
    
  return(plate)
}




#' assign a location on the robot table for a destination or source plate
#'
#' @param plate_names names of plates to be used on the robot
#' @param table table of samples
#' @param dest_or_source is this plate a source or a destination
#' @param identifier digest_id or ligation_id
#'
#' @export
#' @name assign_mek_loc
#' @author Michelle Stuart
#' @examples 
#' source <- assign_mek_loc(dig_plates, source, "source", "digest_id")

assign_mek_loc <- function(plate_names, table, dest_or_source, identifier){
  for (i in 1:nrow(plate_names)){
    if (dest_or_source == "dest"){
      change <- table %>% 
        filter(plate == plate_names$plate[i]) %>% 
        mutate(dest = mek_loc[length(mek_loc)])
      mek_loc <<- mek_loc[1:length(mek_loc)-1]
      table <- change_rows(table, change, identifier)
    }else{
      change <- table %>% 
        filter(plate == plate_names$plate[i]) %>% 
        mutate(source = mek_loc[length(mek_loc)])
      mek_loc <<- mek_loc[1:length(mek_loc)-1]
      table <- change_rows(table, change, identifier)
    }
    
  }
  return(table)
  
  
}


#' find sample id from ligation id
#'
#' @param table_name table containing ligation ids
#'
#' @export
#' @import dplyr
#' @name samp_from_lig
#' @author Michelle Stuart
#' @examples 
#' c5 <- samp_from_lig(genedf)


samp_from_lig <- function(table_name){
  
  if(!exists(lab))
    stop("Error: db connection called 'lab' does not exist, see Michelle for help")
  
  # connect ligation ids to digest ids
  lig <- get_lig() %>% 
    filter(ligation_id %in% table_name$ligation_id) %>% 
    select(ligation_id, digest_id)
  
  # connect digest ids to extraction ids
  dig <- get_dig() %>% 
    filter(digest_id %in% lig$digest_id) %>% 
    select(extraction_id, digest_id)
  
  extr <- get_extr() %>% 
    filter(extraction_id %in% dig$extraction_id) %>% 
    select(extraction_id, sample_id)
  
  mid <- left_join(lig, dig, by = "digest_id")
  samp <- left_join(extr, mid, by = "extraction_id") %>% 
    select(sample_id, ligation_id)
  
  return(samp)
}


#' heatmap - plot a plate map with color
#'
#' @param plate_as_long_table a table of samples and plate locations
#' @param id sample identifier
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @name heatmap
#' @author Michelle Stuart
#' @examples 
#' sample_map <- heatmap(plate_as_long_table)

heatmap <- function(plate_as_long_table, id){
  map <- plate_as_long_table  %>% 
    mutate(row = substr(well, 1, 1),
      row = factor(row, levels = c("H", "G", "F", "E", "D", "C", "B", "A")), 
      col = substr(well, 2, 3),
           col = factor(col, levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))) %>% 
    select(row, col, contains("id"), filter)
  
  plateheatmap <- ggplot(map, aes(x=col, y=row, fill= filter)) + 
    geom_tile()
  
  z <- plateheatmap + 
    geom_text(aes(col, row, label = id), color = "black", size = 4) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank())
  return(z)
  
}


 
# check the work history of those sample_ids
#' Title
#'
#' @param table table_where_ids_are
#' @param id_column column_of_ids - must be sample_id, extraction_id, digest_id, or ligation_id
#' @return a table of lab metadata
#' @export
#'
#' @examples
#' 
work_history <- function(table, id_column){
  if(!exists(lab))
    stop("Error: db connection called 'lab' does not exist, see Michelle for help")
  
  if(id_column == "sample_id"){
    hist <- get_extr() %>% 
      filter(sample_id %in% table$sample_id) %>% 
      select(sample_id, extraction_id, well, plate) %>% 
      rename(extr_well = well, 
             extr_plate == plate)
    
    dig <- get_dig() %>% 
      filter(extraction_id %in% hist$extraction_id) %>% 
      select(extraction_id, digest_id, well, plate) %>% 
      rename(dig_well = well, 
             dig_plate = plate)
    
    hist <- left_join(hist, dig, by = "extraction_id")
    
    
    lig <- get_lig() %>% 
      filter(digest_id %in% hist$digest_id) %>% 
      select(ligation_id, digest_id, well, plate, barcode_num, pool) %>% 
      rename(lig_well = well, 
             lig_plate = plate)
    hist <- left_join(hist, lig, by = "digest_id")
  
    return(hist)
  }  
  
  if(id_column == "extraction_id"){
    hist <- get_extr() %>% 
      filter(extraction_id %in% table$extraction_id) %>% 
      select(sample_id, extraction_id, well, plate) %>% 
      rename(extr_well = well, 
             extr_plate == plate)
    
    dig <- get_dig() %>% 
      filter(extraction_id %in% hist$extraction_id) %>% 
      select(extraction_id, digest_id, well, plate) %>% 
      rename(dig_well = well, 
             dig_plate = plate)
    
    hist <- left_join(hist, dig, by = "extraction_id")
    
    
    lig <- get_lig() %>% 
      filter(digest_id %in% hist$digest_id) %>% 
      select(ligation_id, digest_id, well, plate, barcode_num, pool) %>% 
      rename(lig_well = well, 
             lig_plate = plate)
    
    hist <- left_join(hist, lig, by = "digest_id")
    
    return(hist)
  }  
  
  if(id_column == "digest_id"){
    dig <- get_dig() %>% 
      filter(digest_id %in% table$digest_id) %>% 
      select(extraction_id, digest_id, well, plate) %>% 
      rename(dig_well = well, 
             dig_plate = plate)
    
    hist <- get_extr() %>% 
      filter(extraction_id %in% dig$extraction_id) %>% 
      select(sample_id, extraction_id, well, plate) %>% 
      rename(extr_well = well, 
             extr_plate == plate)
    
    
    hist <- left_join(hist, dig, by = "extraction_id")
    
    
    lig <- get_lig() %>% 
      filter(digest_id %in% hist$digest_id) %>% 
      select(ligation_id, digest_id, well, plate, barcode_num, pool) %>% 
      rename(lig_well = well, 
             lig_plate = plate)
    
    hist <- left_join(hist, lig, by = "digest_id")
    
    return(hist)
  }  
  
  if(id_column == "ligation_id"){
    lig <- get_lig() %>% 
      filter(ligation_id %in% table$ligation_id) %>% 
      select(ligation_id, digest_id, well, plate, barcode_num, pool) %>% 
      rename(lig_well = well, 
             lig_plate = plate)
    
    dig <- get_dig() %>% 
      filter(digest_id %in% lig$digest_id) %>% 
      select(extraction_id, digest_id, well, plate) %>% 
      rename(dig_well = well, 
             dig_plate = plate)
    
    hist <- left_join(dig, lig, by = "digest_id")
    
    extr <- get_extr() %>% 
      filter(extraction_id %in% dig$extraction_id) %>% 
      select(sample_id, extraction_id, well, plate) %>% 
      rename(extr_well = well, 
             extr_plate == plate)
    
    hist <- left_join(hist, extr, by = "extraction_id")
  
    return(hist)
  }  
}


#' Get Lig
#'
#' @param ligation_ids a table that contains ligation_ids
#'
#' @return a table of ligation ids with meta data
#' @export
#'
#' @examples
#' test <- get_lig(table$ligation_id)


get_lig <- function(){
  
  if(!exists(lab))
    stop("Error: db connection called 'lab' does not exist, see Michelle for help")
  
  lig <- lab %>% 
    tbl("ligation") %>% 
    collect() 
  
  return(lig)
}


#' Get Dig
#'
#' @return a table of digests with meta data
#' @export
#'
#' @examples
#' test <- get_dig() %>% select(digest_id, extraction_id)
#' 
get_dig <- function(){
  if(!exists(lab))
    stop("Error: db connection called 'lab' does not exist, see Michelle for help")
  
  dig <- lab %>% 
    tbl("digest") %>% 
    collect()
  
  return(dig)
}

#' Get Extr
#'
#' @return a table of extractions with meta data
#' @export
#'
#' @examples
#' test <- get_extr() %>%  select(sample_id)


get_extr <- function(){
  if(!exists(lab))
    stop("Error: db connection called 'lab' does not exist, see Michelle for help")
  
  extr <- lab %>% 
    tbl("extraction") %>% 
    collect()
  
  return(extr)
}
