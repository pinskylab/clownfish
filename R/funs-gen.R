# read_genepop ####
#' read a genepop generated for ONE population
#' @export
#' @name read_genepop
#' @author Michelle Stuart
#' @param x = filename
#' @examples 
#' genedf <- read_genepop("data/seq17_03_58loci.gen")

# only works if Genepop has loci listed in 2nd line (separated by commas), has individual names separated from genotypes by a comma, and uses spaces between loci
# return a data frame: col 1 is individual name, col 2 is population, cols 3 to 2n+2 are pairs of columns with locus IDs
read_genepop <-  function(filename){
  # get all of the data
  if(readLines(filename, n=3)[[3]] != "pop"){
  dat <- read.table(filename, skip = 2, sep = " ", stringsAsFactors = F, colClasses = "character") 
  }else{
    dat <-  read.table(filename, skip = 3, sep = " ", stringsAsFactors = F, colClasses = "character")
  }
  
  # get the header info
  info <- readLines(filename, n = 2) 
  
  # define the loci names
  loci <- unlist(strsplit(info[2], split=','))    
  
  # rename the dat columns
  names(dat) <- c("names", loci)
  
  # if there is a comma, remove it
  dat$names <- str_replace(dat$names, ",", "")
  
  return(dat)	
  
}

# lig_from_samp ####
#' views all of the fish recaptured at a given site
#' @export
#' @name lig_from_samp
#' @author Michelle Stuart
#' @param x = list of sample_ids
#' @examples 
#' fish <- lig_from_samp(c("APCL13_516", "APCL13_517"))

lig_from_samp <- function(sample_ids){
  
  lab <- read_db("Laboratory")
  
  extr <- lab %>% 
    tbl("extraction") %>% 
    filter(sample_id %in% sample_ids) %>% 
    select(sample_id, extraction_id) %>% 
    collect()
  
  dig <- lab %>% 
    tbl("digest") %>% 
    filter(extraction_id %in% extr$extraction_id) %>%
    select(extraction_id, digest_id) %>% 
    collect()
  
  lig <- lab %>% 
    tbl("ligation") %>% 
    filter(digest_id %in% dig$digest_id) %>%
    select(ligation_id, digest_id) %>% 
    collect()
  
  mid <- left_join(extr, dig, by = "extraction_id")
  lig <- left_join(mid, lig, by = "digest_id") 
  
  return(lig)
}




# check_id_match ####
#' compare individuals that appear to be genetic mark recaptures
#' @export
#' @name check_id_match
#' @author Michelle Stuart
#' @param x = table_name
#' @examples 
#' temp <- check_id_match(row_number, table_of_matches)

check_id_match <- function(row_number, table_of_matches){
  x <- table_of_matches %>% 
    filter(first_id == table_of_matches$first_id[row_number])
  first <- x %>% 
    select(contains("first"))
  second <- x %>% 
    select(contains("second"))
  names(first) <- substr(names(first), 7, 20)
  names(second) <- substr(names(second), 8, 25)
  y <- rbind(first, second) %>% 
    distinct() %>% 
    mutate(year = substr(sample_id, 5,6), 
           size = as.numeric(size)) %>% 
    arrange(year)
  return(y)
}

# create_genid ####
#' add a gen_id to a sample
#' @export
#' @name create_genid
#' @author Michelle Stuart
#' @param x = table_name
#' @examples 
#' updated_table <- create_genid(table_of_candidates)

create_genid <- function(table_of_candidates){
  leyte <- write_db("Leyte")
  fish <- leyte %>% 
    tbl("clownfish") %>% 
    collect()
  dbDisconnect(leyte)
  
  max_gen <- fish %>% 
  summarise(max = max(gen_id, na.rm = T)) %>% 
  collect()

table_of_candidates <- table_of_candidates %>% 
  mutate(gen_id = ifelse(sample_id %in% table_of_candidates$sample_id, max_gen$max + 1, gen_id))

return(table_of_candidates)
  }



# change_db_gen_id ####
#' change the gen_id field in the database for all samples in a table
#' @export
#' @name change_db_gen_id
#' @author Michelle Stuart
#' @param x = table_of_candidates
#' @examples 
#' change_db_gen_id(table_of_candidates)

change_db_gen_id <- function(table_of_candidates){
  # backup db
  ley <- write_db("Leyte")
  fish <- ley %>% 
    tbl("clownfish") %>% 
    collect()
  write_csv(fish, path = paste0("../db_backups/", Sys.time(), "_clownfish_db.csv"))
  
  # make change
  fish <- fish %>%
    mutate(gen_id = ifelse(sample_id %in% table_of_candidates$sample_id, table_of_candidates$gen_id, gen_id))
  
  # write change
  dbWriteTable(ley, "clownfish", fish, row.names = F, overwrite = T)
  dbDisconnect(ley)
}



