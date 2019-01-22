# The clownfish package for all of your clownfish project importing and wrangling needs.

### One useful function is fish_anem_dive which will create a selectable, filterable table of just the fish you are looking for in seconds.  

Example:
```r
fish_2018 <- fish_anem_dive() %>% 
  filter(grepl("2018", date)) %>% 
  select(fish_table_id, sample_id, gen_id, recap, tag_id, anem_obs)
```

### Another useful function is sample_latlon which generates a table of sample_ids with associated map coordinates.  Here we are using the fish_table generated above

Example:
```r
fish_2018 <- filter(fish_2018, !is.na(sample_id))

coords <- sample_latlon(fish2018$sample_id)
```