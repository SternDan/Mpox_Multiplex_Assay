importELISA <- function(experimentIn, fileIn, dateIn, plateIn, isotypeIn, skipIn){
  import(fileIn,
         skip = skipIn) %>%
    pivot_longer(cols = -`<>`, names_to = "col", values_to = "data") %>%
    mutate(experiment = experimentIn,
           filename = fileIn,
           date = as_datetime(dateIn, format = "%y%m%d"),
           well = paste(`<>`, col, sep = ""),
           plate = plateIn, 
           isotype = isotypeIn) %>%
    dplyr::select(experiment, filename, date, well, plate, isotype, data)
}