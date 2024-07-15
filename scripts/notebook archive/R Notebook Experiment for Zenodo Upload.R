# zenodo large file upload
# didn't get this to work.. too large? loads into desktop memory and creates memory liability

library(httr) 
library(dplyr)
library(purrr)

# get your token here
# https://zenodo.org/account/settings/applications/
token <- "5rEkP2xggF4cixlTwiOW5BlQ3snYEGbixDi815R8GaGPjmt7BaGo8xEje0af"
deposit_id <- 12193924 # fill in UI form to get this number
file <- "C://Users//Greig//Downloads//HRDPS_2017.zip"

bucket <- GET(paste0("https://www.zenodo.org/api/deposit/depositions/",deposit_id), 
              add_headers(Authorization = paste("Bearer", token)),
              encode = 'json') %>% 
  content(as = "text") %>% 
  jsonlite::fromJSON() %>% 
  pluck("links") %>% 
  pluck("bucket") %>% 
  gsub("https://zenodo.org/api/files/","",.)
print(bucket)

PUT(url = paste0("https://www.zenodo.org/api/files/",bucket,"/",file,"?access_token=",token), 
    body = upload_file(file) # path to your local file
) %>% 
  content(as = "text")
