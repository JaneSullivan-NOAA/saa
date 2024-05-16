# size at age
1. Install package  
`remotes::install_github('BenWilliams-NOAA/saa')`  
`library(saa)`
2. Pull in example data  
`data("age_data")`  
`data("length_data")`  
`data("age_error")`  
3. run model  
`saa_waa(age_data, length_data, age_error, len_bins=15:45, rec_age=2)`  
4. profit?? 
