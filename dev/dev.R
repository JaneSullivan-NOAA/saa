library(usethis)
library(devtools)
use_r('sizeage')
use_r('data')
use_r('age_error')
use_package('dplyr')
use_package('tidytable')
use_package('RTMB')
use_pipe()
library(dplyr)
library(RTMB)

age_data = read.csv("C:/Users/Ben.Williams/Work/assessments/northern_rockfish/2022/data/raw/goa_ts_saa_age_data.csv")
age_error = read.csv("C:/Users/Ben.Williams/Work/assessments/northern_rockfish/2022/data/output/ae_model.csv")

length_data = read.csv("C:/Users/Ben.Williams/Work/assessments/northern_rockfish/2022/data/raw/goa_ts_saa_length_data.csv")
len_bins = 15:45

usethis::use_data(length_data)
usethis::use_data(age_error)
vbl(age_data, length_data, len_bins)
devtools::document()
rm(list = c("vbl"))
vbl::age_data

data("age_data")
data("length_data")
data('age_error')
data("len_bins")
# saa(age_data, length_data, len_bins)
saa_waa(age_data, length_data, age_error, len_bins=15:45, rec_age=2)
use_r('models')
use_r('utils')
use_r('saa_waa')
