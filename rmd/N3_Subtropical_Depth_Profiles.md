N3\_Subtropical\_Depth\_Profiles
================
Nicholas Baetge
12/28/2021

# Intro

Here, we plot the profiles from the NAAMES station, N3S3 .

``` r
library(tidyverse)
library(lubridate)
library(hms)
library(zoo) 
library(oce)  
```

    ## Error in get(genname, envir = envir) : object 'testthat_print' not found

``` r
library(ggpubr)
library(patchwork)
```

# Import Data

These casts have been excluded from further analysis (see Float and Ship
T-S analysis):

``` r
bf <- read_rds("~/GITHUB/naames_multiday/Input/bottle_data.rds") %>% 
  filter(Cruise == "AT38", Station %in% c(3, 3.5, 4), z <= 200, !plot_date %in% c("Sep 8 03:08",  "Sep 8 15:30", "Sep 11 03:07")) %>% 
  mutate(name = ifelse(Cruise == "AT34", "NAAMES 2", "NAAMES 3"), 
         plot_name = paste(name, "Station", Station),
         plot_sep = ifelse(Station == "6", "Subpolar", "Subtropical"))
  
ctd <- read_rds("~/GITHUB/naames_multiday/Input/ctd_data.rds") %>% 
   filter(Cruise == "AT38", Station %in% c(3, 3.5, 4), z <= 200, !plot_date %in% c("Sep 8 03:08", "Sep 8 15:30", "Sep 11 03:07"))  %>% 
  mutate(name = ifelse(Cruise == "AT34", "NAAMES 2", "NAAMES 3"), 
         plot_name = paste(name, "Station", Station),
         plot_sep = ifelse(Station == "6", "Subpolar", "Subtropical"))

npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
   filter(Cruise == "AT38", Station %in% c(3, 3.5, 4), z <= 200)  %>% 
  mutate(name = ifelse(Cruise == "AT34", "NAAMES 2", "NAAMES 3"), 
         plot_name = paste(name, "Station", Station),
         plot_sep = ifelse(Station == "6", "Subpolar", "Subtropical"))

phyto <- read_rds("~/GITHUB/naames_multiday/Input/phyto_data.rds") %>% 
   filter(Cruise == "AT38", station %in% c("S3", "S3.5", "S4"), depth <= 200) %>% 
  mutate(name = ifelse(Cruise == "AT34", "NAAMES 2", "NAAMES 3"), 
         plot_name = paste(name, "Station", station),
         plot_sep = ifelse(station == "S6", "Subpolar", "Subtropical"))
```

## Plot Profiles

### Chl

### Phyto Cells

### NPP

### N

### DOC

### AOU

### BactA

### Leu

## Combine Plots

![](N3_Subtropical_Depth_Profiles_files/figure-gfm/combine%20plots-1.png)<!-- -->
