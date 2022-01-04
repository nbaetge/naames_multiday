S6\_Depth\_Profiles
================
Nicholas Baetge
2/15/2021

# Intro

Here, we plot the profiles from the multiday NAAMES station, N3S6 .These
casts have different T-S plots than all the float profiles and 7 of the
ship profiles: plot\_date %in% c(“Sep 14 15:33”, “Sep 15 03:04”, “Sep 15
15:20”, “Sep 16 03:04”, “Sep 16 04:50”, “Sep 16 07:26”). We’ll exclude
these from these analyses.

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
T-S analysis): \!plot\_date %in% c(“May 24 02:30”, “May 27 06:07”)

``` r
bf <- read_rds("~/GITHUB/naames_multiday/Input/bottle_data.rds") %>% 
  filter(Cruise == "AT38", Station == 6, !plot_date %in% c("Sep 14 15:33", "Sep 15 03:04", "Sep 15 15:20", "Sep 16 03:04", "Sep 16 04:50", "Sep 16 07:26"), z <= 200) 
  
ctd <- read_rds("~/GITHUB/naames_multiday/Input/ctd_data.rds") %>% 
   filter(Cruise == "AT38", Station == 6, !plot_date %in% c("Sep 14 15:33", "Sep 15 03:04", "Sep 15 15:20", "Sep 16 03:04", "Sep 16 04:50", "Sep 16 07:26"), z <= 200) 

npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
   filter(Cruise == "AT38", Station == 6, !plot_date %in% c("Sep 14 15:33", "Sep 15 03:04", "Sep 15 15:20", "Sep 16 03:04", "Sep 16 04:50", "Sep 16 07:26"), z <= 200) 

phyto <- read_rds("~/GITHUB/naames_multiday/Input/phyto_data.rds") %>% 
   filter(Cruise == "AT38", station == "S6", !plot_date %in% c("Sep 14 15:33", "Sep 15 03:04", "Sep 15 15:20", "Sep 16 03:04", "Sep 16 04:50", "Sep 16 07:26"), depth <= 200) 
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

![](S6_Depth_Profiles_files/figure-gfm/combine%20plots-1.png)<!-- -->
