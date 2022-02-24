S4\_Depth\_Profiles\_CORE
================
Nicholas Baetge
2/25/2021

# Intro

Here, we plot the profiles from the multiday NAAMES station, N2S4.

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
  filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 
  
ctd <- read_rds("~/GITHUB/naames_multiday/Input/ctd_data.rds") %>% 
   filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 

npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
   filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 

phyto <- read_rds("~/GITHUB/naames_multiday/Input/phyto_data.rds") %>% 
   filter(Cruise == "AT34", station == "S4", !plot_date %in% c("May 24 02:30", "May 27 06:07"), depth > 0, depth <= 200) 

int_bf <- read_rds("~/GITHUB/naames_multiday/Output/integrated_bf.rds") %>% 
  filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07")) 
```

## Plot Profiles

### Chl

### Phyto Cells

### NPP

### N

### DOC

### TDAA

### AOU

### BactA

### Leu

## Combine Plots

![](S4_Depth_Profiles_Core_files/figure-gfm/combine%20plots-1.png)<!-- -->
