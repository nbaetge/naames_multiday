N3\_Depth\_Profiles
================
Nicholas Baetge
12/28/2021

# Intro

Here, we plot the profiles from the NAAMES stations N3 S3, 3.5, 4, and
6.

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
  filter(Cruise == "AT38", Station %in% c(3, 3.5, 4, 6), z > 0, z <= 200, !plot_date %in% c("Sep 8 03:08",  "Sep 8 15:30", "Sep 11 03:07", "Sep 14 15:33", "Sep 15 03:04", "Sep 15 15:20", "Sep 16 03:04", "Sep 16 04:50", "Sep 16 07:26")) %>% 
  mutate(name = ifelse(Cruise == "AT34", "NAAMES 2", "NAAMES 3"), 
         plot_name = paste(name, "Station", Station),
         plot_sep = ifelse(Station == "6", "Subpolar", "Subtropical"))
  
ctd <- read_rds("~/GITHUB/naames_multiday/Input/ctd_data.rds") %>% 
   filter(Cruise == "AT38", Station %in% c(3, 3.5, 4, 6), z > 0, z <= 200, !plot_date %in% c("Sep 8 03:08", "Sep 8 15:30", "Sep 11 03:07", "Sep 14 15:33", "Sep 15 03:04", "Sep 15 15:20", "Sep 16 03:04", "Sep 16 04:50", "Sep 16 07:26"))  %>% 
  mutate(name = ifelse(Cruise == "AT34", "NAAMES 2", "NAAMES 3"), 
         plot_name = paste(name, "Station", Station),
         plot_sep = ifelse(Station == "6", "Subpolar", "Subtropical"))

npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
   filter(Cruise == "AT38", Station %in% c(3, 3.5, 4, 6), z > 0, z <= 200, !plot_date %in% c("Sep 8 03:08",  "Sep 8 15:30", "Sep 11 03:07", "Sep 14 15:33", "Sep 15 03:04", "Sep 15 15:20", "Sep 16 03:04", "Sep 16 04:50", "Sep 16 07:26"))  %>% 
  mutate(name = ifelse(Cruise == "AT34", "NAAMES 2", "NAAMES 3"), 
         plot_name = paste(name, "Station", Station),
         plot_sep = ifelse(Station == "6", "Subpolar", "Subtropical"))

phyto <- read_rds("~/GITHUB/naames_multiday/Input/phyto_data.rds") %>% 
   filter(Cruise == "AT38", station %in% c("S3", "S3.5", "S4", "S6"), depth > 0, depth <= 200, !plot_date %in% c("Sep 8 03:08",  "Sep 8 15:30", "Sep 11 03:07", "Sep 14 15:33", "Sep 15 03:04", "Sep 15 15:20", "Sep 16 03:04", "Sep 16 04:50", "Sep 16 07:26")) %>% 
  mutate(station = gsub("S", "", station),
         station = as.character(station),
         name = ifelse(Cruise == "AT34", "NAAMES 2", "NAAMES 3"), 
         plot_name = paste(name, "Station", station),
         plot_sep = ifelse(station == "S6", "Subpolar", "Subtropical"))

custom.shapes <- c("NAAMES 2 Station 4" = 21, "NAAMES 3 Station 3" = 22, "NAAMES 3 Station 3.5" = 23, "NAAMES 3 Station 4" = 24, "NAAMES 3 Station 6" = 25) 
```

``` r
n2_bf <- read_rds("~/GITHUB/naames_multiday/Input/bottle_data.rds") %>% 
  filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 

n2_phyto <- read_rds("~/GITHUB/naames_multiday/Input/phyto_data.rds") %>% 
   filter(Cruise == "AT34", station == "S4", !plot_date %in% c("May 24 02:30", "May 27 06:07"), depth > 0, depth <= 200) 
  
n2_ctd <- read_rds("~/GITHUB/naames_multiday/Input/ctd_data.rds") %>% 
   filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 


n2_npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
   filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 
```

``` r
# bf %>% filter(z <= 57) %>% select(leu_incorp) %>% summary
```

``` r
# npp %>% filter(z <= 57) %>% select(npp) %>% mutate(sd = sd(npp, na.rm = T)) %>% 
#   select(sd) %>% unique()
```

``` r
 # n2_npp %>% filter(z <= 77) %>% select(npp) %>% summary
```

``` r
# n2_npp %>% filter(z <= 77) %>% select(npp) %>% mutate(sd = sd(npp, na.rm = T)) %>% 
#   select(sd) %>% unique()
```

``` r
# npp %>% filter(z >= 57) %>% select(npp) %>% summary
```

``` r
# npp %>% filter(z >= 57) %>% select(npp) %>% mutate(sd = sd(npp, na.rm = T)) %>% 
#   select(sd) %>% unique()
```

``` r
 # n2_bf %>% filter(z >= 77, Date == "2016-05-26") %>% select(leu_incorp) %>% summary
```

``` r
# n2_bf %>% filter(z >= 77, Date == "2016-05-26") %>% select(leu_incorp) %>%  mutate(sd = sd(leu_incorp, na.rm = T)) %>%
#   select(sd) %>% unique()
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

![](N3_Depth_Profiles_files/figure-gfm/combine%20plots-1.png)<!-- -->
