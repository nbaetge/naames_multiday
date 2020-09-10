S6\_Depth\_Profiles
================
Nicholas Baetge
8/7/2020

# Intro

Here, we plot the profiles from the multiday NAAMES station, N3S6.

``` r
library(tidyverse) 
library(rmarkdown)
library(knitr)
library(readxl)
library(data.table) 
library(scales)
library(zoo)
library(oce)
library(patchwork)
#rmarkdown tables
library(stargazer)
library(pander)
#stat tests
library(lmtest)
library(lmodel2)
library(rstatix)
library(ggpubr)
#for odv type plots
library(lubridate)
library(reshape2)
library(MBA)
library(mgcv)
```

``` r
custom_theme <- function() {
  theme_test(base_size = 30) %+replace%
    theme(legend.position = "right",
          legend.spacing.x = unit(0.5,"cm"),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA)) 
}

custom.colors <- c("AT39" = "#377EB8", "AT34" = "#4DAF4A", "AT38" = "#E41A1C", "AT32" = "#FF7F00", "Temperate" = "#A6CEE3", "Subpolar" = "#377EB8", "Subtropical" = "#FB9A99", "GS/Sargasso" = "#E41A1C", "Early Spring" = "#377EB8", "Late Spring" = "#4DAF4A","Early Autumn" = "#E41A1C", "Summer" = "#E41A1C", "Late Autumn" = "#FF7F00", "Gv2_2019" = "#377EB8", "WOA18_MN" = "#4DAF4A", "WOA18_AN" = "#E41A1C")

levels = c("GS/Sargasso", "Subtropical", "Temperate", "Subpolar",  "AT39-6", "AT34", "AT38", "AT32","South", "North", "Early Spring", "Late Spring","Early Autumn",  "Summer", "Late Autumn", "Gv2_2019", "WOA18_MN", "WOA18_AN","Nov", "Nov sd", "Dec", "Dec sd", "Jan", "Jan sd", "Feb", "Feb sd", "Mar", "Mar sd", "Apr", "Apr sd",  "Cruise", "ARGO")

bar.colors <- c("100 m" = "white", "CM" = "#4DAF4A",  "PAM" = "#377EB8")

odv.colors <- c("#feb483", "#d31f2a", "#ffc000", "#27ab19", "#0db5e6", "#7139fe", "#d16cfa")
```

# Import Data

``` r
data <- read_rds("~/GITHUB/naames_multiday/Output/processed_data.rds") %>% 
  filter(Cruise == "AT34" & Station == 4 | Cruise == "AT38" & Station == 6) %>% 
  mutate_at(vars(contains("tdaa"), Asp:Lys), funs(. / 10^3)) %>% #nM to mmol/m3
  mutate(tdaa_yield = round((tdaa_c/doc)*100, 2)) 
```

    ## Warning: funs() is soft deprecated as of dplyr 0.8.0
    ## Please use a list of either functions or lambdas: 
    ## 
    ##   # Simple named list: 
    ##   list(mean = mean, median = median)
    ## 
    ##   # Auto named with `tibble::lst()`: 
    ##   tibble::lst(mean, median)
    ## 
    ##   # Using lambdas
    ##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## This warning is displayed once per session.

``` r
npp <- read_rds("~/GITHUB/naames_multiday/Input/Z_resolved_model_NPP.rds") %>% 
  rename(z = depth,
         npp = NPP) %>% 
  mutate(npp = round((npp * 10^3)/12)) %>% 
  filter(Cruise == "AT34" & Station == 4 | Cruise == "AT38" & Station == 6) 

ctd <-  read_rds("~/GITHUB/naames_multiday/Input/ctd/deriv_naames_ctd.rds") %>% 
  rename(lat = "Latitude [degrees_north]",
         z = bin_depth) %>% 
  mutate(bin = round(lat, 1),
         Date = ymd_hm(`yyyy-mm-ddThh:mm:ss.sss`),
         Date = as_date(round_date(Date, unit = "day"))) %>% 
  filter(Cruise == "AT34" & Station == 4 | Cruise == "AT38" & Station == 6)

casts <- data %>% 
  filter(Cruise == "AT38") %>% 
  distinct(CampCN) %>% 
  as_vector()
```

# Station 6

## Plot MLDs

<img src="S6_Depth_Profiles_files/figure-gfm/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

## Plot T/S

<img src="S6_Depth_Profiles_files/figure-gfm/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

## Plot Profiles

### Chl

### Phyto Cells

### PhytoC

### NPP

#### Set 1

<img src="S6_Depth_Profiles_files/figure-gfm/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

### N

### DOC

### AOU

#### Set 2

<img src="S6_Depth_Profiles_files/figure-gfm/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

### BactC

### BCD

### µ

#### Set 3

<img src="S6_Depth_Profiles_files/figure-gfm/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />
