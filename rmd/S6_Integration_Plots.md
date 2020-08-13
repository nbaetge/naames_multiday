S6\_Integration\_Plots
================
Nicholas Baetge
8/12/2020

# Intro

Here, we plot the integrated from the multiday NAAMES station, N3S6.

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

custom.colors <- c("AT39" = "#377EB8", "AT34" = "#4DAF4A", "AT38" = "#E41A1C", "AT32" = "#FF7F00", "Temperate" = "#A6CEE3", "Subpolar" = "#377EB8", "Subtropical" = "#FB9A99", "GS/Sargasso" = "#E41A1C", "Early Spring" = "#377EB8", "Late Spring" = "#4DAF4A","Early Autumn" = "#E41A1C", "Summer" = "#E41A1C", "Late Autumn" = "#FF7F00", "0 - 100 m" = "#377EB8", "100 - 200 m" = "#4DAF4A", "200 - 300 m" = "#E41A1C")

levels = c("GS/Sargasso", "Subtropical", "Temperate", "Subpolar",  "AT39-6", "AT34", "AT38", "AT32","South", "North", "Early Spring", "Late Spring","Early Autumn",  "Summer", "Late Autumn", "Gv2_2019", "WOA18_MN", "WOA18_AN","Nov", "Nov sd", "Dec", "Dec sd", "Jan", "Jan sd", "Feb", "Feb sd", "Mar", "Mar sd", "Apr", "Apr sd",  "Cruise", "ARGO")

bar.colors <- c("100 m" = "white", "CM" = "#4DAF4A",  "PAM" = "#377EB8")

odv.colors <- c("#feb483", "#d31f2a", "#ffc000", "#27ab19", "#0db5e6", "#7139fe", "#d16cfa")
```

# Import Data

``` r
data <- read_rds("~/GITHUB/naames_multiday/Output/processed_data.rds") %>% 
  filter(Cruise == "AT34" & Station == 4 | Cruise == "AT38" & Station == 6) %>% 
  filter(Cruise == "AT38") %>% 
  mutate(time = ymd_hms(datetime),
         interv = interval(first(time), time),
         dur = as.duration(interv),
         days = as.numeric(dur, "days")) %>% 
  select(Cruise, Station, Date, time:days, bcd.100:npp.300) %>% 
  distinct() %>% 
  mutate_at(vars(contains(c("phyc", "npp", "bc"))), function(x) (x / 1000)) 
```

# Pivot data

``` r
pivot_phyc_data <- data %>% 
  select(Cruise:days, contains("phyc")) %>% 
  pivot_longer(phyc.100:phyc.300, names_to = "depth_interval", names_prefix = "phyc.", values_to = "phyc") %>% 
  mutate(depth_interval = ifelse(depth_interval == 100, "0 - 100 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "100 - 200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200 - 300 m", depth_interval))


pivot_npp_data <- data %>% 
  select(Cruise:days, contains("npp")) %>% 
  pivot_longer(npp.100:npp.300, names_to = "depth_interval", names_prefix = "npp.", values_to = "npp") %>% 
  mutate(depth_interval = ifelse(depth_interval == 100, "0 - 100 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "100 - 200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200 - 300 m", depth_interval))


pivot_bc_data <- data %>% 
  select(Cruise:days, contains("bc.")) %>% 
  pivot_longer(bc.100:bc.300, names_to = "depth_interval", names_prefix = "bc.", values_to = "bc") %>% 
  mutate(depth_interval = ifelse(depth_interval == 100, "0 - 100 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "100 - 200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200 - 300 m", depth_interval))


pivot_bcd_data <- data %>% 
  select(Cruise:days, contains("bcd")) %>% 
  pivot_longer(bcd.100:bcd.300, names_to = "depth_interval", names_prefix = "bcd.", values_to = "bcd") %>% 
  mutate(depth_interval = ifelse(depth_interval == 100, "0 - 100 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "100 - 200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200 - 300 m", depth_interval))

pivoted <- left_join(pivot_phyc_data, pivot_npp_data) %>% 
  left_join(., pivot_bc_data) %>% 
  left_join(., pivot_bcd_data) %>% 
  rename(`Depth Interval` = depth_interval)
```

# Plot Data

### Phyto Carbon

### NPP

### Bact Carbon

### BCD

<img src="S6_Integration_Plots_files/figure-gfm/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

# BCD:NPP

<img src="S6_Integration_Plots_files/figure-gfm/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />
