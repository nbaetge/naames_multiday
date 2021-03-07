S4\_Integration\_Plots\_Core
================
Nicholas Baetge
2/25/2021

# Intro

Here, we plot the integrated from the multiday NAAMES station, N2S4.

``` r
library(tidyverse)
library(lubridate)
library(patchwork)
library(ggpubr)
```

``` r
custom.colors <- c("AT39" = "#377EB8", "AT34" = "#4DAF4A", "AT38" = "#E41A1C", "AT32" = "#FF7F00", "Temperate" = "#A6CEE3", "Subpolar" = "#377EB8", "Subtropical" = "#FB9A99", "GS/Sargasso" = "#E41A1C", "Early Spring" = "#377EB8", "Late Spring" = "#4DAF4A","Early Autumn" = "#E41A1C", "Late Autumn" = "#FF7F00",  "0-75 m" = "#A50026", "75-200 m" = "#313695", "Euphotic" = "#A50026", "Upper Mesopelagic" = "#313695", "Mesopelagic" = "#d16cfa")

levels = c("GS/Sargasso", "Subtropical", "Temperate", "Subpolar",  "AT39-6", "AT34", "AT38", "AT32","South", "North", "Early Spring", "Late Spring","Early Autumn",  "Late Autumn", "0-75 m", "75-200 m", "200-300 m", "Euphotic", "Upper Mesopelagic", "Mesopelagic")

odv.colors <- c("#feb483", "#d31f2a", "#ffc000", "#27ab19", "#0db5e6", "#7139fe", "#d16cfa")
```

# Import Data

``` r
data <- read_rds("~/GITHUB/naames_multiday/Output/processed_data.rds") %>% 
  filter(Cruise == "AT34" & Station == 4) %>% 
  mutate(time = ymd_hms(datetime),
         interv = interval(first(time), time),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         eddy = ifelse(Date != "2016-05-27", "Core", "Periphery")) %>% 
  filter(eddy == "Core") %>% 
  select(Cruise, Station, Date, ez, sd_ez, eddy, time:days, bcd.ez:aou_c.300) %>% 
  distinct() %>% 
  mutate(ou.ez = (aou_c.ez - first(aou_c.ez)),
         ou.200 = (aou_c.200 - first(aou_c.200)),
         ou.300 = (aou_c.300 - first(aou_c.300))) %>% 
   mutate_at(vars(contains(c("tdaa", "bcd", "bp", "phyc"))), funs(. / 10^3)) #nM to mmol/m3 
```

# Pivot data

``` r
pivot_phyc_data <- data %>% 
  select(Cruise:days,  phyc.ez, phyc.200, phyc.300) %>% 
  pivot_longer(phyc.ez:phyc.300, names_to = "depth_interval", names_prefix = "phyc.", values_to = "phyc") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_bc_data <- data %>% 
  select(Cruise:days,  bc.ez, bc.200, bc.300) %>% 
  pivot_longer(bc.ez:bc.300, names_to = "depth_interval", names_prefix = "bc.", values_to = "bc") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_npp_data <- data %>% 
  select(Cruise:days, contains("npp")) %>% 
  pivot_longer(npp.ez:npp.300, names_to = "depth_interval", names_prefix = "npp.", values_to = "npp") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_bcd_data <- data %>% 
  select(Cruise:days, bcd.ez, bcd.200, bcd.300) %>% 
  pivot_longer(bcd.ez:bcd.300, names_to = "depth_interval", names_prefix = "bcd.", values_to = "bcd") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_doc_data <- data %>% 
  select(Cruise:days, contains("doc")) %>% 
  pivot_longer(doc.ez:doc.300, names_to = "depth_interval", names_prefix = "doc.", values_to = "doc") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_aou_data <- data %>% 
  select(Cruise:days, contains("aou")) %>% 
  pivot_longer(aou_c.ez:aou_c.300, names_to = "depth_interval", names_prefix = "aou_c.", values_to = "aou_c") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_ou_data <- data %>% 
  select(Cruise:days, contains("ou")) %>% 
  select(-contains("aou")) %>% 
  pivot_longer(ou.ez:ou.300, names_to = "depth_interval", names_prefix = "ou.", values_to = "ou") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivoted <- left_join(pivot_phyc_data, pivot_npp_data) %>% 
  left_join(., pivot_bc_data) %>% 
  left_join(., pivot_bcd_data) %>% 
  left_join(., pivot_doc_data) %>% 
  left_join(., pivot_aou_data) %>% 
  left_join(., pivot_ou_data) %>% 
  filter(depth_interval %in% c("Euphotic", "Upper Mesopelagic"))
```

# Plot Data

### PhyC

### BC

### NPP

### BCD

### DOC

### AOU

![](S4_Integration_Plots_Core_files/figure-gfm/combine%20plots-1.png)<!-- -->
