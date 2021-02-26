S4 Integration Plots
================
Nicholas Baetge
8/12/2020

# Intro

Here, we plot the integrated from the multiday NAAMES station, N2S4.

``` r
library(tidyverse)
library(lubridate)
library(patchwork)
library(ggpubr)
```

``` r
custom.colors <- c("AT39" = "#377EB8", "AT34" = "#4DAF4A", "AT38" = "#E41A1C", "AT32" = "#FF7F00", "Temperate" = "#A6CEE3", "Subpolar" = "#377EB8", "Subtropical" = "#FB9A99", "GS/Sargasso" = "#E41A1C", "Early Spring" = "#377EB8", "Late Spring" = "#4DAF4A","Early Autumn" = "#E41A1C", "Late Autumn" = "#FF7F00",  "0-75 m" = "#A50026", "75-200 m" = "#313695")

levels = c("GS/Sargasso", "Subtropical", "Temperate", "Subpolar",  "AT39-6", "AT34", "AT38", "AT32","South", "North", "Early Spring", "Late Spring","Early Autumn",  "Late Autumn", "0-75 m", "75-200 m", "200-300 m")

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
  select(Cruise, Station, Date, eddy, time:days, bcd.75:aou.300) %>% 
  distinct() %>% 
  mutate(ou.75 = (aou.75 - first(aou.75)),
         ou.200 = (aou.200 - first(aou.200)),
         ou.300 = (aou.300 - first(aou.300))) %>% 
   mutate_at(vars(contains(c("tdaa", "bcd", "npp", "phyc"))), funs(. / 10^3)) #nM to mmol/m3
```

# Pivot data

``` r
pivot_phyc_data <- data %>% 
  select(Cruise:days,  phyc.75, phyc.200, phyc.300) %>% 
  pivot_longer(phyc.75:phyc.300, names_to = "depth_interval", names_prefix = "phyc.", values_to = "phyc") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_bc_data <- data %>% 
  select(Cruise:days,  bc.75, bc.200, bc.300) %>% 
  pivot_longer(bc.75:bc.300, names_to = "depth_interval", names_prefix = "bc.", values_to = "bc") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_phyto_data <- data %>% 
  select(Cruise:days,  phyto.75, phyto.200, phyto.300) %>% 
  pivot_longer(phyto.75:phyto.300, names_to = "depth_interval", names_prefix = "phyto.", values_to = "phyto") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_ba_data <- data %>% 
  select(Cruise:days, contains("ba")) %>% 
  pivot_longer(ba.75:ba.300, names_to = "depth_interval", names_prefix = "ba.", values_to = "ba") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_npp_data <- data %>% 
  select(Cruise:days, contains("npp")) %>% 
  pivot_longer(npp.75:npp.300, names_to = "depth_interval", names_prefix = "npp.", values_to = "npp") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_bcd_data <- data %>% 
  select(Cruise:days, bcd.75, bcd.200, bcd.300) %>% 
  pivot_longer(bcd.75:bcd.300, names_to = "depth_interval", names_prefix = "bcd.", values_to = "bcd") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_doc_data <- data %>% 
  select(Cruise:days, contains("doc")) %>% 
  pivot_longer(doc.75:doc.300, names_to = "depth_interval", names_prefix = "doc.", values_to = "doc") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_tdaa_data <- data %>% 
  select(Cruise:days, contains("tdaa")) %>% 
  pivot_longer(tdaa.75:tdaa.300, names_to = "depth_interval", names_prefix = "tdaa.", values_to = "tdaa") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_aou_data <- data %>% 
  select(Cruise:days, contains("aou")) %>% 
  pivot_longer(aou.75:aou.300, names_to = "depth_interval", names_prefix = "aou.", values_to = "aou") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_ou_data <- data %>% 
  select(Cruise:days, contains("ou")) %>% 
  select(-contains("aou")) %>% 
  pivot_longer(ou.75:ou.300, names_to = "depth_interval", names_prefix = "ou.", values_to = "ou") %>% 
  mutate(depth_interval = ifelse(depth_interval == 75, "0-75 m", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "75-200 m", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivoted <- left_join(pivot_phyc_data, pivot_npp_data) %>% 
  left_join(., pivot_phyto_data) %>% 
  left_join(., pivot_bc_data) %>% 
  left_join(., pivot_ba_data) %>% 
  left_join(., pivot_bcd_data) %>% 
  left_join(., pivot_doc_data) %>% 
  left_join(., pivot_tdaa_data) %>% 
  left_join(., pivot_aou_data) %>% 
  left_join(., pivot_ou_data) %>% 
  filter(depth_interval %in% c("0-75 m", "75-200 m"))
```

# Plot Data

### BCD

### DOC

### TDAA

### AOU

![](S4_Integration_Plots_files/figure-gfm/combine%20plots-1.png)<!-- -->
