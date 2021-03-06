S4\_Depth\_Profiles\_CORE
================
Nicholas Baetge
2/25/2021

# Intro

Here, we plot the profiles from the multiday NAAMES station, N2S4.

``` r
library(tidyverse) 
library(lubridate)
library(patchwork)
library(ggpubr)
library(viridis)
```

``` r
custom.colors <- c("AT39" = "#377EB8", "AT34" = "#4DAF4A", "AT38" = "#E41A1C", "AT32" = "#FF7F00", "Temperate" = "#A6CEE3", "Subpolar" = "#377EB8", "Subtropical" = "#FB9A99", "GS/Sargasso" = "#E41A1C", "Early Spring" = "#377EB8", "Late Spring" = "#4DAF4A","Early Autumn" = "#E41A1C", "Late Autumn" = "#FF7F00")

levels = c("GS/Sargasso", "Subtropical", "Temperate", "Subpolar",  "AT39-6", "AT34", "AT38", "AT32","South", "North", "Early Spring", "Late Spring","Early Autumn",  "Late Autumn")

# matlab.colors2 <- c("#feb483", "#d31f2a", "#ffc000", "#27ab19", "#0db5e6", "#7139fe", "#d16cfa")
odv.colors <- c( "#d31f2a", "#ffc000", "#27ab19", "#0db5e6", "#7139fe", "#d16cfa")

matlab.colors2 <- c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")
matlab.colors22 <- c("#A50026", "#D73027", "#F46D43", "#FDAE61",  "#ABD9E9", "#74ADD1", "#4575B4", "#313695")
```

# Import Data

``` r
data <- read_rds("~/GITHUB/naames_multiday/Output/processed_data.rds") %>% 
  filter(Cruise == "AT34" & Station == 4) %>% 
  # mutate_at(vars(contains("tdaa"), Asp:Lys), funs(. / 10^3)) %>% #nM to mmol/m3
  mutate(time = ymd_hms(datetime),
         interv = interval(first(time), time),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         eddy = ifelse(Date != "2016-05-27", "Core", "Periphery"),
         tdaa_yield = round((tdaa_c/doc)*100, 2)) %>% 
  filter(z <= 200) %>% 
  filter(eddy == "Core")

npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
  mutate(interv = interval(ymd("2016-05-24"), Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         eddy = ifelse(Date != "2016-05-27" & Station != 6, "Core", "Periphery")) %>% 
  filter(Cruise == "AT34" & Station == 4)  %>% 
  filter(z <= 200) %>% 
  filter(eddy == "Core")

ctd <-  read_rds("~/GITHUB/naames_multiday/Input/ctd_data.rds") %>% 
  mutate(interv = interval(ymd("2016-05-24"), Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         eddy = ifelse(Date != "2016-05-27" & Station != 6, "Core", "Periphery")) %>% 
  filter(Cruise == "AT34" & Station == 4) %>% 
  filter(z <= 200) %>% 
  filter(eddy == "Core")

casts <- data %>% 
  filter(Cruise == "AT34", Date != "2016-05-27") %>% 
  distinct(CampCN) %>% 
  as_vector()
```

# Station 4

## Plot MLDs

## Plot Profiles

### Temp

### Chl

### Phyto Cells

### NPP

### N

### DOC

### TDAA

### TDAA Yield

### AOU

### BactA

### Leu

![](S4_Depth_Profiles_Core_files/figure-gfm/combine%20plots-1.png)<!-- -->
