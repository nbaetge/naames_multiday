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
  select(Cruise, Station, Date, ez, sd_ez, eddy, time:days, bcd.ez:aou.300) %>% 
  distinct() %>% 
  mutate(ou.ez = (aou.ez - first(aou.ez)),
         ou.200 = (aou.200 - first(aou.200)),
         ou.300 = (aou.300 - first(aou.300))) %>% 
   mutate_at(vars(contains(c("tdaa", "bcd", "bp", "phyc"))), funs(. / 10^3)) %>% #nM to mmol/m3 
  mutate( tdaaY.ez = (tdaa.ez/doc.ez) * 100,
         tdaaY.200 = (tdaa.200/doc.200) * 100)


data2 <- read_rds("~/GITHUB/naames_multiday/Output/processed_data.rds") %>% 
  filter(Cruise == "AT34" & Station == 4) %>% 
  mutate(time = ymd_hms(datetime),
         interv = interval(first(time), time),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         eddy = ifelse(Date != "2016-05-27", "Core", "Periphery")) %>% 
  filter(eddy == "Core") %>% 
  select(Cruise, Station, Date, time:days, z, contains("leu")) %>% 
  distinct() %>% 
  filter(between(z, 100, 200)) %>% 
  drop_na(leu_incorp) %>% 
  group_by(days) %>% 
  summarize_at(vars(leu_incorp), list(mean = mean, sd = sd))

data3 <- read_rds("~/GITHUB/naames_multiday/Output/processed_data.rds") %>% 
  filter(Cruise == "AT38" & Station == 6) %>% 
  mutate(time = ymd_hms(datetime),
         interv = interval(first(time), time),
         dur = as.duration(interv),
         days = as.numeric(dur, "days")) %>% 
  select(Cruise, Station, Date, time:days, z, contains("leu")) %>% 
  distinct() %>% 
  filter(between(z, 100, 200)) %>% 
  drop_na(leu_incorp) %>% 
  summarize_at(vars(leu_incorp), list(mean = mean, sd = sd))
```

# Pivot data

``` r
pivot_phyto_data <- data %>% 
  select(Cruise:days,  phyto.ez, phyto.200, phyto.300) %>% 
  pivot_longer(phyto.ez:phyto.300, names_to = "depth_interval", names_prefix = "phyto.", values_to = "phyto") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_ba_data <- data %>% 
  select(Cruise:days,  ba.ez, ba.200, ba.300) %>% 
  pivot_longer(ba.ez:ba.300, names_to = "depth_interval", names_prefix = "ba.", values_to = "ba") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_bp_data <- data %>% 
  select(Cruise:days,  bp.ez, bp.200, bp.300) %>% 
  pivot_longer(bp.ez:bp.300, names_to = "depth_interval", names_prefix = "bp.", values_to = "bp") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_npp_data <- data %>% 
  select(Cruise:days, contains("npp")) %>% 
  pivot_longer(npp.ez:npp.300, names_to = "depth_interval", names_prefix = "npp.", values_to = "npp") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_chl_data <- data %>% 
  select(Cruise:days, contains("chl")) %>% 
  pivot_longer(chl.ez:chl.300, names_to = "depth_interval", names_prefix = "chl.", values_to = "chl") %>% 
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

pivot_tdaa_data <- data %>% 
  select(Cruise:days, contains("tdaa")) %>% 
  pivot_longer(tdaa.ez:tdaa.300, names_to = "depth_interval", names_prefix = "tdaa.", values_to = "tdaa") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivot_tdaaY_data <- data %>% 
  select(Cruise:days, contains("tdaaY")) %>% 
  pivot_longer(tdaaY.ez:tdaaY.200, names_to = "depth_interval", names_prefix = "tdaaY.", values_to = "tdaaY") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))


pivot_aou_data <- data %>% 
  select(Cruise:days, contains("aou")) %>% 
  pivot_longer(aou.ez:aou.300, names_to = "depth_interval", names_prefix = "aou.", values_to = "aou") %>% 
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


pivot_n_data <- data %>% 
  select(Cruise:days, contains("n.")) %>% 
  pivot_longer(n.ez:n.300, names_to = "depth_interval", names_prefix = "n.", values_to = "n") %>% 
  mutate(depth_interval = ifelse(depth_interval == "ez", "Euphotic", depth_interval),
         depth_interval = ifelse(depth_interval == 200, "Upper Mesopelagic", depth_interval),
         depth_interval = ifelse(depth_interval == 300, "200-300 m", depth_interval))

pivoted <- left_join(pivot_phyto_data, pivot_npp_data) %>% 
  left_join(., pivot_ba_data) %>% 
  left_join(., pivot_bcd_data) %>% 
  left_join(., pivot_doc_data) %>% 
  left_join(., pivot_tdaa_data) %>% 
  left_join(., pivot_tdaaY_data) %>% 
  left_join(., pivot_aou_data) %>% 
  left_join(., pivot_ou_data) %>% 
  left_join(., pivot_n_data) %>% 
  left_join(., pivot_chl_data) %>% 
  filter(depth_interval %in% c("Euphotic", "Upper Mesopelagic"))
```

``` r
fc <- pivoted %>% 
  arrange(depth_interval, days) %>% 
  group_by(depth_interval) %>% 
  mutate(phyto_fc = (phyto - first(phyto)) / first(phyto),
         npp_fc = (npp - first(npp)) / first(npp),
         ba_fc = (ba - first(ba)) / first(ba),
         bcd_fc = (bcd - first(bcd)) / first(bcd),
         doc_fc = (doc - first(doc)) / first(doc),
         tdaa_fc = (tdaa - first(tdaa)) / first(tdaa),
         aou_fc = (aou - first(aou)) / first(aou)) %>% 
  select(Date, time, days, depth_interval, phyto:aou_fc) %>% 
  select(-c(tdaaY.ez, tdaaY.200))


fc %>% 
  group_by(depth_interval) %>% 
  drop_na(doc) %>% 
  summarise_at(vars(doc), list(mean = mean, sd = sd))
```

    ## # A tibble: 2 x 3
    ##   depth_interval     mean    sd
    ## * <chr>             <dbl> <dbl>
    ## 1 Euphotic           53.4 0.160
    ## 2 Upper Mesopelagic  53.0 0.365

# Plot Data

### Phyto

### BC

### Chl

### BCD

### DOC

### N+N

### TDAAy

## TDAA

### AOU

![](S4_Integration_Plots_Core_files/figure-gfm/combine%20plots-1.png)<!-- -->

# Plots for powerpoint

``` r
library(officer)
```

    ## Warning: package 'officer' was built under R version 4.0.2

``` r
p1a <- phyto.plot + theme_classic2(16)
p1b <- chl.plot + theme_classic2(16)
p1c <- aou.plot + theme_classic2(16) + labs(x = "Days") 

p1d <- ba.plot + theme_classic2(16)
p1e <- bcd.plot + theme_classic2(16)
p1f <- doc.plot + labs(x = "Days") + theme(legend.position = "top") + theme_classic2(16)

p1g <- tdaa.plot + theme_classic2(16)
p1h <- n.plot + labs(x = "Days") + theme_classic2(16) 
  
p1 <-  p1a + p1b + p1d + p1e + p1h + p1c + p1f  +  p1g + plot_layout(guides = 'collect', ncol = 4) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22)) 


# initialize PowerPoint slide
officer::read_pptx() %>%
  # add slide ----
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(p1, ph_location(width = 14, height = 6)) %>%
  
  # export slide 
  base::print(
    target = "~/Desktop/Dissertation/MS_N2S4/Presentations/integrated.pptx"
    )
```

``` r
p1a <- phyto.plot + theme_classic2(16)
p1b <- chl.plot + theme_classic2(16)

p1c <- doc.plot + labs(x = "Days") + theme(legend.position = "top") + theme_classic2(16)
p1d <- tdaa.plot + theme_classic2(16)

p1e <- ba.plot + labs(x = "Days") + theme_classic2(16)
p1f <- bcd.plot + labs(x = "Days") + theme_classic2(16)


  
p1 <-  p1a + p1b + p1c + p1d + p1e + p1f  + plot_spacer() + guide_area() + plot_layout(guides = 'collect', ncol = 4) 


# initialize PowerPoint slide
officer::read_pptx() %>%
  # add slide ----
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(p1, ph_location(width = 14, height = 6)) %>%
  
  # export slide 
  base::print(
    target = "~/Desktop/Dissertation/MS_N2S4/Presentations/integrated2.pptx"
    )
```
