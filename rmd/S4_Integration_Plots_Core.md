S4\_Integration\_Plots\_Core
================
Nicholas Baetge
2/25/2021

# Intro

Here, we plot the integrated from the multiday NAAMES station, N2S4.

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

``` r
bf <- read_rds("~/GITHUB/naames_multiday/Output/integrated_bf.rds") %>% 
  filter(Cruise == "AT34", !plot_date %in% c("May 24 02:30", "May 27 06:07"))

phyto <- read_rds("~/GITHUB/naames_multiday/Output/integrated_phyto.rds") %>% 
  filter(station == "S4", !plot_date %in% c("May 24 02:30", "May 27 06:07")) %>% 
  pivot_longer(c(phyto_ez, phyto_mz), names_to = "dh", values_to = "phyto" ) %>% 
  select(date, time, plot_date, dh, phyto) %>% 
  mutate(dh = ifelse(dh == "phyto_ez", "Euphotic", "Upper Mesopelagic")) %>% 
  left_join(., read_rds("~/GITHUB/naames_multiday/Output/integrated_phyto.rds") %>% 
              filter(station == "S4", !plot_date %in% c("May 24 02:30", "May 27 06:07")) %>%
              pivot_longer(c(depthNerr_ez, depthNerr_mz), names_to = "dh", values_to = "err" ) %>% 
              select(date, time, plot_date, dh, err) %>% 
              mutate(dh = ifelse(dh == "depthNerr_ez", "Euphotic", "Upper Mesopelagic"))) %>% 
  mutate(datetime = ymd_hms(paste(date, time)))
```

    ## Joining, by = c("date", "time", "plot_date", "dh")

``` r
aou <- read_rds("~/GITHUB/naames_multiday/Output/integrated_aou.rds") %>% 
  filter(Cruise == "AT34", !plot_date %in% c("May 24 02:30", "May 27 06:07"))
```

# Plot Data

## Phyto

``` r
phyto.plot <- phyto %>% 
  ggplot(aes(x = datetime, y = phyto, color = dh, fill = dh)) +
  geom_line(size = 0.5, alpha = 0.7) +
  geom_errorbar(aes(ymin = phyto - err, ymax = phyto + err),  width = 5000) +
  geom_point(color = "black", shape = 21, size = 3, alpha = 0.7) +
  labs(x = "", y = expression(paste("Phytoplankton, cells L"^"-1")), color = "", fill = "") +
  scale_fill_viridis_d(end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  guides(color = F, fill = F, shape = F, linetype = F) +
  scale_x_datetime( date_labels = c("%b %e", "%H:%M"), date_breaks = "12 hours") +
  theme_linedraw(base_size = 16) 
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

## Chl

``` r
chl.plot <- bf %>% 
  select(datetime, chl_ez, chl_mz, chl_depthNerr_ez, chl_depthNerr_mz) %>% 
  drop_na(chl_ez) %>% 
  pivot_longer(c(chl_ez, chl_mz), names_to = "dh", values_to = "chl" ) %>% 
  select(datetime, dh, chl) %>% 
  mutate(dh = ifelse(dh == "chl_ez", "Euphotic", "Upper Mesopelagic")) %>% 
  left_join(., bf %>% 
              select(datetime, chl_ez, chl_mz, chl_depthNerr_ez, chl_depthNerr_mz) %>% 
              drop_na(chl_ez) %>% 
              pivot_longer(c(chl_depthNerr_ez, chl_depthNerr_mz), names_to = "dh", values_to = "err" ) %>% 
              select(datetime, dh, err) %>% 
              mutate(dh = ifelse(dh == "chl_depthNerr_ez", "Euphotic", "Upper Mesopelagic"))) %>% 
  ggplot(aes(x = datetime, y = chl, color = dh, fill = dh)) +
  geom_line(size = 0.5, alpha = 0.7) +
  geom_errorbar(aes(ymin = chl - err, ymax = chl + err),  width = 5000) +
  geom_point(color = "black", shape = 21, size = 3, alpha = 0.7) +
  labs(x = "", y = expression(paste("Chl ",italic("a"), ", µg L"^"-1")), color = "", fill = "") +
  scale_fill_viridis_d(end = 0.7) +
  scale_color_viridis_d(end = 0.7) +
  guides(color = F, fill = F, shape = F, linetype = F) +
  scale_x_datetime( date_labels = c("%b %e", "%H:%M"), date_breaks = "12 hours") +
  theme_linedraw(base_size = 16) 
```

    ## Joining, by = c("datetime", "dh")

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

## BA

## BCD

## DOC

## N+N

## TDAA

## AOU

![](S4_Integration_Plots_Core_files/figure-gfm/combine%20plots-1.png)<!-- -->
