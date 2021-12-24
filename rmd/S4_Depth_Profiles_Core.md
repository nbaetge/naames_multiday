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
library(oce)
```

    ## Error in get(genname, envir = envir) : object 'testthat_print' not found

``` r
library(PlotSvalbard)
library(hms)
library(readxl)
```

# Import Data

``` r
data <- read_rds("~/GITHUB/naames_multiday/Output/processed_data.rds") %>% 
  filter(Cruise == "AT34" & Station == 4) %>% 
  mutate_at(vars(contains("tdaa"), Asp:Lys), funs(. / 10^3)) %>% #nM to mmol/m3
  mutate(time = ymd_hms(datetime),
         plot_date = paste(month(datetime, label = T), day(datetime), format(parse_date_time(datetime, c('ymd_hms', 'HM')), '%H:%M')),
         interv = interval(first(time), time),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         eddy = ifelse(Date != "2016-05-27", "Core", "Periphery"),
         tdaa_yield = round((tdaa_c/doc)*100, 2)) %>% 
  filter(z <= 200) 
  # filter(eddy == "Core")

npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
  mutate(interv = interval(ymd("2016-05-24"), Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         eddy = ifelse(Date != "2016-05-27" & Station != 6, "Core", "Periphery")) %>% 
  filter(Cruise == "AT34" & Station == 4)  %>% 
  filter(z <= 200) %>% 
  mutate(plot_date = paste(month(Date, label = T), day(Date)))
  # filter(eddy == "Core")

ctd <- read_rds("~/GITHUB/naames_multiday/Input/master/deriv_naames_ctd.rds") %>%
    rename(lat = "Latitude [degrees_north]",
         z = bin_depth) %>%
  mutate(bin = round(lat, 1),
         Date = ymd_hm(`yyyy-mm-ddThh:mm:ss.sss`),
         Date = as_date(round_date(Date, unit = "day")),
         datetime = ymd_hm(`yyyy-mm-ddThh:mm:ss.sss`),
        plot_date = paste(month(datetime, label = T), day(datetime), format(parse_date_time(datetime, c('ymd_hms', 'HM')), '%H:%M')),
         aou_c = deriv_aou_umol_l * 0.72) %>%
  rename(aou = deriv_aou_umol_l)

casts <- data %>% 
  filter(Cruise == "AT34") %>% 
  # filter(Cruise == "AT34", Date != "2016-05-27") %>% 
  distinct(CampCN) %>% 
  as_vector()

phyto <- read_excel("~/GITHUB/naames_multiday/Input/Behrenfeld_Flow_Cytometry_NAAMES2.xlsx", sheet = "data") %>% 
  mutate(phyto = prochlorococcus_abun + synechococcus_abun + picoeukaryote_abun + nanoeukaryote_abun) %>% 
  filter(station == "S4", profile == "T", depth <= 200) %>% 
  select(-profile) %>% 
  mutate(time = as_hms(time), date = ymd(date)) %>% 
  group_by(date, time, station, lat, lon) %>% 
  group_modify(~ add_row(., depth = c(0, 5, 10))) %>% 
  fill(date, time) %>% 
  arrange(date, time, depth) %>% 
  ungroup() %>% 
  group_by(date, time, depth) %>% 
  fill(prochlorococcus_abun:phyto) %>% 
  ungroup() %>% 
  distinct() %>% 
  mutate(upper10 = ifelse(depth <= 10, T, F)) %>% 
  group_by(date, time, upper10) %>% 
  fill(prochlorococcus_abun:phyto, .direction = "up") %>% 
  ungroup() %>% 
  group_by(date, time) %>%
  mutate(phyto = zoo::na.approx(phyto, na.rm = F)) %>%
  fill(phyto, .direction = "up") %>% 
  ungroup() %>% 
  drop_na(phyto) %>% 
  mutate(plot_date = paste(month(date, label = T), day(date), format(parse_date_time(time, c('hms', 'HM')), '%H:%M')))
  

int.phyto <- phyto %>% 
  filter(depth <= 75) %>% 
  group_by(date, time) %>% 
  mutate(phyto.ez = integrateTrapezoid(depth, phyto, type = "A")) %>% 
  select(date:lon, phyto.ez) %>% 
  distinct() %>% 
  left_join(., phyto %>% 
              filter(depth > 75) %>% 
              group_by(date, time) %>% 
              mutate(phyto.mz = integrateTrapezoid(depth, phyto, type = "A")) %>% 
              select(date:lon, phyto.mz) %>% 
              distinct())

all.phyto <- left_join(phyto, int.phyto) %>% 
  write_rds("~/GITHUB/naames_multiday/Input/phyto_data.rds")
```

# Station 4

## T-S

``` r
ts <- ctd %>% 
  filter(CampCN %in% c(casts), z <= 300, SCN < 13) %>%
  select(Date, plot_date, SCN, ave_temp_c, pres_db, ave_sal_psu) %>% 
  mutate(potT = swTheta(salinity = ave_sal_psu, temperature = ave_temp_c, pressure = pres_db)) %>%
     mutate(interv = interval(ymd("2016-05-24"), Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days")) %>%
  mutate_at(vars(Date), as.character) %>% 
  distinct() %>% 
  mutate_at(vars(days), as.character)

ts.plot <- ts_plot(ts, temp_col = "potT", sal_col = "ave_sal_psu", color = "plot_date", symbol_size = 3, symbol_shape = 16, symbol_alpha = 0.7, color_isopyc = "black", WM = NULL) +
  scale_color_viridis_d() +
  theme_classic2(18) +
  scale_x_continuous(name = "Practical Salinity", sec.axis = sec_axis(~.,name = expression(paste("Potential Density, kg m"^-3)))) +
  scale_y_continuous(name = "Potential Temperature, ˚C") +
  labs(color = "Day") +
   theme(panel.spacing.x = unit(1, "cm"),
        axis.text.x = element_text(angle = 0),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.4, "cm"),
        legend.position = c(0.8, 0.2),
        legend.text = element_text(size = 12),
        legend.background = element_rect(size = 0.2, linetype = "solid", color = "black"))
```

## Plot MLDs

## Plot Profiles

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

## Combine Plots

![](S4_Depth_Profiles_Core_files/figure-gfm/combine%20plots-1.png)<!-- -->

# Plots for powerpoint

``` r
# library(officer)
```

``` r
# p1a <- ts.plot
# p1b <- mld.plot 
# p1c <- n.plot 
# p1d <- chl.plot 
# p1e <- phyto.plot + theme_classic2(16)
# p1f <- npp.plot + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
# p1g <- aou.plot + theme_classic2(16)
# p1h <- ba.plot + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
# p1i <- leu.plot + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
# p1j <- doc.plot + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
# p1k <- tdaa.plot + theme_classic2(16)
#   
# p1 <-  p1a + p1b +  p1c + guide_area() + p1d + p1e + p1f + p1g  + p1h + p1i + p1j  + p1k + plot_annotation(tag_levels = "a") &
#   theme(plot.tag = element_text(size = 22)) 
# 
# p2 <- p1a + p1b 
# 
# p3 <- (p1c + p1d + p1f) / (p1j + p1i + p1h) + plot_layout(guides = 'collect') 
#   
# # initialize PowerPoint slide
# officer::read_pptx() %>%
#   # add slide ----
#   officer::add_slide() %>%
#   # specify object and location of object 
#   officer::ph_with(p1, ph_location(width = 16, height = 10)) %>%
#   
#   officer::add_slide() %>%
#   # specify object and location of object 
#   officer::ph_with(p2, ph_location(width = 8, height = 4)) %>%
#   
#    officer::add_slide() %>%
#   # specify object and location of object 
#   officer::ph_with(p3, ph_location(width = 16, height = 8)) %>%
#   
#   # export slide 
#   base::print(
#     target = "~/Desktop/Dissertation/MS_N2S4/Presentations/profiles.pptx"
#     )
```

``` r
# all_data <- read_rds("~/GITHUB/naames_multiday/Output/processed_data.rds") %>% 
#   select(Cruise, Station, Latitude, Longitude, CampCN, mld,z, bp, sd_bp ) %>% drop_na(bp) %>% 
#   filter(z >= 100) %>% 
#   filter(Cruise == "AT38")
# 
# summary(all_data)
# sd(all_data$bp)
```
