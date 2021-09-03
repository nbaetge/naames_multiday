S6\_Depth\_Profiles
================
Nicholas Baetge
2/15/2021

# Intro

Here, we plot the profiles from the multiday NAAMES station, N3S6.

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
  filter(Cruise == "AT38" & Station == 6) %>% 
  # mutate_at(vars(contains("tdaa"), Asp:Lys), funs(. / 10^3)) %>% #nM to mmol/m3
  mutate(time = ymd_hms(datetime),
         interv = interval(first(time), time),
         dur = as.duration(interv),
         days = as.numeric(dur, "days"),
         tdaa_yield = round((tdaa_c/doc)*100, 2)) %>% 
  filter(z <= 200) 

npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
  mutate(interv = interval(ymd("2017-09-13"), Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days")) %>% 
  filter(Cruise == "AT38" & Station == 6) %>% 
  filter(z <= 200) 

ctd <- read_rds("~/GITHUB/naames_multiday/Input/master/deriv_naames_ctd.rds") %>%
    rename(lat = "Latitude [degrees_north]",
         z = bin_depth) %>%
  mutate(bin = round(lat, 1),
         Date = ymd_hm(`yyyy-mm-ddThh:mm:ss.sss`),
         Date = as_date(round_date(Date, unit = "day")),
         aou_c = deriv_aou_umol_l * 0.72) %>%
  rename(aou = deriv_aou_umol_l)

casts <- data %>% 
  filter(Cruise == "AT38" & Station == 6) %>% 
  distinct(CampCN) %>% 
  as_vector()
```

# Station 6

# T-S

``` r
ts <- ctd %>% 
  filter(CampCN %in% c(casts), z <= 200, SCN < 13) %>%
  select(Date, SCN, ave_temp_c, pres_db, ave_sal_psu) %>% 
  mutate(potT = swTheta(salinity = ave_sal_psu, temperature = ave_temp_c, pressure = pres_db)) %>%
     mutate(interv = interval(ymd("2017-09-13"), Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days")) %>%
  mutate_at(vars(Date), as.character) %>% 
  distinct() %>% 
  mutate_at(vars(days), as.character)

ts.plot <- ts_plot(ts, temp_col = "potT", sal_col = "ave_sal_psu", color = "days", symbol_size = 3, symbol_shape = 16, symbol_alpha = 0.7, color_isopyc = "black") +
  scale_color_viridis_d() +
  theme_classic2(18) +
  scale_x_continuous(name = "Practical Salinity", sec.axis = sec_axis(~.,name = expression(paste("Potential Density, kg m"^-3)))) +
  scale_y_continuous(name = "Potential Temperature, ˚C") +
  labs(color = "Day") +
   theme(panel.spacing.x = unit(1, "cm"),
        axis.text.x = element_text(angle = 0),
        legend.title = element_text(size = 14),
        legend.key.size = unit(0.4, "cm"),
        legend.position = c(0.9, 0.35),
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

### AOU

### BactA

### Leu

## Combine Plots

![](S6_Depth_Profiles_files/figure-gfm/combine%20plots-1.png)<!-- -->

# Plots for powerpoint

``` r
library(officer)
```

    ## Warning: package 'officer' was built under R version 4.0.2

``` r
p1a <- ts.plot + theme_classic2(16)
p1b <- mld.plot + theme_classic2(16)
p1c <- n.plot + theme_classic2(16)
p1d <- chl.plot + guides(fill = guide_colourbar(barheight = 8, barwidth = 6, frame.colour = "black", frame.linewidth = 2,ticks.colour = "black", ticks.linewidth = 1), color = F) +  theme_classic2(16)
p1e <- phyto.plot + theme_classic2(16)
p1f <- npp.plot + theme_classic2(16)
p1g <- aou.plot + theme_classic2(16)
p1h <- ba.plot + theme_classic2(16) 
p1i <- leu.plot + theme_classic2(16)
p1j <- doc.plot + theme_classic2(16)
  
p1 <-  p1a + p1b +  p1c + guide_area() + p1d + p1e + p1f + p1g  + p1h + p1i + p1j  +  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22)) 
  
  
# initialize PowerPoint slide
officer::read_pptx() %>%
  # add slide ----
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(p1, ph_location(width = 16, height = 10)) %>%
  
  # export slide 
  base::print(
    target = "~/Desktop/Dissertation/MS_N2S4/Presentations/s6_profiles.pptx"
    )
```

    ## Warning: colourbar guide needs continuous scales.

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning -
    ## Inf
    
    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning -
    ## Inf
