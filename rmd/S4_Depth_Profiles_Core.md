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
  mutate_at(vars(contains("tdaa"), Asp:Lys), funs(. / 10^3)) %>% #nM to mmol/m3
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

ctd <- read_rds("~/GITHUB/naames_multiday/Input/master/deriv_naames_ctd.rds") %>%
    rename(lat = "Latitude [degrees_north]",
         z = bin_depth) %>%
  mutate(bin = round(lat, 1),
         Date = ymd_hm(`yyyy-mm-ddThh:mm:ss.sss`),
         Date = as_date(round_date(Date, unit = "day")),
         aou_c = deriv_aou_umol_l * 0.72) %>%
  rename(aou = deriv_aou_umol_l)

casts <- data %>% 
  filter(Cruise == "AT34", Date != "2016-05-27") %>% 
  distinct(CampCN) %>% 
  as_vector()
```

# Station 4

# T-S

``` r
ts <- ctd %>% 
  filter(CampCN %in% c(casts), z <= 300, SCN < 13) %>%
  select(Date, SCN, ave_temp_c, pres_db, ave_sal_psu) %>% 
  mutate(potT = swTheta(salinity = ave_sal_psu, temperature = ave_temp_c, pressure = pres_db)) %>%
     mutate(interv = interval(ymd("2016-05-24"), Date),
         dur = as.duration(interv),
         days = as.numeric(dur, "days")) %>%
  mutate_at(vars(Date), as.character) %>% 
  distinct() 

ts.plot <- ts_plot(ts, temp_col = "potT", sal_col = "ave_sal_psu", color = "days", symbol_size = 2, symbol_shape = 16, color_isopyc = "black") +
  scale_color_viridis() +
  theme_classic2(18) +
  scale_x_continuous(name = "Practical Salinity", sec.axis = sec_axis(~.,name = expression(paste("Density, kg m"^-3)))) +
  scale_y_continuous(name = "Potential Temperature, ˚C") +
  guides(colour = F) 
```

    ## Scale for 'x' is already present. Adding another scale for 'x', which will
    ## replace the existing scale.

    ## Scale for 'y' is already present. Adding another scale for 'y', which will
    ## replace the existing scale.

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

# Plots for powerpoint

``` r
library(officer)
```

    ## Warning: package 'officer' was built under R version 4.0.2

``` r
p1a <- ts.plot + theme_classic2(16)
p1b <- mld.plot + theme_classic2(16) +  geom_segment(aes( y = 75, yend = 75, x = -Inf, xend = Inf), colour = "black", linetype = 3, size = 2 ) 

p1c <- n.plot + labs(x = "Depth, m") + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
p1d <- chl.plot +  labs(x = "") + guides(fill = guide_colourbar(barheight = 8, barwidth = 3, frame.colour = "black", frame.linewidth = 2,ticks.colour = "black", ticks.linewidth = 1), color = F) +  theme_classic2(16) + geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
p1e <- phyto.plot + theme_classic2(16)
p1f <- npp.plot + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
p1g <- aou.plot + theme_classic2(16)
p1h <- ba.plot + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
p1i <- leu.plot + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
p1j <- doc.plot + theme_classic2(16) +  geom_segment(aes( x = 75, xend = 75, y = -Inf, yend = Inf), colour = "black", linetype = 3, size = 2 ) 
p1k <- tdaa.plot + theme_classic2(16)
  
p1 <-  p1a + p1b +  p1c + guide_area() + p1d + p1e + p1f + p1g  + p1h + p1i + p1j  + p1k + plot_layout(guides = 'collect') + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 22)) 

p2 <- p1a + p1b 

p3 <- (p1c + p1d + p1f) / (p1j + p1i + p1h) + plot_layout(guides = 'collect') 
  
# initialize PowerPoint slide
officer::read_pptx() %>%
  # add slide ----
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(p1, ph_location(width = 16, height = 10)) %>%
  
  officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(p2, ph_location(width = 8, height = 4)) %>%
  
   officer::add_slide() %>%
  # specify object and location of object 
  officer::ph_with(p3, ph_location(width = 16, height = 8)) %>%
  
  # export slide 
  base::print(
    target = "~/Desktop/Dissertation/MS_N2S4/Presentations/profiles.pptx"
    )
```
