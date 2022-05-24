S4\_Depth\_Profiles\_CORE
================
Nicholas Baetge
2/25/2021

# Intro

Here, we plot the profiles from the multiday NAAMES station, N2S4.

``` r
library(tidyverse)
library(lubridate)
library(hms)
library(zoo) 
library(oce)  
library(ggpubr)
library(patchwork)
```

# Import Data

These casts have been excluded from further analysis (see Float and Ship
T-S analysis): \!plot\_date %in% c(“May 24 02:30”, “May 27 06:07”)

``` r
all <- read_rds("~/GITHUB/naames_multiday/Input/bottle_data.rds") %>% 
  filter(Cruise == "AT34", z >= 0, z <= 200) 
 
bf <- read_rds("~/GITHUB/naames_multiday/Input/bottle_data.rds") %>% 
  filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 
  
ctd <- read_rds("~/GITHUB/naames_multiday/Input/ctd_data.rds") %>% 
   filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 

npp <- read_rds("~/GITHUB/naames_multiday/Input/npp_data.rds") %>% 
   filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07"), z > 0, z <= 200) 

phyto <- read_rds("~/GITHUB/naames_multiday/Input/phyto_data.rds") %>% 
   filter(Cruise == "AT34", station == "S4", !plot_date %in% c("May 24 02:30", "May 27 06:07"), depth > 0, depth <= 200) 

int_bf <- read_rds("~/GITHUB/naames_multiday/Output/integrated_bf.rds") %>% 
  filter(Cruise == "AT34", Station == 4, !plot_date %in% c("May 24 02:30", "May 27 06:07")) 
```

## Plot Profiles

### Chl

``` r
bf %>% 
  mutate(dh = ifelse(z <= 75, "ez", "mz")) %>% 
  drop_na(chl) %>% 
  filter( CampCN %in% c(53, 67)) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(chl ~ CampCN, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh    estimate estimate1 estimate2 statistic p.value parameter conf.low
    ##   <chr>    <dbl>     <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>
    ## 1 ez      -0.648     0.505     1.15      -2.44  0.0586      5.02  -1.33  
    ## 2 mz       0.160     0.435     0.275      4.09  0.0170      3.75   0.0486
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

### Phyto Cells

``` r
phyto %>% 
  mutate(dh = ifelse(depth <= 75, "ez", "mz")) %>% 
  drop_na(phyto) %>% 
  filter( plot_date %in% c("May 24 06:10", "May 27 14:56")) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(phyto ~ plot_date, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh     estimate estimate1 estimate2 statistic p.value parameter   conf.low
    ##   <chr>     <dbl>     <dbl>     <dbl>     <dbl>   <dbl>     <dbl>      <dbl>
    ## 1 ez    -32837000 14590600  47427600      -2.12 0.101        4.00 -75774502.
    ## 2 mz      4592000 14500667.  9908667.      9.66 0.00176      3.21   3134656.
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

### NPP

``` r
npp %>% 
  mutate(dh = ifelse(z <= 75, "ez", "mz")) %>% 
  drop_na(npp) %>% 
  filter( plot_date %in% c("May 24", "May 27")) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(npp ~ plot_date, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh    estimate estimate1 estimate2 statistic  p.value parameter conf.low
    ##   <chr>    <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>    <dbl>
    ## 1 ez      0.165     1.25    1.08         0.988 3.26e- 1      91.9  -0.167 
    ## 2 mz      0.0538    0.0544  0.000571     7.29  3.24e-11     124.    0.0392
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

### N

``` r
bf %>% 
  mutate(dh = ifelse(z <= 75, "ez", "mz")) %>% 
  drop_na(n) %>% 
  filter( CampCN %in% c(52, 61)) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(n ~ CampCN, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh    estimate estimate1 estimate2 statistic p.value parameter conf.low
    ##   <chr>    <dbl>     <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>
    ## 1 ez       0.202      5.28      5.08      3.09  0.0313      4.49   0.0281
    ## 2 mz      -0.157      5.30      5.46     -2.06  0.169       2.10  -0.469 
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

### DOC

``` r
bf %>% 
  mutate(dh = ifelse(z <= 75, "ez", "mz")) %>% 
  drop_na(doc) %>% 
  filter( CampCN %in% c(52, 61)) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(doc ~ CampCN, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh    estimate estimate1 estimate2 statistic p.value parameter conf.low
    ##   <chr>    <dbl>     <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>
    ## 1 ez       0.380      53.6      53.2      1.68  0.141       6.51  -0.164 
    ## 2 mz       0.433      53.5      53.1      4.91  0.0390      2      0.0539
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

### TDAA

``` r
bf %>% 
  mutate(dh = ifelse(z <= 75, "ez", "mz")) %>% 
  drop_na(tdaa) %>% 
  filter( CampCN %in% c(52, 61)) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(tdaa ~ CampCN, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh    estimate estimate1 estimate2 statistic p.value parameter conf.low
    ##   <chr>    <dbl>     <dbl>     <dbl>     <dbl>   <dbl>     <dbl>    <dbl>
    ## 1 ez     -0.0406     0.156     0.197    -2.16   0.0634      7.96  -0.0841
    ## 2 mz      0.0319     0.177     0.145     0.895  0.421       4.00  -0.0671
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

### AOU

``` r
ctd %>% 
  mutate(dh = ifelse(z <= 75, "ez", "mz")) %>% 
  drop_na(aou) %>% 
  filter( plot_date %in% c("May 24 04:11", "May 27 20:16")) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(aou ~ plot_date, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh    estimate estimate1 estimate2 statistic  p.value parameter conf.low
    ##   <chr>    <dbl>     <dbl>     <dbl>     <dbl>    <dbl>     <dbl>    <dbl>
    ## 1 ez        5.57      16.9      11.3      7.57 1.03e-10      71.2     4.10
    ## 2 mz       -1.02      17.4      18.4    -34.9  1.45e-97     247.     -1.08
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

### BactA

``` r
bf %>% 
  mutate(dh = ifelse(z <= 75, "ez", "mz")) %>% 
  drop_na(ba) %>% 
  filter( CampCN %in% c(52, 61)) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(ba ~ CampCN, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh       estimate  estimate1   estimate2 statistic p.value parameter  conf.low
    ##   <chr>       <dbl>      <dbl>       <dbl>     <dbl>   <dbl>     <dbl>     <dbl>
    ## 1 ez    -433600000  1193000000 1626600000      -2.15 0.0919       4.38   -9.75e8
    ## 2 mz    -794333333.  710000000 1504333333.     -6.68 0.00708      2.97   -1.18e9
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

### Leu

``` r
bf %>% 
  mutate(dh = ifelse(z <= 75, "ez", "mz")) %>% 
  drop_na(bp) %>% 
  filter( CampCN %in% c(52, 61)) %>% 
  group_by(dh) %>% do(broom::tidy(t.test(bp ~ CampCN, data = .)))
```

    ## # A tibble: 2 × 11
    ## # Groups:   dh [2]
    ##   dh    estimate estimate1 estimate2 statistic     p.value parameter conf.low
    ##   <chr>    <dbl>     <dbl>     <dbl>     <dbl>       <dbl>     <dbl>    <dbl>
    ## 1 ez       -30.7      21.6      52.3    -17.7  0.000000108      8.00    -34.7
    ## 2 mz       -15.9      19.2      35.2     -3.78 0.0205           3.89    -27.8
    ## # … with 3 more variables: conf.high <dbl>, method <chr>, alternative <chr>

![](S4_Depth_Profiles_Core_files/figure-gfm/bp%20plot2-1.png)<!-- -->

## Combine Plots

![](S4_Depth_Profiles_Core_files/figure-gfm/combine%20plots-1.png)<!-- -->

![](S4_Depth_Profiles_Core_files/figure-gfm/combine%20plots2-1.png)<!-- -->
