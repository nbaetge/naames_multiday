Process Data
================
Nicholas Baetge
8/4/2020

# Intro

Here, data are processed such that depth-resolved bacterial carbon and
specific growth rates are estimated. Integrated and depth-normalized
values are also calculated for the 0-100, 100-200, and 200-300 m bins.

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
    theme(legend.position = "top",
          legend.spacing.x = unit(0.5,"cm"),
          legend.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA)) 
}

custom.colors <- c("AT39" = "#377EB8", "AT34" = "#4DAF4A", "AT38" = "#E41A1C", "AT32" = "#FF7F00", "Temperate" = "#A6CEE3", "Subpolar" = "#377EB8", "Subtropical" = "#FB9A99", "GS/Sargasso" = "#E41A1C", "Early Spring" = "#377EB8", "Late Spring" = "#4DAF4A","Early Autumn" = "#E41A1C", "Summer" = "#E41A1C", "Late Autumn" = "#FF7F00", "Gv2_2019" = "#377EB8", "WOA18_MN" = "#4DAF4A", "WOA18_AN" = "#E41A1C")

levels = c("GS/Sargasso", "Subtropical", "Temperate", "Subpolar",  "AT39-6", "AT34", "AT38", "AT32","South", "North", "Early Spring", "Late Spring","Early Autumn",  "Summer", "Late Autumn", "Gv2_2019", "WOA18_MN", "WOA18_AN","Nov", "Nov sd", "Dec", "Dec sd", "Jan", "Jan sd", "Feb", "Feb sd", "Mar", "Mar sd", "Apr", "Apr sd",  "Cruise", "ARGO")

bar.colors <- c("100 m" = "white", "CM" = "#4DAF4A",  "PAM" = "#377EB8")
```

# Import Data

``` r
bf <- read_rds("~/GITHUB/naames_multiday/Input/export_ms/processed_bf.2.2020.rds") %>%  
  select(Cruise:degree_bin, CampCN, Z_MLD, Target_Z, interp_DOC, DOC_sd, interp_O2_Winkler, O2_Winkler_sd, interp_N_N, N_N_sd, interp_Chl_a_Fluor,interp_Pro_Influx:interp_Nano_Influx, interp_BactProd_C, BactProd_C_sd, interp_BactAbund, BactAbund_sd) %>%
  rename(bp = interp_BactProd_C, 
         sd_bp = BactProd_C_sd, 
         ba = interp_BactAbund,
         sd_ba = BactAbund_sd,
         o2 = interp_O2_Winkler,
         sd_o2 = O2_Winkler_sd,
         n = interp_N_N,
         sd_n = N_N_sd,
         doc = interp_DOC,
         sd_doc = DOC_sd, 
         chl = interp_Chl_a_Fluor,
         pro = interp_Pro_Influx,
         syn = interp_Syn_Influx,
         pico = interp_Pico_Influx,
         nano = interp_Nano_Influx) %>% 
  mutate(degree_bin = ifelse(Cruise == "AT34" & Station == 4, 48, degree_bin),
         phyto = pro + syn + pico + nano,
         rel.pro = pro/phyto,
         rel.syn = syn/phyto,
         rel.pico = pico/phyto,
         rel.nano = nano/phyto) %>% 
  mutate_all(~ replace(., is.nan(.), NA)) %>% 
  left_join(., read_csv("~/GITHUB/naames_multiday/Input/phytoC.csv") %>% 
              select(Cruise, Station, CampCN, Target_Z, PhytoC, sd_PhytoC)) %>% 
  #convert phytoC from ug C / L to µmol C / m3
  mutate(PhytoC = PhytoC * (10^3/12),
         sd_PhytoC = sd_PhytoC * (10^3/12))

ctd <- read_rds("~/GITHUB/naames_multiday/Input/ctd/deriv_naames_ctd.rds")

bge <- read_rds("~/GITHUB/naames_multiday/Input/bioav_ms/processed_bge.rds")

#units for npp are mg C / d, we'll convert to µmol C / d
npp <- read_rds("~/GITHUB/naames_multiday/Input/Z_resolved_model_NPP.rds") %>% 
  rename(z = depth,
         npp = NPP) %>% 
  mutate(npp = round((npp * 10^3)/12))
```

# Interpolate Phytoplankton Carbon

``` r
#split the df by CampCN  
bf.list <- split(bf, bf$CampCN)


#create a function that will linearly interpolate each VOI according to the depth intervals of the casts 
interpolate.func <- function(casper) {
to_interpolate.df <- casper %>% 
  select(Target_Z:ncol(.)) %>% 
  zoo(., order.by = .$Target_Z) 
interp_PhytoC <- as.numeric(na.approx(to_interpolate.df$PhytoC, na.rm = F))
Target_Z <- to_interpolate.df$Target_Z
interpolations.df <- data.frame(Target_Z, interp_PhytoC)
}

#apply function to list 
interpolations.list <- lapply(bf.list, interpolate.func)

#save the list as a data frame 
interpolations.df <- plyr::ldply(interpolations.list, data.frame) %>% 
  rename(., CampCN = .id) %>% 
  group_by(CampCN) %>% 
  fill(interp_PhytoC, .direction = "up") %>% 
  ungroup()

#combine the interpolated and non-interpolated data frames
interpolations.df$CampCN <- as.numeric(interpolations.df$CampCN)
interpolated.df <- left_join(bf, interpolations.df) %>% 
  select(-PhytoC) %>% 
  rename(phyc = interp_PhytoC,
         sd_phyc = sd_PhytoC)
```

# Estimate BCD and Specific Growth Rates

We’ll estimate BCD for each cast of each station of all the cruises
using bacterial production and the seasonal averages of BGE. Station
specific BGEs are not used because they could not be computed for N2S4,
where ∆DOC was NR. We’ll also estimate specific growth rates (mew) by
first converting bacterial abundance to bacterial carbon using seasonal
averages of CCFs (from POC taken from initial population) and then
dividing bacterial production by bacterial carbon.

Because there are no BGEs and CCFs associated with AT32, we will apply
those from AT38 to that cruise.

BCD and specific growth rates will be integrated for the 0-100 m, the
100-200, and the 200-300 m boxes.

``` r
#calculate bcd and specific growth rates
bcd <- interpolated.df %>% 
  full_join(.,   bge %>%
              select(Season, season_bge,  ave_initial_ccf) %>%
              distinct()) %>% 
  mutate(Season2 = ifelse(Season %in% c("Early Autumn", "Late Autumn"), "Autumn", "Spring")) %>% 
  group_by(Season2) %>%
  fill(season_bge:ave_initial_ccf, .direction = "updown") %>% 
  ungroup() %>% 
  select(-Season2) %>% 
  mutate_at(vars(ave_initial_ccf), round) %>% 
  rename(ccf = ave_initial_ccf) %>% 
  mutate(###########BCD##################
         #units are  µmol C / m3 / d 
         #bp in nmol C / L / d is equivalent to  µmol C / m3 / d 
         bcd = round(bp/season_bge),
         
         ##########Bacterial Carbon###########
         # bact carbon in µmol C/m3
         bc = ba * (10^3) * ccf * (1/10^15) * (1/12) * 10^6,
         sd_bc = sd_ba * (10^3) * ccf * (1/10^15) * (1/12) * 10^6,
        
          ##########Specific growth rates###########
         mew = bp/bc,
         mew.error = sqrt(  ((sd_bp/bp)^2) + ((sd_bc/bc)^2)  )
         ) %>% 
  rename(z = Target_Z,
         mld = Z_MLD)
```

# Integrations

``` r
integ_100 <- bcd %>% 
  group_by(CampCN) %>% 
  filter(z <= 100) %>% 
  mutate(bcd.100 = integrateTrapezoid(z, bcd, type = "A"),
         bc.100 = integrateTrapezoid(z, bc, type = "A"),
         doc.100 = integrateTrapezoid(z, doc, type = "A"),
         phyc.100 = integrateTrapezoid(z, phyc, type = "A")) %>% 
  mutate_at(vars(contains(".100")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".100")), funs(./100)) %>% 
  select(Cruise, Station, CampCN, contains(".100")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>%
  distinct() %>% 
  ungroup()
```

    ## Warning: funs() is soft deprecated as of dplyr 0.8.0
    ## Please use a list of either functions or lambdas: 
    ## 
    ##   # Simple named list: 
    ##   list(mean = mean, median = median)
    ## 
    ##   # Auto named with `tibble::lst()`: 
    ##   tibble::lst(mean, median)
    ## 
    ##   # Using lambdas
    ##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## This warning is displayed once per session.

``` r
integ_200 <- bcd %>% 
  group_by(CampCN) %>% 
  filter(between(z, 100, 200)) %>% 
  mutate(bcd.200 = integrateTrapezoid(z, bcd, type = "A"),
         bc.200 = integrateTrapezoid(z, bc, type = "A"),
         doc.200 = integrateTrapezoid(z, doc, type = "A"),
         phyc.200 = integrateTrapezoid(z, phyc, type = "A")) %>% 
  mutate_at(vars(contains(".200")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".200")), funs(./100)) %>% 
  select(Cruise, Station, CampCN, contains(".200")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>% 
  distinct() %>% 
  ungroup()

integ_300 <- bcd %>% 
  group_by(CampCN) %>% 
  filter(between(z, 200, 300)) %>% 
  mutate(bcd.300 = integrateTrapezoid(z, bcd, type = "A"),
         bc.300 = integrateTrapezoid(z, bc, type = "A"),
         doc.300 = integrateTrapezoid(z, doc, type = "A"),
         phyc.300 = integrateTrapezoid(z, phyc, type = "A")) %>% 
  mutate_at(vars(contains(".300")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".300")), funs(./100)) %>% 
  select(Cruise, Station, CampCN, contains(".300")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>% 
  distinct() %>% 
  ungroup()

npp_100 <- npp %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(z <= 100) %>% 
  mutate(npp.100 = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".100")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".100")), funs(./100)) %>% 
  select(Cruise, Station, Date, contains(".100"))  %>%
  distinct() %>% 
  ungroup()

npp_200 <- npp %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(between(z, 100, 200)) %>% 
  mutate(npp.200 = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".200")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".200")), funs(./100)) %>% 
  select(Cruise, Station, Date, contains(".200"))  %>%
  distinct() %>% 
  ungroup()

npp_300 <- npp %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(between(z, 201, 300)) %>% 
  mutate(npp.300 = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".300")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".300")), funs(./100)) %>% 
  select(Cruise, Station, Date, contains(".300"))  %>%
  distinct() %>% 
  ungroup()



integrated.df <- left_join(bcd, integ_100) %>% 
  left_join(., integ_200) %>% 
  left_join(., integ_300) %>% 
  left_join(., npp_100) %>% 
  left_join(., npp_200) %>% 
  left_join(., npp_300)
```

# Save Data

``` r
saveRDS(integrated.df, "~/GITHUB/naames_multiday/Output/processed_data.rds")
```
