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
library(lubridate)
library(zoo)
library(oce)
```

# Import Data

``` r
bf <- read_rds("~/GITHUB/naames_multiday/Input/master/processed_bf.8.2020.rds") %>%  
  select(Cruise:degree_bin, CampCN, Z_MLD, Target_Z, interp_DOC, DOC_sd, interp_TDAA, interp_tdaa_c, TDAA_sd, Asp:Lys, interp_O2_Winkler, O2_Winkler_sd, interp_N_N, N_N_sd,  interp_Chl_a_Fluor,interp_Pro_Influx:interp_Nano_Influx, interp_BactProd_C, BactProd_C_sd, interp_BactAbund, BactAbund_sd) %>%
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
         tdaa = interp_TDAA,
         tdaa_c = interp_tdaa_c,
         sd_tdaa = TDAA_sd,
         chl = interp_Chl_a_Fluor,
         pro = interp_Pro_Influx,
         syn = interp_Syn_Influx,
         pico = interp_Pico_Influx,
         nano = interp_Nano_Influx) %>% 
  mutate(degree_bin = ifelse(Cruise == "AT34" & Station == 4, 48, degree_bin),
         Station = ifelse(Station == "1A", 0, Station),
         phyto = pro + syn + pico + nano,
         rel.pro = pro/phyto,
         rel.syn = syn/phyto,
         rel.pico = pico/phyto,
         rel.nano = nano/phyto) %>% 
  mutate_all(~ replace(., is.nan(.), NA)) %>% 
  mutate_at(vars(Station), as.numeric) %>% 
  left_join(., read_csv("~/GITHUB/naames_multiday/Input/master/phytoC.csv") %>% 
              select(Cruise, Station, CampCN, Target_Z, PhytoC, sd_PhytoC)) %>% 
  #convert phytoC from ug C / L to µmol C / m3
  mutate(PhytoC = PhytoC * (10^3/12),
         sd_PhytoC = sd_PhytoC * (10^3/12)) 
saveRDS(bf, "~/GITHUB/naames_multiday/Input/bottle_data.rds")

ctd <- read_rds("~/GITHUB/naames_multiday/Input/master/deriv_naames_ctd.rds") %>% 
    rename(lat = "Latitude [degrees_north]",
         z = bin_depth) %>% 
  mutate(bin = round(lat, 1),
         Date = ymd_hm(`yyyy-mm-ddThh:mm:ss.sss`),
         Date = as_date(round_date(Date, unit = "day")))
#lol, measured surface par but not par
saveRDS(ctd, "~/GITHUB/naames_multiday/Input/ctd_data.rds")

bge <- read_rds("~/GITHUB/naames_multiday/Input/master/processed_bge_summary.rds") %>% 
  select(Season, Station, ave_bge_station, ave_bge_cruise, ave_ccf_station, ave_ccf_cruise) %>% 
  distinct() %>% 
  filter(!Station %in% c("S2RD", "S2RF")) %>% 
  mutate_at(vars(Station), as.numeric)
saveRDS(bge, "~/GITHUB/naames_multiday/Input/bge_data.rds") 

#units for npp are mg C / d, we'll convert to µmol C / d
npp <- read_rds("~/GITHUB/naames_multiday/Input/master/Z_resolved_model_NPP.rds") %>% 
  rename(z = depth,
         npp = NPP) %>% 
  mutate(npp = round((npp * 10^3)/12))
saveRDS(npp, "~/GITHUB/naames_multiday/Input/npp_data.rds")
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
using bacterial production and the station averages of BGE. Otherwise,
they are calculated using the cruise average bge. We’ll also estimate
specific growth rates (mew) by first converting bacterial abundance to
bacterial carbon using station or seasonal averages of CCFs (from POC
taken from initial population) and then dividing bacterial production by
bacterial carbon.

Because there are no BGEs and CCFs associated with AT32, we will apply
those from AT38 to that cruise.

BCD and specific growth rates will be integrated for the 0-75 m, the
75-200 m, and the 200-300 m boxes.

``` r
#calculate bcd and specific growth rates
bcd <- interpolated.df %>% 
  full_join(.,   bge) %>% 
  mutate(Season2 = ifelse(Season %in% c("Early Autumn", "Late Autumn"), "Autumn", "Spring")) %>% 
  group_by(Season2) %>%
  fill(c(ave_bge_cruise, ave_ccf_cruise), .direction = "updown") %>% 
  ungroup() %>% 
  select(-Season2) %>% 
  mutate(bge = ifelse(!is.na(ave_bge_station), ave_bge_station, ave_bge_cruise), 
         ccf = ifelse(!is.na(ave_ccf_station), ave_ccf_station, ave_ccf_cruise)) %>% 
 mutate_at(vars(ccf), round) %>% 
  mutate(###########BCD##################
         #units are  µmol C / m3 / d 
         #bp in nmol C / L / d is equivalent to  µmol C / m3 / d 
         bcd = round(bp/bge),
         
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
integ_75 <- bcd %>% 
  group_by(CampCN) %>% 
  filter(z <= 75) %>% 
  mutate(bcd.75 = integrateTrapezoid(z, bcd, type = "A"),
         bc.75 = integrateTrapezoid(z, bc, type = "A"),
         doc.75 = integrateTrapezoid(z, doc, type = "A"),
         phyc.75 = integrateTrapezoid(z, phyc, type = "A"),
         tdaa.75 = integrateTrapezoid(z, tdaa_c, type = "A")) %>% 
  mutate_at(vars(contains(".75")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".75")), funs(./75)) %>% 
  select(Cruise, Station, CampCN, contains(".75")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>%
  distinct() %>% 
  ungroup()


integ_100 <- bcd %>% 
  group_by(CampCN) %>% 
  filter(between(z, 75, 100)) %>% 
  mutate(bcd.100 = integrateTrapezoid(z, bcd, type = "A"),
         bc.100 = integrateTrapezoid(z, bc, type = "A"),
         doc.100 = integrateTrapezoid(z, doc, type = "A"),
         phyc.100 = integrateTrapezoid(z, phyc, type = "A"),
         tdaa.100 = integrateTrapezoid(z, tdaa_c, type = "A")) %>% 
  mutate_at(vars(contains(".100")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".100")), funs(./25)) %>% 
  select(Cruise, Station, CampCN, contains(".100")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>% 
  distinct() %>% 
  ungroup()

integ_300 <- bcd %>% 
  group_by(CampCN) %>% 
  filter(between(z, 100, 300)) %>% 
  mutate(bcd.300 = integrateTrapezoid(z, bcd, type = "A"),
         bc.300 = integrateTrapezoid(z, bc, type = "A"),
         doc.300 = integrateTrapezoid(z, doc, type = "A"),
         phyc.300 = integrateTrapezoid(z, phyc, type = "A"),
         tdaa.300 = integrateTrapezoid(z, tdaa_c, type = "A")) %>% 
  mutate_at(vars(contains(".300")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".300")), funs(./200)) %>% 
  select(Cruise, Station, CampCN, contains(".300")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>% 
  distinct() %>% 
  ungroup()

npp_75 <- npp %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(z <= 75) %>% 
  mutate(npp.75 = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".75")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".75")), funs(./75)) %>% 
  select(Cruise, Station, Date, contains(".75"))  %>%
  distinct() %>% 
  ungroup()

npp_100 <- npp %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(between(z, 75, 100)) %>% 
  mutate(npp.100 = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".100")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".100")), funs(./25)) %>% 
  select(Cruise, Station, Date, contains(".100"))  %>%
  distinct() %>% 
  ungroup()

npp_300 <- npp %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(between(z, 100, 300)) %>% 
  mutate(npp.300 = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".300")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".300")), funs(./200)) %>% 
  select(Cruise, Station, Date, contains(".300"))  %>%
  distinct() %>% 
  ungroup()



integrated.df <- left_join(bcd, integ_75) %>% 
  left_join(., integ_100) %>% 
  left_join(., integ_300) %>% 
  left_join(., npp_75) %>% 
  left_join(., npp_100) %>% 
  left_join(., npp_300)
```

# Save Data

``` r
saveRDS(integrated.df, "~/GITHUB/naames_multiday/Output/processed_data.rds")
```
