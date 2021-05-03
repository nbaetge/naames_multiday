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

    ## Error in get(genname, envir = envir) : object 'testthat_print' not found

We’ll use a CCF of 5 from N2S4 surface and deep experiments to convert
bact abund to bact
c

# Import Data

``` r
bf <- read_rds("~/GITHUB/naames_multiday/Input/master/processed_bf.8.2020.rds") %>%  
  select(Cruise:degree_bin, CampCN, Z_MLD, EZD, Target_Z, interp_DOC, DOC_sd, interp_TDAA, interp_tdaa_c, TDAA_sd, Asp:Lys, interp_O2_Winkler, O2_Winkler_sd, interp_N_N, N_N_sd,  interp_Chl_a_Fluor,interp_Pro_Influx:interp_Nano_Influx, interp_BactProd, BactProd_sd, interp_BactProd_C, BactProd_C_sd, interp_BactAbund, BactAbund_sd) %>% 
  group_by(CampCN) %>% 
  mutate(surf_bala = ifelse(Target_Z == 5, Bala, NA),
         surf_gaba = ifelse(Target_Z == 5, GABA, NA),
         surf_Val = ifelse(Target_Z == 5, Val, NA),
         surf_Lys = ifelse(Target_Z == 5, Lys, NA),
         surf_Ser = ifelse(Target_Z == 5, Ser, NA),
         surf_Leu = ifelse(Target_Z == 5, Leu, NA),
         surf_Arg = ifelse(Target_Z == 5, Arg, NA),
         surf_Tyr = ifelse(Target_Z == 5, Tyr, NA)) %>% 
  fill(c(surf_bala:surf_Tyr), .direction = "updown") %>% 
  ungroup() %>%
  mutate(Bala = ifelse(Target_Z == 0, surf_bala, Bala),
         GABA = ifelse(Target_Z == 0, surf_gaba, GABA),
         Val = ifelse(Target_Z == 0, surf_Val, Val),
         Lys = ifelse(Target_Z == 0, surf_Lys, Lys),
         Ser = ifelse(Target_Z == 0, surf_Ser, Ser),
         Leu = ifelse(Target_Z == 0, surf_Leu, Leu),
         Arg = ifelse(Target_Z == 0, surf_Arg, Arg),
         Tyr = ifelse(Target_Z == 0, surf_Tyr, Tyr)) %>% 
  select(-c(surf_bala:surf_Tyr)) %>% 
  rename(bp = interp_BactProd_C, 
         sd_bp = BactProd_C_sd, 
         leu_incorp = interp_BactProd,
         sd_leu_incorp = BactProd_sd,
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
         rel.nano = nano/phyto,
         rel.bala = Bala/tdaa,
         rel.gaba = GABA/tdaa,
         rel.balagaba = (Bala + GABA)/tdaa,
         rel.fresh.amino = (Val + Lys + Ser + Arg + Leu + Tyr)/tdaa) %>% 
  mutate_all(~ replace(., is.nan(.), NA)) %>% 
  mutate_at(vars(Station), as.numeric) %>% 
  left_join(., read_csv("~/GITHUB/naames_multiday/Input/master/phytoC.csv") %>% 
              select(Cruise, Station, CampCN, Target_Z, PhytoC, sd_PhytoC)) %>% 
  #convert phytoC from ug C / L to µmol C / m3
  mutate(PhytoC = PhytoC * (10^3/12),
         sd_PhytoC = sd_PhytoC * (10^3/12)) %>% 
  group_by(Cruise, Station) %>% 
  mutate(ez = mean(EZD, na.rm = T),
         sd_ez = sd(EZD, na.rm = T)
         ) %>% 
  ungroup() 

saveRDS(bf, "~/GITHUB/naames_multiday/Input/bottle_data.rds")

ctd <- read_rds("~/GITHUB/naames_multiday/Input/master/deriv_naames_ctd.rds") %>% 
    rename(lat = "Latitude [degrees_north]",
         z = bin_depth) %>% 
  mutate(bin = round(lat, 1),
         Date = ymd_hm(`yyyy-mm-ddThh:mm:ss.sss`),
         Date = as_date(round_date(Date, unit = "day")),
         aou_c = deriv_aou_umol_l * 0.72) %>% 
  rename(aou = deriv_aou_umol_l)

#lol, measured surface par but not par from ctd
saveRDS(ctd, "~/GITHUB/naames_multiday/Input/ctd_data.rds")

bge <- read_rds("~/GITHUB/naames_multiday/Input/master/processed_BCD.rds") %>% 
  select(Season, Station, station_bge, cruise_bge) %>% 
  filter(Season %in% c("Late Spring", "Early Autumn")) %>% 
  distinct()
saveRDS(bge, "~/GITHUB/naames_multiday/Input/bge_data.rds") 

#units for npp are  mg C / m3 / d...conv to to mmol C / m3 / d 
npp <- read_rds("~/GITHUB/naames_multiday/Input/master/Z_resolved_model_NPP.rds") %>% 
  rename(z = depth,
         npp = NPP) %>% 
  mutate(npp = (npp)/12)
saveRDS(npp, "~/GITHUB/naames_multiday/Input/npp_data.rds")
```

## note on C-OPS PAR

not enough PAR data to estimate the depth corresponding to the 0.1%
light level for NAAMES 2 Station 4….for instance, on May 24, surf par
was on average 1977.154:

\-0.1% \* 1977.154 = 1.977154 -PAR at the deepest depth of the C-OPS
cast was 4.582724 at 93 meters

so in the absence of the necessary amount of data, i’ll call upper mz
100-200
m

``` r
par_s4_24 <- read_csv("~/GITHUB/naames_multiday/Input/master/PAR/NAAMES2_20160524_142911_C-OPS_PAR.csv") %>% 
  mutate(bin = round(Depth)) %>% 
  group_by(bin) %>% 
  mutate(ave_par = mean(PAR)) %>% 
  ungroup() %>% 
  mutate(surf_par = ifelse(bin == 0, ave_par, NA)) %>% 
  fill(surf_par) %>% 
  mutate(perc1 = 0.01*surf_par,
         perc01 = 0.001*surf_par) %>% 
  add_row() %>% 
  fill(surf_par:perc01) %>% 
  mutate(PAR = ifelse(is.na(PAR), perc1, PAR)) %>% 
  add_row() %>% 
  fill(surf_par:perc01) %>% 
  mutate(PAR = ifelse(is.na(PAR), perc01, PAR)) %>% 
  arrange(-PAR) %>% 
  mutate(bin = ifelse(is.na(bin), na.approx(bin), bin)) %>% 
  filter(PAR == perc1 | PAR == perc01)
```

    ## Parsed with column specification:
    ## cols(
    ##   Depth = col_double(),
    ##   PAR = col_double()
    ## )

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

# Fill 0m data for CTD parmaters (AOU)

``` r
#split the df by CampCN  
ctd.list <- split(ctd, ctd$CampCN)

add.func <- function(morty){
  morty[nrow(morty) + 1,] <- NA
  morty$z[is.na(morty$z)] <- 0
  rick <- morty %>% 
    fill(., Cruise:id, .direction = c("updown")) %>% 
    arrange(CampCN, z) %>% 
    fill(., pres_db:Date, .direction = c("updown")) %>% 
    fill(., aou, .direction = "up") 
   }

#apply function to list 
added.list <- lapply(ctd.list, add.func)

#save the list as a data frame 
ctd.df <- plyr::ldply(added.list, data.frame) %>% 
  select(-.id)  
```

# Estimate BCD and BC

We’ll estimate BCD for each cast of each station of all the cruises
using bacterial production and the station averages of BGE. Otherwise,
they are calculated using the cruise average bge. We’ll also convert
bacterial abundance to bacterial carbon for N2S4 using station average
of CCFs (from POC taken from initial population).

BCD and specific growth rates will be integrated for the 0-100 m, the
100-200 m, and the 200-300 m boxes.

``` r
#calculate bcd and specific growth rates
bcd <- interpolated.df %>% 
  # full_join(.,   bge) %>% 
  mutate(###########BCD##################
         #units are  µmol C / m3 / d 
         #bp in nmol C / L / d is equivalent to  µmol C / m3 / d 
         bcd = ifelse(Cruise == "AT34" & Station == 4, round(bp/0.24), NA),
         ba = ba * 10^3,
         bc = ba * (5/12) / (10^12),
         bc = ifelse(Cruise == "AT34" & Station == 4, bc, NA),
          # bc = ifelse(Cruise == "AT38" & Station == 6, ba * (43/12) / (10^12), bc),
         phyto = phyto * 10^3,
         sd_ba = sd_ba * 10^3) %>% 
  rename(z = Target_Z,
         mld = Z_MLD)


ezd <- bcd %>%
  select(Cruise, Station, Date, ez) %>%
  distinct()
```

# Integrations

``` r
integ_data <- bcd %>% 
  filter(Cruise == "AT34" & Station == 4 | Cruise == "AT38" & Station == 6)

npp_integ <- npp %>%
  filter(Cruise == "AT34" & Station == 4 | Cruise == "AT38" & Station == 6) %>% 
  left_join(., ezd)

aou_integ <- ctd.df %>% 
  filter(Cruise == "AT34" & Station == 4 | Cruise == "AT38" & Station == 6) %>% 
  left_join(., ezd)


integ_ez <- integ_data %>% 
  group_by(CampCN) %>% 
  filter(z <= 75) %>% 
  mutate(bcd.ez = integrateTrapezoid(z, bcd, type = "A"),
         ba.ez = integrateTrapezoid(z, ba, type = "A"),
         bp.ez = integrateTrapezoid(z, bp, type = "A"),
         bc.ez = integrateTrapezoid(z, bc, type = "A"),
         doc.ez = integrateTrapezoid(z, doc, type = "A"),
         chl.ez = integrateTrapezoid(z, chl, type = "A"),
         phyc.ez = integrateTrapezoid(z, phyc, type = "A"),
         phyto.ez = integrateTrapezoid(z, phyto, type = "A"),
         tdaa.ez = integrateTrapezoid(z, tdaa, type = "A"),
         tdaa_c.ez = integrateTrapezoid(z, tdaa_c, type = "A"),
         bala.ez = integrateTrapezoid(z, Bala, type = "A"),
         gaba.ez = integrateTrapezoid(z, GABA, type = "A")) %>% 
  mutate_at(vars(contains(".ez")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".ez")), funs(./75)) %>% 
  select(Cruise, Station, CampCN, contains(".ez")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>%
  distinct() %>% 
  ungroup()


integ_200 <- integ_data %>% 
  group_by(CampCN) %>% 
  filter(z <= 200, z >= 100) %>% 
  mutate(bcd.200 = integrateTrapezoid(z, bcd, type = "A"),
         ba.200 = integrateTrapezoid(z, ba, type = "A"),
         bp.200 = integrateTrapezoid(z, bp, type = "A"),
         bc.200 = integrateTrapezoid(z, bc, type = "A"),
         doc.200 = integrateTrapezoid(z, doc, type = "A"),
         chl.200 = integrateTrapezoid(z, chl, type = "A"),
         phyc.200 = integrateTrapezoid(z, phyc, type = "A"),
         phyto.200 = integrateTrapezoid(z, phyto, type = "A"),
         tdaa.200 = integrateTrapezoid(z, tdaa, type = "A"),
         tdaa_c.200 = integrateTrapezoid(z, tdaa_c, type = "A"),
         bala.200 = integrateTrapezoid(z, Bala, type = "A"),
         gaba.200 = integrateTrapezoid(z, GABA, type = "A")) %>% 
  mutate_at(vars(contains(".200")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".200")), funs(./(100) )) %>% 
  select(Cruise, Station, CampCN, contains(".200")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>% 
  distinct() %>% 
  ungroup()

integ_300 <- integ_data %>% 
  group_by(CampCN) %>% 
  filter(between(z, 200, 300)) %>% 
  mutate(bcd.300 = integrateTrapezoid(z, bcd, type = "A"),
         ba.300 = integrateTrapezoid(z, ba, type = "A"),
         bp.300 = integrateTrapezoid(z, bp, type = "A"),
         bc.300 = integrateTrapezoid(z, bc, type = "A"),
         doc.300 = integrateTrapezoid(z, doc, type = "A"),
         chl.300 = integrateTrapezoid(z, chl, type = "A"),
         phyc.300 = integrateTrapezoid(z, phyc, type = "A"),
         phyto.300 = integrateTrapezoid(z, phyto, type = "A"),
         tdaa.300 = integrateTrapezoid(z, tdaa, type = "A"),
         tdaa_c.300 = integrateTrapezoid(z, tdaa_c, type = "A"),
         bala.300 = integrateTrapezoid(z, Bala, type = "A"),
         gaba.300 = integrateTrapezoid(z, GABA, type = "A")) %>% 
  mutate_at(vars(contains(".300")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".300")), funs(./100)) %>% 
  select(Cruise, Station, CampCN, contains(".300")) %>% 
  # mutate(
  # ######### Estimate BCD:NPP #########
  #        bcd.npp = int.bcd/int.NPP * 100) %>% 
  distinct() %>% 
  ungroup()

npp_ez <- npp_integ %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(z <= 75) %>% 
  mutate(npp.ez = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".ez")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".ez")), funs(./75)) %>% 
  select(Cruise, Station, Date, contains(".ez"))  %>%
  distinct() %>% 
  ungroup()

npp_200 <- npp_integ %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(z <= 200, z >= 100) %>%  
  mutate(npp.200 = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".200")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".200")), funs(./(100))) %>% 
  select(Cruise, Station, Date, contains(".200"))  %>%
  distinct() %>% 
  ungroup()

npp_300 <- npp_integ %>% 
  group_by(Cruise, Station, Date) %>% 
  filter(between(z, 200, 300)) %>% 
  mutate(npp.300 = integrateTrapezoid(z, npp, type = "A")) %>% 
  mutate_at(vars(contains(".300")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".300")), funs(./100)) %>% 
  select(Cruise, Station, Date, contains(".300"))  %>%
  distinct() %>% 
  ungroup()


aou_ez <- aou_integ %>% 
  group_by(Cruise, Station, CampCN) %>% 
  filter(z <= 75) %>% 
  mutate(aou.ez = integrateTrapezoid(z, aou , type = "A")) %>% 
  mutate_at(vars(contains(".ez")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".ez")), funs(./75)) %>% 
  select(Cruise, Station, CampCN, contains(".ez"))  %>%
  distinct() %>% 
  ungroup()

aou_200 <- aou_integ %>% 
  group_by(Cruise, Station, CampCN) %>%  
  filter(z <= 200, z >= 100) %>% 
  mutate(aou.200 = integrateTrapezoid(z, aou , type = "A")) %>% 
  mutate_at(vars(contains(".200")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".200")), funs(./(100))) %>% 
  select(Cruise, Station, CampCN, contains(".200"))  %>%
  distinct() %>% 
  ungroup()

aou_300 <- aou_integ %>% 
  group_by(Cruise, Station, CampCN) %>% 
  filter(between(z, 200, 300)) %>% 
  mutate(aou.300 = integrateTrapezoid(z, aou , type = "A")) %>%  
  mutate_at(vars(contains(".300")), round) %>% 
  # depth normalize 
  mutate_at(vars(contains(".300")), funs(./100)) %>% 
  select(Cruise, Station, CampCN, contains(".300"))  %>%
  distinct() %>% 
  ungroup()



integrated.df <- left_join(bcd, integ_ez) %>% 
  left_join(., integ_200) %>% 
  left_join(., integ_300) %>% 
  left_join(., npp_ez) %>% 
  left_join(., npp_200) %>% 
  left_join(., npp_300) %>% 
  left_join(., aou_ez) %>% 
  left_join(., aou_200) %>% 
  left_join(., aou_300)
```

# Save Data

``` r
saveRDS(integrated.df, "~/GITHUB/naames_multiday/Output/processed_data.rds")
```
