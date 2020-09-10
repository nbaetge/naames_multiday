#####PLOT WITH CUSTOM TBL####

#we can plot using out custom table

#Here's a crude way to look at specific taxa within a treatment
custom_plot.data <- custom_tbl_merge %>% 
  filter(Treatment %in% c("Control", "NOB", "AOA")) %>% 
  filter(Replicate == "a")

custom_plot <- custom_plot.data %>% 
  ggplot(aes(x= Days, y = custom_plot.data$`Campylobacteria _ Campylobacterales _ Arcobacteraceae _ Arcobacter` )) +
  geom_point(aes(color = Treatment)) +
  facet_wrap( ~ Bottle)

custom_plot


#one thing you can try is combining your cell count data with the relative abundance data...
