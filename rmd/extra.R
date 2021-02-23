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


#####Z Heat Relabund #####

# Custom Z-score Table (Top 10 taxa)

#calculate z-scores (value - mean / sd)

value.relabund <- relabund %>% 
  select(Date, z, `Proteobacteria _ Alphaproteobacteria _ SAR11_clade`:`Proteobacteria _ Gammaproteobacteria _ Salinisphaerales`)

mean.relabund <- relabund  %>% 
  select(Date, z, `Proteobacteria _ Alphaproteobacteria _ SAR11_clade`:`Proteobacteria _ Gammaproteobacteria _ Salinisphaerales`) %>%
  group_by(z) %>%
  summarise_at(vars(`Proteobacteria _ Alphaproteobacteria _ SAR11_clade`:`Proteobacteria _ Gammaproteobacteria _ Salinisphaerales`), mean, na.rm = T) %>% 
  ungroup() %>%
  setNames(paste0(names(.), ".ave" )) %>% 
  mutate(z = z.ave) %>% 
  select(z, everything(), -z.ave)

sd.relabund <- relabund  %>% 
  select(Date, z, `Proteobacteria _ Alphaproteobacteria _ SAR11_clade`:`Proteobacteria _ Gammaproteobacteria _ Salinisphaerales`) %>%
  group_by(z) %>% 
  summarise_at(vars(`Proteobacteria _ Alphaproteobacteria _ SAR11_clade`:`Proteobacteria _ Gammaproteobacteria _ Salinisphaerales`), sd, na.rm = T) %>% 
  ungroup() %>% 
  setNames(paste0(names(.), ".sd" )) %>% 
  mutate(z = z.sd) %>% 
  select(z, everything(), -z.sd)



zscore <- left_join(value.relabund, mean.relabund, by = "z") %>% 
  pivot_longer(cols = -c(Date,z)) %>% 
  arrange(Date, z, name) %>% 
  mutate(name = gsub(".ave", "", name)) %>% 
  group_by(Date, z, name) %>% 
  mutate(diff = dplyr::first(value) - dplyr::last(value)) %>% 
  select(-value) %>% 
  distinct() %>% 
  pivot_wider(names_from = name, values_from = diff) %>% 
  ungroup() %>% 
  left_join(., sd.relabund, by = "z") %>% 
  pivot_longer(cols = -c(Date,z)) %>% 
  arrange(Date, z, name) %>% 
  mutate(name = gsub(".sd", "", name)) %>% 
  group_by(Date, z, name) %>% 
  mutate(zscore = dplyr::first(value)/dplyr::last(value)) %>% 
  select(-value) %>% 
  distinct() %>% 
  pivot_wider(names_from = name, values_from = zscore) 

# ASV Depth Profiles

sar11 <- relabund %>% 
  ggplot(aes(x = z, y = `Proteobacteria _ Alphaproteobacteria _ SAR11_clade` * 100, group = interaction(Date, CampCN))) +
  geom_line(aes(linetype = eddy, color = Date), size = 0.7) +
  geom_point(aes(shape = eddy, fill = Date), size = 6, color = "black", stroke = 1, alpha = 0.7) + 
  labs(x = expression(italic(paste(""))), y = expression(italic(paste(""))), colour = "", fill = "", linetype = "Eddy Location", shape = "Eddy Location") +
  scale_x_reverse(breaks = pretty_breaks(), expand = c(0,0)) +
  coord_flip() +
  scale_y_continuous(expand = c(0,0)) +
  scale_shape_manual(values = c("Core" = 21, "Outside" = 24)) +
  scale_fill_date(low = "#0db5e6", high = "#d31f2a") +
  scale_color_date(low = "#0db5e6", high = "#d31f2a") +
  guides(fill = guide_colourbar(barheight = 20, barwidth = 2, frame.colour = "black", frame.linewidth = 2,ticks.colour = "black", ticks.linewidth = 1), color = F) +
  custom_theme() +
  theme(panel.spacing.x = unit(1, "cm"),
        axis.text.x = element_text(angle = 0),
        legend.key.size = unit(1, "cm")) +
  ggtitle("SAR11")

# Constrained Ordination

# Above we used unconstrained ordinations (NMDS) to show relationships between samples in low dimensions. We can use a constrained ordination to see how environmental variables are associated with these changes in community composition. We constrain the ordination axes to linear combinations of environmental variables. We then plot the environmental scores onto the ordination

# Remove data points with missing metadata
ps_not_na <- n2s4 %>%
  subset_samples(
    !is.na(mld) &
      !is.na(npp) & 
      !is.na(doc) &
      !is.na(ba) &
      !is.na(n) & 
      !is.na(aou) & 
      !is.na(beamt) & 
      !is.na(bcd) & 
      !is.na(o2) & 
      !is.na(phyc) &
      !is.na(days) & 
      !is.na(fl) & 
      !is.na(dens) & 
      !is.na(sal) & 
      !is.na(temp) &
      eddy == "Core"
  )

bray <- phyloseq::distance(ps_not_na, method = "bray")

# CAP ordinate
cap_ord <- ordinate(ps_not_na, method = "CAP", distance = bray, formula = ~ mld + doc + n  + bcd + ba  + fl  + z + npp +  aou + beamt + phyc + days
)






# CAP plot
cap.plot <- plot_ordination(ps_not_na, ordination = cap_ord,  axes = c(1,2)) + 
  stat_ellipse(aes(color = days, group = interaction(days), linetype = eddy), type = "t") +
  geom_point(aes(fill = days, shape = z_interv), alpha = 0.6, stroke = 2, size = 4) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_gradientn(colors = odv.colors) +
  scale_color_gradientn(colors = odv.colors) +
  theme_classic2(base_size = 16) +
  guides(fill = guide_colorbar(title = "Days"), shape = guide_legend(title = "Depth Interval"), linetype = F, color = F) +
  ggtitle("NAAMES 2 Station 4 Eddy Core")
#removing one of the plotting layers (there are points within points)

cap.plot$layers <- cap.plot$layers[-1]



# Now add the environmental variables as arrows
arrowmat <- vegan::scores(cap_ord, display = "bp")

# Add labels, make a data.frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)

# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

# Make a new graphic
cap.graph <- cap.plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "gray", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

# Do a permutational ANOVA on constrained axes used in ordination

anova(cap_ord)

