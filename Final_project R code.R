library(sf)       # For management of sf data
library(dplyr)    # For data manipulation
library(spdep)    # Functions for creating spatial weight, spatial analysis, probability mapping
library(DCluster) # Package with functions for spatial cluster analysis
library(tmap)     # For mapping
library(ggplot2)  # For making various plots
## The package below are used for bivariate maps (The tutorial I referred used)
library(viridis)  # Viridis color scale
library(cowplot)  # stack ggplots
library(tidyr)    # Manipulate string data

# Load Attribute Data-----------
health.raw <- read.csv("Data/HIV-2018_Raw.csv")
summary(health.raw)
edu.raw <- read.csv("Data/Education-2018_Raw.csv")
summary(edu.raw)
tw <- read.csv("Data/TW_HK_MC.csv")

# Manipulate data----------
# Create EDUCATION
edu <- edu.raw %>% 
  mutate(EDUCATION = round((高中+大学专科+大学本科+研究生)/合计 * 100,2))
# Clean health dataset and add new variables
health <- health.raw %>% 
  mutate(CASE = round(POP * (RATE/10))) %>% 
  left_join(edu[, c("PROVINCE", "EDUCATION")], by = "PROVINCE") 
# Clean tw dataset
tw$RATE <- round(tw$CASE * 10 / tw$POP,2)
# Add back to health
health[33:35,c("PROVINCE", "ADM1_PCODE", "POP", "RATE", "CASE")] <- tw[,c("PROVINCE", "ADM1_PCODE", "POP", "RATE", "CASE")]

# Summary and histograms (results are hided)
summary(health)
# hist(health$RATE)  # Right-skewed
# hist(health$EDUCATION)  # Almost normal

# Load Geography Data--------
province <- st_read("Data/chn_adm_ocha_2020_shp/chn_admbnda_adm1_ocha_2020.shp") %>% 
  dplyr::select(ADM1_EN, ADM1_ZH, ADM1_PCODE) # remove variables that will not be used
# Summary and check data (results are hided)
summary(province)
class(province)
st_crs(province)

# Combine Attribute Data with Geodata
hiv <- province %>% left_join(health, by = "ADM1_PCODE") 

summary(health)

# I checked HEALTHCARE here, because, the national level was also included in health dataset, which for this variable is the total number of all the healthcare facilities in China
# hist(hiv$HEALTHCARE)  # Normal? 

rate0 <- health %>% 
  dplyr::filter(RATE < 10) %>%  # Exclude provinces with rate higher than 10
  dplyr::filter(!ADM1_PCODE == "CN") %>%  # The first line is national level, so we exclude it here
  dplyr::select(RATE, POP, CASE)
summary(rate0$RATE)

# Descriptive map by just mapping the incidence rate
tm_shape(hiv) +
  tm_fill('RATE',
          style="fixed",
          breaks = c(-Inf, 1.84, 2.30, 3.06, 6, +Inf),
          palette = 'Oranges',
          alpha = 0.8,
          title = 'Raw Rate (per 100,000 py)',
          legend.hist = T)+ # include a tiny histgram here to see if data are grouped reasonably
  tm_borders(lwd = 0.8) +
  tm_credits('Source: China Health Statistical Yearbook, 2019',
             size = 0.7,
             position = c('LEFT','BOTTOM')) + 
  tm_layout(main.title = 'Incidence Rate of AIDS by Province in China, 2018',
            main.title.size = 1.2,
            inner.margins = c(.05,0.05,0.025,0.025),
            legend.title.size = 0.9 ,
            bg.color = 'lavender',
            legend.outside = T) +
  tm_compass(position = c('LEFT','TOP'), 
             text.size = 0.75,
             color.dark = "grey35")

# Map for EDUCATION
tm_shape(hiv) +
  tm_fill('EDUCATION',
          style="fisher",
          palette = 'Blues',
          alpha = 0.85,
          title = '% of Population',
          legend.hist = T)+
  tm_borders(lwd = 0.8) +
  tm_credits('Source: China Statistical Yearbook, 2019',
             size = 0.7,
             position = c('LEFT','BOTTOM')) + 
  tm_layout(main.title = 'The Percentage of Population with High School Degree or above by Province in China, 2018',
            main.title.size = 1,
            inner.margins = c(.05,0.05,0.025,0.025),
            legend.title.size = 1.2,
            bg.color = 'linen',
            legend.hist.size = 1.4,
            legend.outside = T) +
  tm_compass(position = c('LEFT','TOP'), 
             text.size = 0.75,
             color.dark = "grey35")

# Map for HEALTHCARE
tm_shape(hiv) +
  tm_fill('HEALTHCARE',
          style="fisher",
          palette = 'Greens',
          alpha = 0.9,
          title = 'Counts',
          legend.hist = T)+
  tm_borders(lwd = 0.8) +
  tm_credits('Source: China Health Statistical Yearbook, 2019',
             size = 0.7,
             position = c('LEFT','BOTTOM')) + 
  tm_layout(main.title = 'Number of Healthcare Institutions by Province in China, 2018',
            main.title.size = 1.2,
            inner.margins = c(.05,0.05,0.025,0.025),
            legend.title.size = 1.3,
            bg.color = 'papayawhip',
            legend.hist.size = 1.2,
            legend.outside = T) +
  tm_compass(position = c('LEFT','TOP'), 
             text.size = 0.75,
             color.dark = "grey35")

# Calculate the expected count
# Calculate the nation-wide average rate 
avg_rate <- sum(hiv$CASE) / sum(hiv$POP) * 10

# Calculate SMR and expected case count
hiv <- hiv %>% 
  mutate(smr = round(RATE / avg_rate, 2),
         expected = (avg_rate / 10) * POP)
# hist(hiv$smr) # Check the distribution

# Test for homogeneity
achisq.test(CASE ~ offset(log(expected)), 
            data = hiv,
            model = 'poisson', # because AIDS is a rare disease, we use poisson distribution here
            R = 499)
# Chi-square = 61939.44 , p-value = 0.002

# Poisson probability map
x <- probmap(n = hiv$CASE, x = hiv$POP, 
             alternative = 'greater')
hiv$pmap <- x$pmap

# Combine the p-value result with SMR data
# Group significant regions together
pv <- hiv %>%
  mutate(pmap.pv = ifelse(smr > 1 & pmap < 0.05, 1, 0)) %>%
  group_by(pmap.pv) %>%
  summarise() %>%
  filter(pmap.pv == 1)

# SMR map with probability, using raw rate
tm_shape(hiv) +
  tm_fill('smr',
          style="fisher",
          palette = '-RdYlBu',
          alpha = 0.8,
          title = 'Std. Rate Ratio')+
  tm_borders(lwd = 0.8) +
  tm_credits('Provinces with higher than expected risk (p<0.05) highlighted with red borders',
             size = 2,
             position = c('LEFT','BOTTOM')) + 
  tm_layout(main.title = 'Rate Ratio of AIDS by Province in China, 2018',
            main.title.size = 1.1,
            inner.margins = c(.05,0.05,0.025,0.025),
            legend.text.size = 0.8,
            legend.title.size = 1.2,
            bg.color = 'lavender') +
  tm_shape(pv) +
  tm_borders(lwd = 1.5, col = 'indianred4') +
  # Add chi-square test result to the plot, make the plot more informative
  tm_credits("Chi-square = 61939.44 , p-value = 0.002",
             size = 1.2,
             col = "rosybrown4",
             position = c('LEFT','TOP')) 

# Define neighbors
qnb <- poly2nb(hiv)
# Because Hainan and Taiwan are islands, they do not share boundaries with other provinces, 
# we decide to manually add neighbors for them according to administrative division 
qnb[[34]] <- as.integer(30) # Hainan: Guangdong
qnb[[27]] <- as.integer(31) # Taiwan: Fujian
q_listw <- nb2listw(qnb, style = 'W')

# Because these rate data are not normally distributed, we choose the randomisation = T option
moran.test(hiv$RATE, 
           listw = q_listw, 
           randomisation = T)
# Moran's I: 0.367, p-value = 0.00025 < 0.05
# There is significant positive spatial autocorrelation

# Use linear regression to account for the difference in population-at-risk
# Create weights for each observation
wts <- hiv$POP / sum(hiv$POP) * 34
# Create a linear regression model
reg1 <- lm(RATE ~ 1, 
           data = hiv,
           weights = wts)

lm <- localmoran.exact(reg1, nb = qnb)
lm <- print(lm)
hiv$pvalExact <- lm[,3]

# create lagged local raw_rate - in other words the average of the queen neighbors value
hiv$lag <- lag.listw(q_listw, var = hiv$RATE)

# Create a new dataset that includes standardized values, and then creates a new variable
# 'lm_quad' which takes on the above categorical values.
hiv_lm <- hiv %>%
  mutate(raw_std = as.numeric(scale(RATE)), # scale means standardize to mean 0, 1 SD
         lag_std = as.numeric(scale(lag)),
         lm_quad = factor(case_when(  # All of this is assigning labels based on values
           raw_std >= 0 & lag_std >= 0 & pvalExact < 0.05 ~ 'High-High',
           raw_std <= 0 & lag_std <= 0 & pvalExact < 0.05 ~ 'Low-Low',
           raw_std <= 0 & lag_std >= 0 & pvalExact < 0.05 ~ 'Low-High',
           raw_std >= 0 & lag_std <= 0 & pvalExact < 0.05 ~ 'High-Low',
           pvalExact >= 0.05 ~ 'Non-significant'),
           levels = c('High-High','Low-Low','Low-High','High-Low','Non-significant')))

# Plot
tm_shape(hiv_lm) +
  tm_fill('lm_quad',
          style = 'cat',
          palette = c("salmon3", "lightskyblue3", "palegreen4", "mediumpurple", "lightyellow"),
          alpha = 0.8,
          title = 'Cluster category') +
  tm_borders(lwd = 0.9) +
  tm_layout(main.title = 'LISA Cluster Map for Incidence Rate of AIDS in China, 2018',
            main.title.size = 1.1,
            inner.margins = c(.05,0.05,0.025,0.025),
            legend.text.size = 0.8,
            legend.title.size = 1.2,
            bg.color = 'lavender') +
  tm_credits("Moran's I: 0.367, p-value = 0.00025",
          size = 1.2,
          col = "rosybrown4",
          position = c('LEFT','TOP')) 

## Bivariate map 1-----------
# Create the Bivariate Color Scale-------
# Create 3 buckets for rate
quantiles_rate0 <- hiv %>%
  pull(RATE) %>%
  quantile(probs = seq(0, 1, length.out = 4))

quantiles_rate <- c(0.96, 2.58, 5.60, 17.47) # Manually set the breaks to distinguish outliers of rates

# Create 3 buckets for EDUCATION
quantiles_edu <- hiv %>%
  pull(EDUCATION) %>%
  quantile(probs = seq(0, 1, length.out = 4),na.rm = T)

# Create color scale that encodes two variables
# red for RATE and blue for EDUCATION
bivariate_color_scale <- tibble(
  "3 - 3" = "#3F2949", # high rate, high edu
  "2 - 3" = "#435786",
  "1 - 3" = "#4885C1", # low rate, high edu
  "3 - 2" = "#77324C",
  "2 - 2" = "#806A8A", # medium rate, medium edu
  "1 - 2" = "#89A1C8",
  "3 - 1" = "#AE3A4E", # high rate, low edu
  "2 - 1" = "#BC7C8F",
  "1 - 1" = "#CABED0" # low rate, low edu
) %>%
  gather("group", "fill")

# Join Color Codes to the Data---------
# cut into groups defined above and join RATE
hiv1 <- hiv %>% 
  mutate(
    rate_quantiles = cut(RATE,breaks = quantiles_rate,include.lowest = TRUE),
    edu_quantiles = cut(EDUCATION,breaks = quantiles_edu, include.lowest = TRUE),
    # by pasting the factors together as numbers we match the groups defined
    # in the tibble bivariate_color_scale
    group = paste(as.numeric(rate_quantiles), "-",as.numeric(edu_quantiles))) %>%
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale, by = "group")

# Draw the Map------------
# Give hex value to missing variables
hiv1[is.na(hiv1$fill), "fill"] <- "#CCCCCC"
# Mapping
map <- ggplot(hiv1) +
  # use the "alpha hack" (as the "fill" aesthetic is already taken)
  scale_alpha(name = "",
              range = c(0.6, 0),
              guide = F) + # suppress legend
  # color municipalities according to their RATE / EDUCATION combination
  geom_sf(aes(fill = fill),
    # use thin white stroke for provinces
    color = "white",size = 0.2) +
  # as the sf object municipality_prod_geo has a column with name "fill" that
  # contains the literal color as hex code for each municipality, we can use
  # scale_fill_identity here
  scale_fill_identity() +
  # add titles
  labs(x = NULL,
       y = NULL,
       title = "Map of AIDS Incidence Rate and % of Population with Higher Education",
       subtitle = "China, 2018") +
  # add the theme
  theme_map() +
  theme(plot.title = element_text(color = "black", face = "bold", 
                                  size = 12, vjust = 1),
        plot.subtitle = element_text(color = "black",  
                                     size = 12, vjust = 1))

# Draw the Legend--------
# separate the groups
color_scale <- bivariate_color_scale %>% 
  separate(group, into = c("RATE", "EDUCATION"), sep = " - ") %>%
  mutate(RATE = as.factor(RATE),
         EDUCATION = as.factor(EDUCATION))
color_scale <- color_scale %>% 
  mutate(RATE = factor(RATE, 1:3, c("0.96 to 2.58", "2.58 to 5.60", "5.60 to 17.47")),
         EDUCATION= factor(EDUCATION, 1:3, c("13.4% to 24.3%", "24.3% to 29.0%", "29.0% to 62.5%")))

# Draw legend
legend <- ggplot() +
  geom_tile(
    data = color_scale,
    mapping = aes(x = RATE, y = EDUCATION, fill = fill)) +
  scale_fill_identity() +
  labs(x = "Incidence Rate",
       y = "Education %")+
  theme_minimal(base_size=6.5) +
  theme(
    axis.title = element_text(size = 7),
    axis.text.x = element_text(vjust = 0.7, hjust = 0.7, angle = 45)
  ) +
  # quadratic tiles
  coord_fixed()

# Add Missing Legend------------
legend_miss <- ggplot() +
  geom_tile(data = color_scale[1,],
    mapping = aes(x = RATE, y = EDUCATION, fill = "#CCCCCC")) +
  scale_fill_identity() +
  labs( y = "   ",
    x = "Missing")+
  theme_map() +
  theme(
    axis.title = element_text(size = 7)) 

# Combine Map and Legends
ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(legend, 0.0005, 0.001, 0.25, 0.25) +
  draw_plot(legend_miss, 0.20, 0.081, 0.0647, 0.1074)

# create 3 buckets for HEALTHCARE
quantiles_health <- hiv %>%
  pull(HEALTHCARE) %>%
  quantile(probs = seq(0, 1, length.out = 4),na.rm = T)

# create color scale that encodes two variables
# blue for RATE and green for HEALTHCARE
bivariate_color_scale2 <- tibble(
  "3 - 3" = "#2A5A5B", # high rate, high healthcare
  "2 - 3" = "#5A9178",
  "1 - 3" = "#73AE80", # low rate, high healthcare
  "3 - 2" = "#567994",
  "2 - 2" = "#90B2B3", # medium rate, medium healthcare
  "1 - 2" = "#B8D6BE",
  "3 - 1" = "#6C83B5", # high rate, low healthcare
  "2 - 1" = "#B5C0DA",
  "1 - 1" = "#C5DDDE" # low rate, low healthcare
) %>%
  gather("group", "fill")

# Join Color Codes to the Data----------------
# cut into groups defined above and join RATE
hiv2 <- hiv %>% 
  mutate(
    rate_quantiles = cut(RATE,breaks = quantiles_rate,include.lowest = TRUE),
    health_quantiles = cut(HEALTHCARE,breaks = quantiles_health, include.lowest = TRUE),
    # by pasting the factors together as numbers we match the groups defined
    # in the tibble bivariate_color_scale
    group = paste(as.numeric(rate_quantiles), "-",as.numeric(health_quantiles))) %>%
  # we now join the actual hex values per "group"
  # so each municipality knows its hex value based on the his gini and avg
  # income value
  left_join(bivariate_color_scale2, by = "group")

# Draw the Map------------
# Give hex value to missing variables
hiv2[is.na(hiv2$fill), "fill"] <- "#CCCCCC"
# Mapping
map2 <- ggplot(hiv2) +
  # use the "alpha hack" (as the "fill" aesthetic is already taken)
  scale_alpha(name = "",
              range = c(0.6, 0),
              guide = F) + # suppress legend
  # color municipalities according to their RATE / HEALTHCARE combination
  geom_sf(aes(fill = fill),
          # use thin white stroke for provinces
          color = "white",size = 0.2) +
  # as the sf object municipality_prod_geo has a column with name "fill" that
  # contains the literal color as hex code for each municipality, we can use
  # scale_fill_identity here
  scale_fill_identity() +
  # add titles
  labs(x = NULL,
       y = NULL,
       title = "Map of AIDS Incidence Rate and # of Healthcare Institution",
       subtitle = "China, 2018") +
  # add the theme
  theme_map() +
  theme(plot.title = element_text(color = "black", face = "bold", 
                                  size = 12, vjust = 1),
        plot.subtitle = element_text(color = "black",  
                                     size = 12, vjust = 1))

# Draw the Legend-----------
# separate the groups
color_scale2 <- bivariate_color_scale2 %>% 
  separate(group, into = c("RATE", "HEALTHCARE"), sep = " - ") %>%
  mutate(RATE = as.factor(RATE),
         HEALTHCARE = as.factor(HEALTHCARE))
color_scale2 <- color_scale2 %>% 
  mutate(RATE = factor(RATE, 1:3, c("0.96 to 2.58", "2.58 to 5.60", "5.60 to 17.47")),
         HEALTHCARE= factor(HEALTHCARE, 1:3, c("4,450 to 22,691", "22,691 to 35,300", "35,300 to 85,088")))

# Draw legend
legend2 <- ggplot() +
  geom_tile(
    data = color_scale2,
    mapping = aes(x = RATE, y = HEALTHCARE, fill = fill)) +
  scale_fill_identity() +
  labs(x = "Incidence Rate",
       y = "Healthcare Institution #")+
  theme_minimal(base_size = 6.5) +
  theme(
    axis.title = element_text(size = 7),
    axis.text.x = element_text(vjust = 0.7, hjust = 0.7, angle = 45)
  ) +
  # quadratic tiles
  coord_fixed()

# Combine Map and Legend
ggdraw() +
  draw_plot(map2, 0, 0, 1, 1) +
  draw_plot(legend2, 0.0005, 0.001, 0.25, 0.25) +
  draw_plot(legend_miss, 0.20, 0.081, 0.0647, 0.1074)

# Test the collinearity of covariates
overall_case_rate <- sum(hiv$CASE) / sum(hiv$POP)
HIV <- hiv
st_geometry(HIV) <- NULL
vars <- c("HEALTHCARE", "EDUCATION")
pairs(HIV[,vars]) 
# There is no linear relationship between HEALTHCARE and EDUCATION

# Delete Hong kong, Macau and Taiwan in this regression analysis part for missing HEALTHCARE and EDUCATION information
HIV <- hiv[-c(11,14,27),] 

# Create weights indicating the relative contribution of each tract to total population
HIV$wts <- HIV$POP/sum(HIV$POP,na.rm = T)

# Fitting crude/unconditional model
m0 <- lm(RATE ~ 1,
         data = HIV,
         weights = wts)
summary(m0)

m1 <- lm(RATE ~ HEALTHCARE+EDUCATION,
         data = HIV,
         weights = wts)
summary(m1) 
# Both HEALTHCARE and EDUCATION are insignificant

# package MASS for calculating studentized residuals
HIV$m0_resids <- MASS::studres(m0)
HIV$m1_resids <- MASS::studres(m1)

# Add HIV dataset back to geodata to plot the regions with missing value
# 祖国统一万岁！(耶!!)
HIV_join <- HIV %>% 
  dplyr::select(ADM1_PCODE, m0_resids, m1_resids)
st_geometry(HIV_join) <- NULL
map_residual <- province %>% 
  left_join(HIV_join, by="ADM1_PCODE")

# style = 'sd' option in tm_fill(). In contrast to style = 'quantile', 
# this will emphasize extreme outliers (high or low), which are locations
# where the model is particularly poor-fitting

# the residual maps for crude model and conditional model adjusted by HEALTHCARE and EDUCATION
tm_shape(map_residual) + 
  tm_fill(c('m0_resids', 'm1_resids'),
          style = 'sd',
          palette = 'RdYlGn',
          title = 
            c('Residuals on\nUnconditional\nMean Model','Residuals on\nConditional\nMean Model')) +
  tm_borders() +
  tm_layout(legend.position = c('RIGHT','bottom'),
            inner.margins = c(.05,0.05,0.025,0.025),
            legend.format = list(digits = 1),
            legend.title.size = 1 ,
            bg.color = 'lavender',
            main.title = 'Residuals of Average AIDS Rate by Province in China, 2018\nControlled for Education Level and Healthcare Resources',
            main.title.size = 1)

# Creating spatial neighbors
qnb <- poly2nb(HIV)
qnb[[31]] <- as.integer(27) # assign Hainan province, an island in China, 
# to be neighbor of Guangdong province considering close distance and similar demographic characteristics

q_listw <- nb2listw(qnb, style = 'W') # row-standardized weights

# lm.morantest() function
lm.morantest(m0, listw = q_listw)
lm.morantest(m1, listw = q_listw)
# the p-values in both of the moran tests are < 0.05
# indicating there are significant clusterings for the residuals both in crude model and conditional model
