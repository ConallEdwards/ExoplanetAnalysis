library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(tidyverse)
library(psych)


#Load exoplanet dataset into R
exoplanet_df <- read.csv("exoplanets4.csv", stringsAsFactors = TRUE)

#View the first 5 rows of the dataset
head(exoplanet_df)

#View the structure of the dataset
str(exoplanet_df)

#View the summary statistics of the dataset
summary(exoplanet_df)


#Keeping only confirmed exoplanets observed by Kepler
exoplanet_df <- filter(exoplanet_df, default_flag == 1)
exoplanet_df <- filter(exoplanet_df, disc_facility == "Kepler")


#Is there any missing data?
colSums(is.na(exoplanet_df))

#Dropping irrelevant columns and those that have too much missing data
exoplanet_df <- subset(exoplanet_df, 
                       select = -c(pl_name, hostname, pl_masse, pl_insol, default_flag,
                                   disc_facility))




names(exoplanet_df) <- c('num_of_stars', 'num_of_planets', 'discovery_method', 'discovery_year',
                         'orbit_period_days', 'orbit_semi_major_axis',
                         'planet_radius', 'planet_eccentricity', 'planet_equilibrium_temperature', 
                         'stellar_spectral_type', 'stellar_effective_temperature', 'stellar_radius', 
                         'stellar_mass', 'stellar_metallicity', 'stellar_metratio', 'stellar_surface_gravity', 
                         'distance_from_solar_system', 'stellar_gaia_magnitude')



#Imputing missing data 
exoplanet_df$distance_from_solar_system[is.na(exoplanet_df$distance_from_solar_system)] <- mean(exoplanet_df$distance_from_solar_system, na.rm=TRUE)
exoplanet_df$orbit_period_days[is.na(exoplanet_df$orbit_period_days)] <- mean(exoplanet_df$orbit_period_days, na.rm=TRUE)
exoplanet_df$orbit_semi_major_axis[is.na(exoplanet_df$orbit_semi_major_axis)] <- mean(exoplanet_df$orbit_semi_major_axis, na.rm=TRUE)
exoplanet_df$planet_radius[is.na(exoplanet_df$planet_radius)] <- mean(exoplanet_df$planet_radius, na.rm=TRUE)
exoplanet_df$stellar_effective_temperature[is.na(exoplanet_df$stellar_effective_temperature)] <- mean(exoplanet_df$stellar_effective_temperature, na.rm=TRUE)
exoplanet_df$planet_eccentricity[is.na(exoplanet_df$planet_eccentricity)] <- mean(exoplanet_df$planet_eccentricity, na.rm=TRUE)
exoplanet_df$planet_equilibrium_temperature[is.na(exoplanet_df$planet_equilibrium_temperature)] <- mean(exoplanet_df$planet_equilibrium_temperature, na.rm=TRUE)
exoplanet_df$stellar_radius[is.na(exoplanet_df$stellar_radius)] <- medan(exoplanet_df$stellar_radius, na.rm=TRUE)
exoplanet_df$stellar_mass[is.na(exoplanet_df$stellar_mass)] <- mean(exoplanet_df$stellar_mass, na.rm=TRUE)
exoplanet_df$stellar_metallicity[is.na(exoplanet_df$stellar_metallicity)] <- mean(exoplanet_df$stellar_metallicity, na.rm=TRUE)
exoplanet_df$stellar_surface_gravity[is.na(exoplanet_df$stellar_surface_gravity)] <- mean(exoplanet_df$stellar_surface_gravity, na.rm=TRUE)
exoplanet_df$stellar_gaia_magnitude[is.na(exoplanet_df$stellar_gaia_magnitude)] <- mean(exoplanet_df$stellar_gaia_magnitude, na.rm=TRUE)


#Fixing formatting issue in metallicity ratio
exoplanet_df$stellar_metratio <- str_replace(exoplanet_df$stellar_metratio,
                                                      "\\[Fe/H\\[", "\\[Fe/H\\]")
exoplanet_df$stellar_metratio <- str_replace(exoplanet_df$stellar_metratio,
                                             "\\[Me/H\\]", "\\[Fe/H\\]")
exoplanet_df$stellar_metratio <- str_replace(exoplanet_df$stellar_metratio,
                                             "\\[m/H\\]", "\\[M/H\\]")




#Discovery method usage per decade
ggplot(exoplanet_df, aes(y=discovery_method)) + geom_bar() + labs(y = "Discovery Method",
       x = "Number of times method was successfully used",
       title= "A graph showing the methods used to discover exoplanets") +
  facet_grid(discovery_decade~.)



#How far away are these planets from our solar system?
ggplot(exoplanet_df, aes(x=distance_from_solar_system)) + geom_boxplot() +
  labs(x = "Distance (in Parsec)", main = "A boxplot to show the distance the planets are from our solar system")+ 
  facet_grid(discovery_method~discovery_decade)






#Do exoplanet masses have a normal distribution?
ggplot(exoplanet_df, aes(x=Mass)) + geom_boxplot() +
  labs(x = "Mass of Exoplanet (in unit masses of the Earth)",
       title = "The distribution of the mass of exoplanets")

        






#Question 1 - does the Earth have a similar radius to other exoplanets?
summary(exoplanet_df$planet_radius)


#View distribution of the data 
ggplot(exoplanet_df, aes(x=planet_radius)) + 
  geom_histogram(color="darkblue", fill = "lightblue") + 
  labs(x = "Radius of the planet, measured in radius of the Earth",
       y = "The number of exoplanets",
       title = "A histogram to show how exoplanet size is distributed") +
  theme(plot.title = element_text(hjust = 0.5))



wilcox.test(exoplanet_df$planet_radius, 
            mu = 1,
            alternative = "two.sided")











#Question 2 - what factors significantly correlate to the brightness of a stellar?
pairs.panels(exoplanet_df[c("stellar_effective_temperature", "stellar_radius", "stellar_mass", "stellar_metallicity",
                   "num_of_stars", "num_of_planets", "stellar_gaia_magnitude")])





brightness_model <- lm(formula = stellar_gaia_magnitude ~ stellar_effective_temperature + 
                         stellar_radius + stellar_mass + stellar_metallicity +
                         stellar_surface_gravity + distance_from_solar_system +
                         num_of_stars + num_of_planets , data = exoplanet_df )
summary(brightness_model)



exoplanet_df$stellar_density <- exoplanet_df$stellar_mass / ((exoplanet_df$stellar_radius * 2) * 3.1415 * (4/3))



brightness_model2 <- lm(formula = stellar_gaia_magnitude ~ stellar_effective_temperature + 
                         stellar_mass + stellar_radius + stellar_density +
                          stellar_metallicity  +
                         stellar_surface_gravity + distance_from_solar_system +
                         num_of_stars, data = exoplanet_df )
summary(brightness_model2)


boxplot(brightness_model2$residuals, main="A boxplot showing the residuals of the brightness model")












#Question 3 - Is the radius of the Sun similar to that of other stellars? 
summary(exoplanet_df$stellar_radius)

ggplot(exoplanet_df, aes(x=stellar_radius)) + 
  geom_histogram(color="darkblue", fill = "lightblue") + 
  labs(x = "The radius of the stellar, in radii of the sun",
       y = "The number of stellars",
       title = "A histogram to show how stellars radius are distributed") +
  theme(plot.title = element_text(hjust = 0.5))



t.test(x = exoplanet_df$stellar_radius, mu = 1.00)









#Question 4 - Do exoplanets take longer than a year (365.25 days) to orbit around their stellar? 
summary(exoplanet_df$orbit_period_days)

ggplot(exoplanet_df, aes(x=orbit_period_days)) + 
  geom_histogram(color="darkblue", fill = "lightblue") +
  labs(x = "The number of days it takes an exoplanet to orbit its stellar",
       y = "The number of exoplanets",
       title = "A histogram to show how exoplanet orbital periods are distributed")+
  theme(plot.title = element_text(hjust = 0.5))


binom.test(x = count(exoplanet_df$orbit_period_days > 30)[2,2],
           n = 2778,
           alternative = "less")



