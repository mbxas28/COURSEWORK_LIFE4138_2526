#installing and librarying the required packages
install.packages("tidyverse")
library(tidyverse)
library(ggplot2)

#reading the file into R
pumpkin_data <- read.csv("D:/UoN/modules/4138 - Coding/coursework/Pumpkins/pumpkins_datasets/pumpkins_10.csv")

#heaviest pumpkin (variety, location, year)
pumpkin_data %>%
  arrange(desc(weight_lbs)) %>%
  slice(1) %>% #selecting the first answer
  select(variety, city, id)

#creating a new column: weight_kg, that shows weight converted from lbs to kgs
pumpkin_data$weight_kg <- pumpkin_data$weight_lbs * 0.453592

#creating new column: weight_class, that categorises pumpkins by light, medium, heavy
pumpkin_data$weight_class <- ifelse(pumpkin_data$weight_lbs < 942, "light",
                                    ifelse(pumpkin_data$weight_lbs < 1884, "medium",
                                           "heavy"))

#plotting est vs actual weight, coloured by class
ggplot(aes(x = est_weight, y = weight_lbs, colour = weight_class) , data = pumpkin_data)+
  geom_point() + 
  labs(title = "Estimated Weight vs Weight in Pounds") +
  xlab("Estimated Weight in Lbs") +
  ylab("Weight in Lbs")

#filtering for 3 countries 
filtered_pumpkin <- pumpkin_data %>%
  filter(country == "Germany" | country == "France" | country == "Austria")

#exporting/saving file in csv format on laptop
write.csv(filtered_pumpkin, "D:/UoN/modules/4138 - Coding/coursework/assignment work/Pumpkin results/filtered_pumpkin.csv")

#finding mean weight (in lbs) of pumpkins for each country and arranging them in descending order 
filtered_pumpkin %>%
  group_by(country) %>%
  summarise_at(vars(weight_lbs), mean) %>%
  arrange(desc(weight_lbs))

#finding the lowest mean weight of each variety for every country
filtered_pumpkin %>%
  group_by(variety, country) %>%
  summarise_at(vars(weight_lbs), mean) %>%
  arrange(weight_lbs)

#filtering to find mean weight of each variety of pumpkin for each country
filtered_pumpkin %>%
  group_by(country, variety) %>%
  summarise(mean_weight = mean(weight_lbs)) %>%
  arrange(country, mean_weight) 

#making a boxplot for pumpkin weight distributions for each country
ggplot(aes(x = country, y = weight_lbs), data = filtered_pumpkin)+
  geom_boxplot(fill = 'lightblue')+
  labs(title = "Pumpkin Weight Distribution") +
  xlab("Countries") +
  ylab("Weight in Lbs")

#redrawing above plot as a facet plot for varieties
ggplot(aes(x = country, y = weight_lbs, colour = country), data = filtered_pumpkin)+
  geom_boxplot()+
  facet_wrap(~variety) +
  scale_colour_manual(values = c(
    'Austria' = 'turquoise',
    'France' = 'pink',
    'Germany' = 'darkolivegreen'
  )) + 
  labs(title = "Pumpkin Weight Distribution Across Varieties for 3 Countries") +
  xlab("Countries") +
  ylab("Weight in Lbs")

