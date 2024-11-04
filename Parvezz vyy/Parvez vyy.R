library(readxl)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)


data <- read_excel("Ataher.xlsx", sheet = "species1" )


data %>% 
  #filter( Station == "Station1" ) %>% 
  ggplot(aes(x = Length, y = Frequency, color = as.factor(Station )))+
  geom_line(lty = 4, lwd = 0.3)+
  scale_color_discrete(name = "Station")+
  ylim(0, 4) 

# Assuming you have a data frame fish_data with columns Length, Frequency, and Station
# Example data

# Example frequency data
library(dplyr)
library(tidyr)

# Example frequency data with station column
frequency_data <- data.frame(
  Station = c("A", "A", "B", "B", "C", "C"),
  Length = c(20, 25, 30, 35, 40, 45),
  Frequency = c(5, 8, 10, 6, 3, 7)
)

# Expand the data by repeating each length value according to its frequency
normal_data <- frequency_data %>%
  group_by(Station) %>%
  mutate(normal_Length = list(unlist(mapply(rep, Length, Frequency)))) %>%
  unnest(normal_Length)

# View the resulting normal data
View(normal_data)









 data %>% 
  group_by(Station, Length) %>% 
  count(Frequency) %>% 
  ungroup() %>% 
  select(1,2,4) %>% 
   #filter(Station == "Station1") %>% 
  #View()
  ggplot(aes(x = Length, y = n, linetype = Station, shape = Station, color = Station))+
   geom_point(size = 3)+
  geom_line()+
   ylim(0,5)+
   labs(shape = "Station")
  
 

 
 
 
 # Sample data
 # Sample data
 df <- data.frame(
   Station = rep(c("Station A", "Station B", "Station C"), each = 100),  # Example station labels
   Length = runif(300, min = 10, max = 50),  # Example length data
   Frequency = rpois(300, lambda = 20)  # Example frequency data
 )
 
 
 df <- data
 # Define length intervals
 length_intervals <- seq(10, 36, by = 2)
 
 
 # Create new dataframe to store summarized data
 summary_df <- data.frame(Station = character(), Length = numeric(), Frequency = numeric())

 
 
 # Iterate over each station
 for (station in unique(df$Station)) {
   # Subset data for the current station
   station_data <- df[df$Station == station, ]
   
   # Cut length data into intervals and summarize frequency
   freq_summary <- tapply(station_data$Frequency, cut(station_data$Length, breaks = length_intervals), sum)
   
   # Ensure equal number of rows by padding with NA values
   max_length_intervals <- length(length_intervals)
   n_freq_summary <- length(freq_summary)
   if (n_freq_summary < max_length_intervals) {
     freq_summary <- c(freq_summary, rep(NA, max_length_intervals - n_freq_summary))
   }
   
   # Create dataframe for the summarized data
   station_summary <- data.frame(
     Station = rep(station, max_length_intervals),  # Repeat station label
     Length = length_intervals,  # Length intervals
     Frequency = freq_summary  # Summarized frequency
   )
   
   # Combine with summary dataframe
   summary_df <- rbind(summary_df, station_summary)
 }
 
 
 mdata <- na.omit(summary_df)
 
 myplot <- mdata %>% 
  ggplot(aes(x = Length, y = Frequency, linetype = Station, color = Station, shape = Station))+
  geom_line(size = 1)+
  geom_point(size = 3)+
  theme_bw()+
  theme(legend.position = c(0.15, 0.85))+
  labs( y = "Frequency (n)", x = "Length (cm)", linetype = "Stations",
        shape = "Stations", color = "Stations")
 
 
ggsave("lineplot.pdf", myplot, dpi = 500, height = 4.5, width = 8)
  



df <- read_excel("Ataher.xlsx", sheet = "species2")


length_intervals <- seq(10, 36, by = 2)


# Create new dataframe to store summarized data
summary_df <- data.frame(Station = character(), Length = numeric(), Frequency = numeric())



# Iterate over each station
for (station in unique(df$Station)) {
  # Subset data for the current station
  station_data <- df[df$Station == station, ]
  
  # Cut length data into intervals and summarize frequency
  freq_summary <- tapply(station_data$Frequency, cut(station_data$Length, breaks = length_intervals), sum)
  
  # Ensure equal number of rows by padding with NA values
  max_length_intervals <- length(length_intervals)
  n_freq_summary <- length(freq_summary)
  if (n_freq_summary < max_length_intervals) {
    freq_summary <- c(freq_summary, rep(NA, max_length_intervals - n_freq_summary))
  }
  
  # Create dataframe for the summarized data
  station_summary <- data.frame(
    Station = rep(station, max_length_intervals),  # Repeat station label
    Length = length_intervals,  # Length intervals
    Frequency = freq_summary  # Summarized frequency
  )
  
  # Combine with summary dataframe
  summary_df <- rbind(summary_df, station_summary)
}


max(df$Length)
min(df$Length)

mdata <- na.omit(summary_df)

myplot <- mdata %>% 
  ggplot(aes(x = Length, y = Frequency, linetype = Station, color = Station, shape = Station))+
  geom_line(size = 1)+
  geom_point(size = 3)+
  theme_bw()+
  theme(legend.position = c(0.15, 0.85))+
  labs( y = "Frequency (n)", x = "Length (cm)", linetype = "Stations",
        shape = "Stations", color = "Stations")


ggsave("lineplot2.pdf", myplot, dpi = 500, height = 4.5, width = 8)














 
