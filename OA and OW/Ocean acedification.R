install.packages("oce")
install.packages("marelac", 
                 "oceano",
                 "seacarb",
                 "oceanoGraphy",
                 "gsw")


library(oce)
data(ctd)
plot(ctd, which(c(1,2), type = "l", span = 150))

View(ctd)

?oce

library(ggplot2)
library(readxl)
library(tidyverse)
library(dplyr)

setwd()
data <- read_excel("ocean acidification.xlsx", sheet = "Sheet1" )

data %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "value") %>% 
  filter(parameters == "SST") %>% 
  #ggplot(aes(factor(Year), value, fill = Range))+
  ggplot()+
  geom_line(aes(x = factor(Year)))


data %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "value") %>% 
  filter(parameters == "SST") %>% 
  ggplot(aes(factor(Year), value, fill = "#69b2a8" ))+
  geom_boxplot()


data %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "value") %>% 
  filter(parameters == "SST") %>% 
  ggplot(aes(factor(Year), value, fill = Range))+
  geom_col(position = "dodge" )+
lines(predict(values), col = "black", lwd = 2)


#all in a column chart
#data$Year_last_two <- str_sub(data$Year, -2)

data %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "value") %>%
  #filter(parameters == "SST") %>% 
  ggplot(aes(factor(Year), value, fill = parameters))+
  geom_col()+
  facet_wrap(~parameters, scales = "free_y" )+
  scale_x_discrete(breaks = seq(1998, 2023, by = 3))+
  theme(legend.position = "none")

#with max value
data %>% 
  filter(Range == "Maximum") %>%
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  ggplot(aes(factor(Year), values, fill = parameters))+
  geom_col(position = "dodge")+
  facet_wrap(~parameters, scales = "free_y" )+
  theme(legend.position = "none")+
  scale_x_discrete(breaks = seq(2000, 2023, by = 3))
  
#with min
data %>% 
  filter(Range == "Minimum") %>%
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  ggplot(aes(factor(Year), values, fill = parameters))+
  geom_col(position = "dodge")+
  facet_wrap(~parameters, scales = "free_y" )+
  theme(legend.position = "none")+
  scale_x_discrete(breaks = seq(2000, 2023, by = 3))

#line graph

y1 <- data %>% 
  filter(Range == "Maximum") %>%
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values")

y2 <- data %>% 
  filter(Range == "Minimum") %>%
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values")
 
ggplot()+
  geom_line(data = y1, aes(x = Year, y = values,  color = parameters))+
  facet_wrap(~parameters, scales = "free_y" )+
  theme(legend.position = "none")

ggplot()+
  geom_line(data = y2, aes(x = Year, y = values,  color = parameters))+
  facet_wrap(~parameters, scales = "free_y" )+
  theme(legend.position = "none")


data %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  filter(Range == "Minimum" & parameters == "SST" ) %>% 
  ggplot(aes(x = Year, y = values))+
  geom_line()

data %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  filter(parameters == "SST" ) %>% 
  ggplot(aes(x = Year, y = values, col = Range))+
  geom_line()


#all line plot
data %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  ggplot(aes(x = Year, y = values, color = Range ))+
  facet_wrap(~parameters, scales = "free_y")+
  geom_line()




#line graph of all mean max and min
library(readxl)
library(ggplot2)
library(tidyverse)
library(GGally)

data2 <- read_excel("ocean acidification.xlsx", sheet = "Sheet2" )

#all line plot
allline <- data2 %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  ggplot(aes(x = Year, y = values, color = Range ))+
  facet_wrap(~parameters, scales = "free_y")+
  geom_point(size = 1, alpha = 0.5 )+
  geom_line()+
  theme_bw()+
  theme(axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black")
        )+
  labs(x = "Year", y = "Value", color = "Legend" )
  
ggsave("AllLine.png", allline, dpi = 1500, height = 6.8, width = 9)
ggsave("AllLine.pdf", allline, dpi = 1500, height = 6.8, width = 9)

#boxplot all

data2 %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  ggplot(aes(x = factor(Year), y = values, fill = parameters))+
  facet_wrap(~parameters, scales = "free_y")+
  geom_boxplot(show.legend = F,
               outlier.shape = NA)
  #geom_jitter(show.legend = F)
  


#lm line plot 
data2 %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  ggplot(aes(x = Year, y = values, color = Range))+
  facet_wrap(~parameters, scales = "free_y")+
  geom_point(size = 1.5, alpha = 0.5 )+
  geom_smooth(method = "lm", se = F )+
  theme(legend.position = "none")


#lm line plot 
alllm <- data2 %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  ggplot(aes(x = Year, y = values, color = Range))+
  facet_wrap(~parameters, scales = "free_y")+
  geom_point(size = 1.5, alpha = 0.5 )+
  geom_smooth(se = T, show.legend = F )+
  theme(axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"))+
  labs( y = "Value", color = "Legend")
       


ggsave("AllLM.png", alllm, dpi = 1500, height = 6.8, width = 9)
ggsave("AllLM.pdf", alllm, dpi = 1500, height = 6.8, width = 9)













#############################################################
library(readxl)
library(ggplot2)
library(tidyverse)
library(GGally)

data2 <- read_excel("ocean acidification.xlsx", sheet = "Sheet2" )

data2 %>% 
  pivot_longer(cols = 3:15,
               names_to = "parameters",
               values_to = "values") %>% 
  View()


data <- data.frame(
  x = rep(1:10, 30),
  y = c(rnorm(10, mean = 2), rnorm(10, mean = 5), rnorm(10, mean = 8)),
  group = rep(c("A", "B", "C"), each = 10),
  color = rep(c("red", "blue", "green"), each = 10)
)

data %>% 
  View()

# Fit the linear models for each group
lm_model <- lapply(unique(data$group), function(grp) {
  subset_data <- subset(data, group == grp)
  lm(y ~ x, data = subset_data)
})

# Extract the coefficients from each linear model
coefficients <- lapply(lm_model, coef)

# Create the equation text for each group
equation_text <- lapply(coefficients, function(coef) {
  paste0("y = ", round(coef[1], 2), " + ", round(coef[2], 2), "x")
})

# Create the plot with color, fill, and facet_wrap
ggplot(data, aes(x = x, y = y,  fill = color)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ group) +
  geom_text(data = data.frame(group = unique(group), equation_text), 
            aes(label = equation_text), x = Inf, y = -Inf, 
            hjust = 1, vjust = 0, size = 4)

# Create the plot with color, fill, and facet_wrap







  
library(corrplot)

data2 %>% 
  select(3:15) %>% 
  cor() %>% 
  corrplot(method =  "number", type = "full")

p_values <- data2 %>% 
  select(3:15) %>% 
  cor.mtest()

cm <- data2 %>% 
  select(3:15) %>% 
  cor() %>% 
  corrplot(method = "ellipse", p.mat = p_values)
  

corrplot(cm, method = "number", p.mat = p_values, sig.level = 0.05,
         type = "lower", insig = "blank", addCoef.col = "black")


data2 %>% 
  select(3:15) %>% 
  cor() %>% 
  corrplot(method =  "circle", type = "full", t1.col = "black", addCoef.col = "black")

data2 %>% 
  filter(Range == "Maximum") %>% 
  select(3:15) %>% 
  cor() %>% 
  corrplot(method =  "circle", type = "full", t1.col = "black", addCoef.col = "black")


#ggpairs 

cc <- data2 %>% 
  select(2:15)

pairs(cc)
ggpairs(cc, mapping = aes(color = factor(Year)))



#temperature and salinity diagram (T-S diagram)
library(shape)
library(readxl)
library(ggplot2)
library(tidyverse)
library(GGally)
library(marelac)
library(shape)
library(plot3D)

data <- read_excel("ocean acidification.xlsx", sheet = "Sheet2" )

ts <- data %>% 
  select(SST, SSS)

mint <- min(data$SST)
maxt <- max(data$SST)
mins <- min(data$SSS)
maxs <- max(data$SSS)

salC <- seq(from = mins, to = maxs, length.out = 156 )
tempC <- seq(from = mint, to = maxt, length.out = 156)
sigma.c <- outer(salC, tempC, FUN = function(S, t)sw_dens(S = S, t = t) - 1000 )
sigma.c

#png(filename = 'ts_diagram.png', width = 1, res = 500, pointsize = 12, bg = "white", dpi = 300 )
par(mar = c(5,5,4,6))

contour2D(x = salC, y = tempC, z = sigma.c, lwd = 2, main = "General T-S diagram",
col = "black", xlab = expression("Salinity (ppt)"), ylab = expression("Temperature(\u00B0C)"))
temp <- unlist(ts['SST'], use.names = FALSE)
sal <- unlist(ts['SSS'], use.names = FALSE)
sigma_theta <- sw_dens(S = sal, t = temp) - 1000
scatter2D(sal, temp, colvar = sigma_theta, pch = 16, cex = 1.5, add = TRUE,
          clim = range(sigma.c),colkey = FALSE)
colkey(clim = range(sigma.c),dist = 0.005, side= 4,add=TRUE,
       clab = expression("Density(Kg/m^3)"), col.clab = 'black',
       side.clab = 4,line.clab = 2.5,length = 1, width = 0.8,
       col.axis = 'black',col.ticks = 'black',cex.axis = 0.9)

#dev.off()

# another try with another data
ts <- data %>% 
  select(`PO4=`, `NH4+`)

mint <- min(data$`PO4=`)
maxt <- max(data$`PO4=`)
mins <- min(data$`NH4+`)
maxs <- max(data$`NH4+`)

salC <- seq(from = mins, to = maxs, length.out = 156 )
tempC <- seq(from = mint, to = maxt, length.out = 156)
sigma.c <- outer(salC, tempC, FUN = function(S, t)sw_dens(S = S, t = t) - 1000 )
sigma.c

#png(filename = 'ts_diagram.png', width = 1, res = 500, pointsize = 12, bg = "white", dpi = 300 )
par(mar = c(5,5,4,6))

contour2D(x = salC, y = tempC, z = sigma.c, lwd = 2, main = "General T-S diagram",
          col = "black", xlab = expression("Salinity (ppt)"), ylab = expression("Temperature(\u00B0C)"))
temp <- unlist(ts['PO4'], use.names = FALSE)
sal <- unlist(ts['NH4+'], use.names = FALSE)
sigma_theta <- sw_dens(S = sal, t = temp) - 1000
scatter2D(sal, temp, colvar = sigma_theta, pch = 16, cex = 1.5, add = TRUE,
          clim = range(sigma.c),colkey = FALSE)
colkey(clim = range(sigma.c),dist = 0.005, side= 4,add=TRUE,
       clab = expression("Density(Kg/m^3)"), col.clab = 'black',
       side.clab = 4,line.clab = 2.5,length = 1, width = 0.8,
       col.axis = 'black',col.ticks = 'black',cex.axis = 0.9)






# Heatmap  ----------------------------------------------------------------

#1 cluster
#2 heatmap


data <- read_excel("ocean acidification.xlsx", sheet = "Sheet2" )

hm <- data %>% 
  filter(Range == "Maximum") %>% 
  select(3:15) %>% 
  heatmap()
glimpse(hm)

hm <- as.matrix(hm)

heatmap(hm, cexCol = 0.9, scale = "column", col = cm.colors(365))

#We can aviod dendrograma and see the matrix using colv adn rowv arguments

heatmap(hm, cexCol = 0.9, scale = "column", Colv = NA ) #column dendrogram removed
heatmap(hm, cexCol = 0.9, scale = "column", Rowv =  NA ) #row dendrogram removed
heatmap(hm, cexCol = 0.9, scale = "column", Colv = NA, Rowv = NA ) 



datai <- as.matrix(iris[, -5])
glimpse(datai)

heatmap(datai, cexCol = 0.9, scale = "column")





library("PerformanceAnalytics")

cor <- clw %>% 
  select(2:9)

data <- read_excel("ocean acidification.xlsx", sheet = "Sheet2" )

cor <- data %>% 
  filter(Range == "Mean") %>% 
  select(3:15)

chart.Correlation(cor, histogram = T)






