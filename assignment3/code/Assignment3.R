#Call the library that need to use in this code
library(readr)
library(dplyr)
library(ggplot2)
library(knitr)

#Call the data as tsv via read_tsv function as data
data <- read_tsv("raw_data/PlateletHW.tsv")

#Pipeline the data to take absolute to the ADP column and store in df variable
df <- data %>%
  mutate(ADP_abs = abs(ADP))

#Save new cleaned file to clean_data.tsv name in clean_data folder.
write_tsv(df, "clean_data/clean_data.tsv")

#Call clean_data data frame and take log to ADP column
clean_data <- read_tsv("clean_data/clean_data.tsv")
clean_data$ADP_log <- log(clean_data$ADP)

#Create linear model to each SNPs.
log_liner_A <- lm(ADP_log ~ rs4244285, data =clean_data)
log_liner_B <- lm(ADP_log ~ rs4986893, data =clean_data)
log_liner_C <- lm(ADP_log  ~ rs662, data =clean_data)

#Call summary data to see the Min Max and p-values
summary(log_liner_A)

#Use ggplot to show the linear model. I use green color to show rs4244285.
ggplot(clean_data, aes(x = rs4244285, y = ADP_log)) +
  geom_point() +
  geom_smooth(method = "lm", color = "green") +
  labs(title = "Association between ADP_log and rs4244285",
       x = "rs4244285 Genotype (0, 1, 2)",
       y = "ADP-Induced Platelet Aggregation")+theme_minimal()

#Use qqnorm and qqline to show for check the data that this data are normalize already.
qqnorm(log_liner_A$residuals)
qqline(log_liner_A$residuals, col = "blue")


summary(log_liner_B)


#Use ggplot to show the linear model. I use red color to show rs4986893.
ggplot(clean_data, aes(x = rs4986893 , y = ADP_log)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  labs(title = "Association between ADP and rs4986893",
       x = "rs4986893 Genotype (0, 1, 2)",
       y = "ADP-Induced Platelet Aggregation") +theme_minimal()

qqnorm(log_liner_B$residuals)
qqline(log_liner_B$residuals, col = "violetred1")


summary(log_liner_C)

#Use ggplot to show the linear model. I use yelloq color to show rs662.
ggplot(clean_data, aes(x = rs662 , y = ADP_log)) +
  geom_point() +
  geom_smooth(method = "lm", color = "yellow") + 
  labs(title = "Association between ADP and rs662",
       x = "rs662 Genotype (0, 1, 2)",
       y = "ADP-Induced Platelet Aggregation") +theme_minimal()

qqnorm(log_liner_C$residuals)
qqline(log_liner_C$residuals, col = "brown")