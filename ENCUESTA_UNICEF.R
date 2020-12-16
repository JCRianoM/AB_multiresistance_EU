library(dplyr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(rvest)
library(stringr)
library(imputeTS)
library(factoextra)
library(FactoMineR)
library(GGally)
library(ggpubr)
library(ggplot2)
library(ggfortify)
library(ggcorrplot)


h <- read.csv("country_1.csv") %>% set_colnames('Country')

i_d <- read.csv2('income_pc_US.csv', na.strings=c("","NA")) %>% set_colnames(c('Country','2005', '2006',
                                                                               '2007','2008', '2009', 
                                                                               '2010','2011', '2012', 
                                                                               '2013', '2014','2015', 
                                                                               '2016', '2017', 'income_pc'))

d_d <- HDI <- read.csv2('HDI_total.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                   '2007','2008', '2009', 
                                                                                   '2010','2011', '2012', 
                                                                                   '2013', '2014','2015', 
                                                                                   '2016', '2017', 'HDI'))


d_d <- d_d[, c(1, 15)]
i_t <- i_d[, c(1, 15)]

ht <- left_join(h, i_t) %>% left_join(., d_d) %>% as.tibble()

i <- c(2:ncol(ht))
ht[ , i] <- apply(ht[ , i], 2,function(x) as.numeric(as.character(x)))

h_NAs <- ht %>% filter(!is.na(income_pc) & !is.na(HDI))

h_NAs <- h_NAs %>% 
    mutate(., income_levl = ifelse(income_pc <1036, "Low",
                                   ifelse(income_pc %in% 1036:4045, "Lower-middle",
                                          ifelse(income_pc %in% 4046:12535, "Upper-middle", "High"))))

h_NAs %>% group_by(income_levl) %>% count()
h_NAs %>% view()

h_NAs %>% ggplot(aes(HDI)) + 
    scale_x_continuous(trans = "log2",
                       labels = scales::number_format(accuracy = 0.01)) +
    geom_density(aes(fill=factor(income_levl)), 
                 alpha=0.5, 
                 position = "stack") + 
    labs(title="Density plot", 
         subtitle="HDI or Income per capital",
         #caption="Source: mpg",
         x="HDI",
         fill="Region")


