###############################################################################
                    #####################################
                    ####Project multiresistance in EU####
                    #####################################
###############################################################################
### Carga de librerias ###
devtools::install_github("rstudio/addinexamples", type = "source")
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
#########################Tidy and Wrangling database###########################

#### carga de base de datos de resistencia a los antibióticos####
ab_bact_db<- read.table('ECDC_surveillance_data_Antimicrobial_resistance.csv', 
                        header=TRUE, sep = ',', stringsAsFactors = FALSE)
str(ab_bact_db)
options(digits=10)

#ModificaciÃ³n de la base de datos a Tibble, eliminaciÃ³n de columnas no necesarias###
ab_bact_db <- ab_bact_db %>% as_tibble(rm.na=TRUE) %>%
    select(-c('HealthTopic', 'Unit', 'TxtValue')) %>%
    separate('Population', c('Bacteria', 'Antibiotic'), sep = '\\|') %>%
    mutate(NumValue = str_replace_all(NumValue, "-", ""))

# Convertir un columna en nÃºmerica_ en este caso NumValue. 
ab_bact_db <- ab_bact_db %>% mutate(NumValue = as.numeric(NumValue)) #%>%
ab_bact_db %>% filter(is.na(NumValue))


################################################################################################
###Organizando la BD por Indicador y antibiÃ³tico ####BASE DE DATOS ORGANIZADA PARA TRABAJO######
################################################################################################

percent <- function(x) {x/100} # formula para convertir % en frecuencia

# cÃ³digo para limpiar la bd, cambiar nombres de los elementos en las filas para que sean mÃ¡s
# cortos y poder organizarlos mejor; se une el indicador con el AB y luego se pasan las filas a 
# columnas 'spread' de la uniÃ³n entre AB e indicador; finalmente se limpia la base a travÃ©s del
# del control de calidad de 'data quality_completeness age.

tidy_by_indic2 <- ab_bact_db %>%
    mutate(Indicator = recode (Indicator, "R - resistant isolates" = "r_isolates", 
                               "R - resistant isolates, percentage  " = "r_percentage", 
                               "Total tested isolates" = "total_isolates", 
                               "I - susceptible, increased exposure isolates" = "i_isolates" , 
                               "S - susceptible isolates" = "s_isolates", 
                               "Completeness age" = "compl_age", 
                               "Completeness gender" = "compl_gender", 
                               "Penicillin non-wild-type isolates, percentage" =
                                   "penicilin_isolate_percentage")) %>% 
    mutate(Antibiotic = recode(Antibiotic, 
                               'Combined resistance (fluoroquinolones, aminoglycosides and carbapenems)' =
                                   'CR_Fluor;Amino;Carba', 
                               'Combined resistance (third-generation cephalosporin, fluoroquinolones and aminoglycoside)' =
                                   'CR_thr;cepha;fluor;amino', 
                               'Combined resistance (at least three of piperac. and tazob., fluoroq., ceftaz., aminogl. and carbapenems)' =
                                   'CR_atlleast_three')) %>%
    unite(AB_indicator, Antibiotic, Indicator) %>%
    spread(AB_indicator, NumValue) %>% 
    mutate_at(vars(matches('percentage')), percent) %>% 
    drop_na(`Data quality_compl_age`)

################################################################################################

##############################################################################################
#####Scraping para obtener regiones y organizacion de las regiones en la base de datos########
##############################################################################################
url <- 'https://www.worldatlas.com/articles/the-four-european-regions-as-defined-by-the-united-nations-geoscheme-for-europe.html'
tab <- read_html(url) %>% html_nodes("table")

ws <- tab %>% html_table(fill = TRUE)

ws <- data.frame(ws) %>% 
    rbind(c('', '','Cyprus', '')) %>% 
    rbind(c('', 'Ireland', '', '')) %>%
    mutate(Eastern.Europe = recode(Eastern.Europe, 
                                   'Czech Republic' = 'Czechia')) %>%
    rename('Eastern' = 'Eastern.Europe', 
           'Northern' = 'Northern.Europe', 
           'Southern' = 'Southern.Europe', 
           'Western' = 'Western.Europe')

ws <- ws %>% gather(Region, Country, c(,1:4)) %>% 
    as_tibble() %>% 
    mutate_all(na_if,"") %>% 
    drop_na()

################################################################################################

db_multiv2 <- tidy_by_indic2 %>% select(!contains('compl')) ### eliminando "data quality" (completness variables)

### verificando información de la base de datos de multiresistencia ###
mul_R <- db_multiv2 %>% group_by(Bacteria) %>% 
    select(contains('CR_')) %>% 
    select(contains('r_percentage'))
view(mul_R)

mul_R_na <- mul_R[rowSums(is.na(mul_R[2:ncol(mul_R)])) 
                  != ncol(mul_R[2:ncol(mul_R)]), ] %>% 
    rename('R1' = names(.[2]), 
           'R2' = names(.[3]), 
           'R3' = names(.[4]))

x <- mul_R_na %>% group_by(Bacteria) %>% count()
y <- db_multiv2 %>% group_by(Bacteria) %>% count() %>% rename('n2' = 'n')
left_join(x, y) %>% mutate(Change = n-n2) # para saber la cantidad de información se quita si se eliminan todos los NAs

mul_R_na1 <-mul_R_na %>%
    gather(R_type, R_multi, starts_with("R"))

mul_R_na1 %>% ggplot(aes(x=R_multi, fill=Bacteria)) +
    geom_histogram(na.rm = TRUE, alpha=0.5)

### filtro por 3 bacterias y organización únicamente multiresistencia ###
multi_R_td <- db_multiv2 %>% group_by(RegionName, Bacteria, RegionCode, Time) %>% 
    filter(Bacteria == "Escherichia coli" |
               Bacteria == "Klebsiella pneumoniae" |
               Bacteria == "Pseudomonas aeruginosa") %>% 
    select(contains('CR_')) %>% select(contains('r_percentage')) %>% 
    rename('Country' = names(.[1]),
           'Code' = names(.[3]), 
           'Year' = names(.[4]), 
           'r1' = names(.[5]), #r1 = CR_atlleast_three_r_isolates
           'r2' = names(.[6]), #r2 = CR_Fluor;Amino;Carba_r_isolates
           'r3' = names(.[7]))#r1 = CR_thr;cepha;fluor;amino_r_isolates

ec <- multi_R_td %>% filter(Bacteria == "Escherichia coli") %>% select(r3) %>% rename('R_multi' = 'r3')
kp <- multi_R_td %>% filter(Bacteria == "Klebsiella pneumoniae") %>% select(r3) %>% rename('R_multi' = 'r3')
pa <- multi_R_td %>% filter(Bacteria == "Pseudomonas aeruginosa") %>% select(r1) %>% rename('R_multi' = 'r1')

multi_r_td1 <- rbind(ec, kp, pa)

multi_r_td2 <- multi_r_td1 %>% 
    left_join(., ws) %>% 
    relocate(Country, .after = Bacteria) %>% 
    relocate(Region, .after = Country)

###indentificación de NAs###
multi_r_td2 %>% group_by(Year) %>% summarise(na_count = sum(is.na(R_multi)))
multi_r_td2 %>% group_by(Year) %>% summarise(na_count = sum(is.na(R_multi)))

###revisión de base de datos, remoción de años y paises ###
multi_r_td2 %>% group_by(Bacteria) %>% count(Year) %>% spread(Bacteria, n)
multi_r_td2 %>% group_by(Country) %>% count(Year) %>% spread(Year, n) %>% view()
multi_r_td2 %>% group_by(Bacteria) %>% count(Country) %>% spread(Bacteria, n) %>% view()

########################################################################################
### se decide dejar la base desde 2006-2018 y se remueve Slovakia país con menos obs ###
########################################################################################
multi_r_td2 %>% nrow()
multi_r_td3 <- multi_r_td2 %>% filter(Year >= 2006 & !Country == 'Slovakia')

multi_r_td3 %>% group_by(Bacteria) %>% count(Year) %>% spread(Bacteria, n)
multi_r_td3 %>% group_by(Country) %>% count(Year) %>% spread(Year, n) %>% view()
multi_r_td3 %>% group_by(Bacteria) %>% count(Country) %>% spread(Bacteria, n) %>% view()



########################################################################
###Interpolación toda la base de datos de resistencia #################
########################################################################

multi_r_int<- multi_r_td3 %>% group_by(Country, Bacteria) %>%
    mutate_at(vars(ncol(multi_r_td3)), funs(na_interpolation(., option = 'spline'))) %>% 
    ungroup()

multi_r_int <- multi_r_int %>% mutate(R_multi = ifelse(R_multi < 0, 0, R_multi))

multi_r_int %>% group_by(Bacteria) %>% count(Year) %>% spread(Bacteria, n)
multi_r_int %>% group_by(Country) %>% count(Year) %>% spread(Year, n) %>% view()
multi_r_int %>% group_by(Bacteria) %>% count(Country) %>% spread(Bacteria, n) %>% view()

multi_r_int %>% write.csv(file = 'multi_r_int.csv', row.names = F) ## creación de archivo .csv con la BD interpolada



############################################
###gráficas para comprobar interpolación####
############################################
summary(multi_r_td3$R_multi)
summary(multi_r_int$R_multi)

int<- multi_r_int %>% ggplot(aes(x=R_multi, fill=Bacteria)) +
    geom_density(na.rm = TRUE, alpha=0.5) + xlab("R_multi_int") + scale_x_continuous(trans = "log2")

non_int <- multi_r_td3 %>% ggplot(aes(x=R_multi, fill=Bacteria)) + 
    geom_density(na.rm = TRUE, alpha=0.5) + xlab("R_multi_non") + scale_x_continuous(trans = "log2")

ggarrange(int, non_int, labels = c("A", "B"),
    common.legend = TRUE, legend = "bottom",  nrow = 2)

############################################

##############################################################################
### base de datos datos económicos y sociales ################################
##############################################################################
per_cap_year <- read.csv2('per_cap_year.csv') %>% set_colnames(c('Country', 'Year', 'per_cap_US'))


### base per capital $ en salud por pa????s ###
GDP_year <- read.csv2('GDP_year.csv') %>% set_colnames(c('Country', 'Year', 'GDP_health'))


### organizando base de datos de DDD_comunitario###
commun_ddd <- read.csv2('community_DDD.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                       '2007','2008', '2009', 
                                                                                       '2010','2011', '2012', 
                                                                                       '2013', '2014','2015', 
                                                                                       '2016', '2017','2018'))
commun_ddd <- commun_ddd %>% gather(key = 'Year', value = 'DDD_sys_commun', 2:15)
commun_ddd$Year <- as.integer(commun_ddd$Year)

### organizando base de datos de Gasto de Bolsillo###
Out_pocket_exp <- read.csv2('OutOP_%_H_Exp.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                           '2007','2008', '2009', 
                                                                                           '2010','2011', '2012', 
                                                                                           '2013', '2014','2015', 
                                                                                           '2016', '2017'))
Out_pocket_exp <- Out_pocket_exp %>% gather(key = 'Year', value = 'Out_pocket_exp', 2:14)
Out_pocket_exp$Year <- as.integer(Out_pocket_exp$Year)

### organizando base de datos de % de la población rural###
rural_pop <- read.csv2('Rural_pop_%.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                    '2007','2008', '2009', 
                                                                                    '2010','2011', '2012', 
                                                                                    '2013', '2014','2015', 
                                                                                    '2016', '2017', '2018'))
rural_pop <- rural_pop %>% gather(key = 'Year', value = 'rural_pop', 2:15)
rural_pop$Year <- as.integer(rural_pop$Year)

### organizando base de datos de PIB total###
GDP_total <- read.csv2('GDP_total.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                  '2007','2008', '2009', 
                                                                                  '2010','2011', '2012', 
                                                                                  '2013', '2014','2015', 
                                                                                  '2016', '2017', '2018'))
GDP_total <- GDP_total %>% gather(key = 'Year', value = 'GDP_total', 2:15)
GDP_total$Year <- as.integer(GDP_total$Year)

### organizando base de datos de población total por pa????s###
total_pop <- read.csv2('pop_total_country.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                          '2007','2008', '2009', 
                                                                                          '2010','2011', '2012', 
                                                                                          '2013', '2014','2015', 
                                                                                          '2016', '2017', '2018'))

total_pop <- total_pop %>% gather(key = 'Year', value = 'total_pop', 2:15)
total_pop$Year <- as.integer(total_pop$Year)

### indice de desarrollo humano (http://hdr.undp.org/en/data#) ###
HDI <- read.csv2('HDI_total.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                            '2007','2008', '2009', 
                                                                            '2010','2011', '2012', 
                                                                            '2013', '2014','2015', 
                                                                            '2016', '2017', '2018'))

HDI <- HDI %>% gather(key = 'Year', value = 'HDR', 2:15)
HDI$Year <- as.integer(HDI$Year)

### organizando base de datos de % de la población rural###

#### indicador control corrupcion banco mundial ####
ctrl_corrup <- read.csv2('control_corruption.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                             '2007','2008', '2009', 
                                                                                             '2010','2011', '2012', 
                                                                                             '2013', '2014','2015', 
                                                                                             '2016', '2017', '2018'))

ctrl_corrup <- ctrl_corrup %>% gather(key = 'Year', value = 'ctrl_corrup', 2:15)
ctrl_corrup$Year <- as.integer(ctrl_corrup$Year)

#### Efectividad gobernanza banco mundial ####
GOV_effect <- read.csv2('GOV_effectiv.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                      '2007','2008', '2009', 
                                                                                      '2010','2011', '2012', 
                                                                                      '2013', '2014','2015', 
                                                                                      '2016', '2017', '2018'))

GOV_effect <- GOV_effect %>% gather(key = 'Year', value = 'GOV_effect', 2:15)
GOV_effect$Year <- as.integer(GOV_effect$Year)

#### Seguimiento de leyes (rule_of_law) banco mundial ####
rule_law <- read.csv2('rule_law.csv', na.strings=c("","NA")) %>% set_colnames(c('Country', '2005', '2006',
                                                                                '2007','2008', '2009', 
                                                                                '2010','2011', '2012', 
                                                                                '2013', '2014','2015', 
                                                                                '2016', '2017', '2018'))

rule_law <- rule_law %>% gather(key = 'Year', value = 'rule_law', 2:15)
rule_law$Year <- as.integer(rule_law$Year)

#### Seguimiento de leyes (rule_of_law) banco mundial ####
income_pc <- read.csv2('income_pc_US.csv', na.strings=c("","NA")) %>% set_colnames(c('Country','2005', '2006',
                                                                                '2007','2008', '2009', 
                                                                                '2010','2011', '2012', 
                                                                                '2013', '2014','2015', 
                                                                                '2016', '2017', '2018'))

income_pc <- income_pc %>% gather(key = 'Year', value = 'income_pc', 2:15)
income_pc$Year <- as.integer(income_pc$Year)

##########################################################################
############ CONSTRUCCIÓN DE BASE DE DATOS TOTAL VARIABLES ###############
##########################################################################
multi_r_var <- left_join(multi_r_int, commun_ddd, by = c('Country','Year')) %>%
    left_join(., per_cap_year, c('Country','Year')) %>%
    left_join(., GDP_total,  c('Country','Year')) %>% 
    left_join(., GDP_year, c('Country','Year')) %>%
    left_join(., Out_pocket_exp,  c('Country','Year')) %>%
    left_join(., income_pc,  c('Country','Year')) %>%
    left_join(., rural_pop, c('Country','Year')) %>%
    left_join(., HDI, c('Country','Year')) %>%
    left_join(., ctrl_corrup, c('Country','Year')) %>%
    left_join(., GOV_effect, c('Country','Year')) %>%
    left_join(., rule_law, c('Country','Year')) %>%
    left_join(., total_pop, c('Country','Year'))


### pasar caracter a numéricos para poder hacer operaciones con a estos datos####
i <- c(6:ncol(multi_r_var))
multi_r_var[ , i] <- apply(multi_r_var[ , i], 2,function(x) as.numeric(as.character(x)))
str(multi_r_var)

######### NAs verificacion ###########
xxx <- split.data.frame(multi_r_var, multi_r_var$Bacteria) 
yyy <- lapply(xxx, function(x){colSums(is.na(x))})
(zzz <- bind_rows(yyy)) # missing values per country
colSums(zzz)

multi_r_var %>% group_by(Bacteria) %>% count(!is.na(.)) %>% spread(Bacteria, n)
multi_r_var %>% group_by(Country) %>% count(Year) %>% spread(Year, n) %>% view()
multi_r_var %>% group_by(Bacteria) %>% count(Country) %>% spread(Bacteria, n) %>% view()

multi_r_var %>% group_by(Bacteria) %>% filter_all(any_vars(is.na(.))) %>% view()

#######################################################
###interpolación base de datos todas las variables#####
#######################################################

multi_var_int <- multi_r_var %>% group_by(Country, Bacteria) %>%
    mutate_at(vars(7:ncol(multi_r_var)), funs(na_interpolation(., option = 'spline'))) %>% 
    ungroup() 

multi_var_int %>% write.csv(file = 'multi_var_int.csv', row.names = F) ## creación de archivo .csv con la BD interpolada

multi_r_var %>% view()

####################################################
###Gráficas descripción de variable Resistencia#####
####################################################

###Density####

DR_Region <- multi_var_int %>% ggplot(aes(R_multi)) + 
    scale_x_continuous(trans = "log2",
                       labels = scales::number_format(accuracy = 0.01)) +
    geom_density(aes(fill=factor(Region)), 
                 alpha=0.5, 
                 position = "stack") + 
    labs(title="Density plot", 
         subtitle="antibiotic resistance total",
         #caption="Source: mpg",
         x="% resistance",
         fill="Region")

DR_Region

###Boxplot####

BR_region <- multi_var_int %>%
    mutate(Region = reorder(Region, R_multi, FUN = mean)) %>%      # reorder
    ggplot(aes(Bacteria, R_multi, fill = Region)) +    # color by continent
    geom_boxplot(alpha=0.5) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("Resistance (log2)") +
    xlab("") + 
    geom_jitter(width = 0.1, alpha = 0.1)+
    scale_y_continuous(trans = 'log2', 
                       labels = scales::number_format(accuracy = 0.01))
BR_country <- multi_var_int %>%
    mutate(Country = reorder(Country, R_multi, FUN = mean)) %>%      # reorder
    ggplot(aes(Country, R_multi, fill = Bacteria)) +    # color by antibiotic
    geom_boxplot(alpha=0.5) +    
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("Resistance (log2)")+
    xlab("") + #geom_point(show.legend = FALSE) + 
    scale_y_continuous(trans = 'log2', 
                       labels = scales::number_format(accuracy = 0.01))

BR_year <- multi_var_int %>%
    mutate(Year = reorder(Year, R_multi, FUN = mean)) %>%      # reorder
    ggplot(aes(Year, R_multi, fill = Bacteria)) +    # color by antibiotic
    geom_boxplot(alpha=0.5) +    
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    ylab("Resistance (log2)")+
    xlab("") + #geom_point(show.legend = FALSE) + 
    scale_y_continuous(trans = 'log2', 
                       labels = scales::number_format(accuracy = 0.01))

BR_region
BR_country
BR_year


###Barplot###
multi_var_int %>%
    mutate(Country = reorder(Country, R_multi, FUN = mean)) %>%      # reorder
    ggplot(aes(Country, R_multi, fill=Region)) +
    geom_bar(stat="summary", position=position_dodge(), alpha=0.7) + 
    coord_flip()+ scale_fill_brewer(palette="Paired") + theme_minimal()

# devtools::install_github("kassambara/ggcorrplot")

# Correlation matrix

corr <- round(cor(multi_var_int[7:ncol(multi_var_int)-1]), 2)
# Plot
ggcorrplot(corr, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method="circle", 
           colors = c("tomato2", "white", "springgreen3"), 
           title="Correlogram of variables", 
           ggtheme=theme_bw)

#### creación de nueva variable nivel de ingresos segun el banco mundial  ####
#### https://blogs.worldbank.org/opendata/new-world-bank-country-classifications-income-level-2020-2021

names(multi_var_int)
multi_var_int %>% 
    mutate(., income_levl = ifelse(income_pc <1036, "Low",
                                                 ifelse(income_pc %in% 1036:4045, "Lower-middle",
                                                        ifelse(income_pc %in% 4046:12535, "Upper-middle", "High")))) %>%
    group_by(income_levl) %>% count()
    
