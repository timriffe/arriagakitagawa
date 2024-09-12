# ----------------------------------------------------------------------- #
library(tidyverse)
library(coddecomp)
library(ggridges)
library(xtable)
library(scales)
library(ggpubr)
library(mgcv)
# ----------------------------------------------------------------------- #
# data from Sergi
load("data/data_esp_1621.RData")

# initial data preparation
data_5_prepped <- data_esp_1621  |>
  as_tibble() |>
  # redundant category
  filter(year  != "Total") |>
  rename(age = age5) |>
  # turn age and year to numeric 
  mutate(age  = parse_number(as.character(age)),
         year = as.integer(year),
         # group years into 2 periods pre and post covid
         year = ifelse(between(year, 2016, 2019), "2016-2019", "2020-2021"),
         year = as.factor(year),
         # combine secondary education into 1 
         educ = ifelse(str_detect(educ, "Secundaria"), "Secundaria", educ),
         educ = case_when(educ == "Primaria"   ~ "Primary",
                          educ == "Secundaria" ~ "Secondary",
                          educ == "Superior"   ~ "Higher",
                          TRUE                 ~ "Total"),
         causa = str_remove_all(causa, "^[:digit:]+. "))|>
  group_by(sex, educ, year, cause = causa, age) |> 
  # add up data to match new periods and education groups
  summarise(deaths  = sum(deaths),
            pop     = sum(pop),
            .groups = "drop") |> 
  mutate(cause = ifelse(cause ==  "Other diseases", "Congenital",
                        cause),
         cause = ifelse(cause ==  "Other minor causes", 
                        "Other diseases", cause)
         )
