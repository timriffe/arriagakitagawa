
library(tidyverse)
library(coddecomp)
library(ggridges)
library(scales)
library(ggpubr)
library(mgcv)
# ----------------------------------------------------------------------- #
# data from Sergi
load("data/data_esp_1621.RData")

# TR: these data prep/naming steps were repeated throughout, so i took care of them once 
# at the start. And the name is indicative of what is in the file.
data_5_prepped <-
  data_esp_1621  |> 
  as_tibble() |> 
  # redundant category
  filter(year  != "Total") |> 
  rename(age = age5) |> 
  # turn age and year to numeric 
  mutate(age = parse_number(as.character(age)),
         year = as.integer(year),
         # group years into 2 periods pre and post covid
         year = ifelse(between(year, 2016, 2019), "2016-2019", "2020-2021"),
         year = as.factor(year),
         # combine secondary education into 1 
         educ = ifelse(str_detect(educ, "Secundaria"), "Secundaria", educ),
         educ = case_when(educ == "Primaria" ~ "Primary",
                          educ == "Secundaria" ~ "Secondary",
                          educ == "Superior" ~ "Higher",
                          TRUE ~ "Total"),
         causa = str_remove_all(causa, "^[:digit:]+. ")) |>
  group_by(sex, educ, year, cause = causa, age) |> 
  # add up data to match new periods and education groups
  summarise(deaths  = sum(deaths),
            pop     = sum(pop),
            .groups = "drop") 

# ----------------------------------------------------------------------- #
# I will use ths data further for predict() with GAM model results
new_data <- expand_grid(
  age  = seq(35, 100, 1),
  cause = unique(data_esp_1621$causa),
  year  = c("2016-2019", "2020-2021"),
  educ  = unique(data_esp_1621$educ),
  sex   = unique(data_esp_1621$sex),
  pop   = 1) |>
  # these categories are redundant
  filter(year  != "Total") |>
  # for model
  mutate(year = as.factor(year),
         educ = ifelse(str_detect(educ, "Secundaria"), "Secundaria", educ),
         educ = case_when(educ == "Primaria" ~ "Primary",
                          educ == "Secundaria" ~ "Secondary",
                          educ == "Superior" ~ "Higher",
                          TRUE ~ "Total"),
         cause = str_remove_all(cause, "^[:digit:]+. ")) |>
  distinct()
# ----------------------------------------------------------------------- #
# Cleaning and rearranging original data original data 
# TR: Rather than relying on nesting and mapping and stuff, 
# I would have made a function that takes df in df out, which
# would simplify the pipeline greatly.
mxc_single <- data_5_prepped |>
  # nest by sex education groups and causes 
  # unfortunately all together do not work, I tried
  # the confidence intervals in this case are too narrow
  # and variability is very low
  # so we fit separately for sex, educ group and cause
  group_nest(sex, educ, cause) |>
  # add data for predict()
  nest_join(new_data, by = join_by(cause, sex, educ)) |>
  # GAM model with log pop as offset and p-spline on age
  mutate(model = map(
    data,
    ~ gam(
      deaths ~ s(age, bs = "ps") + year + offset(log(pop)),
      data   = .,
      family = quasipoisson
    )
  )) |>
  # here predict() happens
  mutate(pred = map2(.x = model, 
                     .y = new_data, ~ .y %>% # this one is necessary
                       mutate(new = predict(.x, 
                                            newdata = .,
                                            type    = "response"))
  )) |>
  # remove intermediate columns and unnest
  dplyr::select(-c(data:model)) |>
  unnest(pred) |> 
  dplyr::select(-pop) |> 
  rename(mx = new)

mxc_single |> 
  filter(cause == "All", sex != "Total", year == "2016-2019") |> 
  ggplot(aes(x = age, y = mx, color = educ, linetype = sex)) +
  geom_line() +
  scale_y_log10()
