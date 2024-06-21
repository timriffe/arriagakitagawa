# initial data
source("00_initial_data_preparation.R")
# ----------------------------------------------------------------------- #
# I will use this data further for predict() with GAM model results
new_data <- expand_grid(
  age  = seq(35, 100, 1),
  year  = c("2016-2019", "2020-2021"),
  pop   = 1)


# smoothing function
smooth_ungroup <- function(.data, new_data) { 
  
  model <- gam(
    deaths ~ s(age, bs = "ps") + year + offset(log(pop)),
    data   = .data,
    family = quasipoisson)
  
  result <- predict(model, 
                    newdata = new_data, 
                    type = "response") |>
    as_tibble() |>
    mutate(age  = new_data$age,
           year = new_data$year) |>
    rename(mx = value)
  
  return(result)
  
}

# smooth and ungroup separately by sex, education, cause
mxc_single <- data_5_prepped |>
  group_nest(sex, educ, cause) |>
  mutate(data = map(data, ~ 
                      smooth_ungroup(.x, new_data = new_data)
  )) |> 
  unnest(data)



