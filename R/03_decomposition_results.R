# source("R/00_initial_data_preparation.R")
# source("R/01_smoothing_and_ungroupping.R")
tic()
source("R/02_LT_and_average_quantities.R")

mxc_decomp <- mxc_single |>
  filter(sex   != "Total",
         cause != "All",
         educ  != "Total") |>
  pivot_wider(names_from  = sex,
              values_from = mx) |> 
  mutate(mxc_diff = Females - Males) |> 
  full_join(averages, 
            by = join_by(educ, age, year)) |>
  mutate(result = mxc_diff * sensitivity) # 13, sens 12
# ----------------------------------------------------------------------- #
# Kitagawa part: ex by groups
e35_kit <- Lt |>
  filter(sex  != "Total",
         educ !="Total") |> 
  dplyr::select(sex, educ, year, ex) |> 
  pivot_wider(names_from = sex, 
              values_from = ex, 
              names_prefix = "e35_") |> 
  mutate(e35_diff = e35_Females - e35_Males,
         e35_avg = (e35_Females + e35_Males) / 2)
# ----------------------------------------------------------------------- #
# population prevalence data from original source
if (exists("data_5_prepped")){
  prev <- data_5_prepped |> 
    as_tibble() |>
    filter(year  != "Total",
           sex   != "Total",
           educ  !="Total",
           cause == "All",
           age   == 35) |> 
    group_by(sex, educ, year) |> 
    summarise(deaths  = sum(deaths),
              pop     = sum(pop),
              .groups = "drop") |> 
    group_by(sex, educ, year) |> 
    # population by groups
    summarise(pop = sum(pop), .groups = "drop") |> 
    group_by(sex, year) |>
    # scale population
    mutate(prev = pop / sum(pop)) |> 
    ungroup() |> 
    select(-pop) |> 
    group_by(educ, year)  |> 
    pivot_wider(names_from   = sex, 
                values_from  = prev, 
                names_prefix = "st_") 

  prev |> 
    write_csv("data/prev.csv")
} else {
  # otherwise read in the last one saved (also committed)
  prev <- read_csv("data/prev.csv")
}

struct_kit <-
  prev |> 
  mutate(st_diff = st_Females - st_Males,
         st_mean = (st_Females + st_Males) / 2) 
# ----------------------------------------------------------------------- #
# combine for Kitagawa
kit <- left_join(struct_kit, 
                 e35_kit, 
                 by = join_by(educ, year)) |> 
  mutate(e35_component = st_mean * e35_diff,
         st_component  = e35_avg * st_diff)

decomp_total <- kit |> 
  select(educ, year, e35_component) |> 
  right_join(mxc_decomp, by = join_by(educ, year)) |> 
  group_by(educ, year) |> 
  mutate(result_rescaled = result / sum(result) * e35_component)
toc()

write_csv(decomp_total, file = "data/decomp_results.csv")
