# ----------------------------------------------------------------------- #
library(tidyverse)
library(coddecomp)
library(ggridges)
library(scales)
library(ggpubr)
library(mgcv)
library(xtable)
# ----------------------------------------------------------------------- #
# data from Sergi
load("data/data_esp_1621.RData")

# TR: these data prep/naming steps were repeated throughout, so i took care of them once 
# at the start. And the name is indicative of what is in the file.
data_5_prepped <- data_esp_1621  |>
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
         educ = case_when(educ == "Primaria"   ~ "Primary",
                          educ == "Secundaria" ~ "Secondary",
                          educ == "Superior"   ~ "Higher",
                          TRUE                 ~ "Total"),
         causa = str_remove_all(causa, "^[:digit:]+. ")) |>
  group_by(sex, educ, year, cause = causa, age) |> 
  # add up data to match new periods and education groups
  summarise(deaths  = sum(deaths),
            pop     = sum(pop),
            .groups = "drop") 

# ----------------------------------------------------------------------- #
# I will use this data further for predict() with GAM model results
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
         educ = case_when(educ == "Primaria"    ~ "Primary",
                          educ == "Secundaria"  ~ "Secondary",
                          educ == "Superior"    ~ "Higher",
                          TRUE                  ~ "Total"),
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
# ----------------------------------------------------------------------- #
# calculate life expectancy
# TR: note, we have DemoTools installed already, so could just call 
# DemoTools::lt_single_mx(), which will handle closeouts more cleanly.
# I modified the dx and Lx formulas below for easier handling. Note also
# you had All among the causes, and that even if you had causes without
# the total you would want to SUM rather than MEAN. Pipeline modified greatly.

Lt <- mxc_single |>
  filter(cause == "All") |> 
  group_by(sex, year, educ) |>
  mutate(age = as.integer(age)) %>%
  group_by(sex, educ, cause, year) %>% 
  summarise(ex = mx_to_e0(mx, age = age), .groups = "drop")

# we compare this later with our weighted-average e35 values just to see.
e35_total_compare <- Lt |>
  filter(educ == "Total", sex != "Total") |> 
  select(sex, year, ex)
# ----------------------------------------------------------------------- #
# Average male and female overall mortality for each age, year and educ
averages <- mxc_single %>%
  filter(cause == "All",
         sex   != "Total",
         educ  != "Total") %>%
  mutate(year = as.character(year),
         age = as.integer(age)) %>% 
  dplyr::select(-cause) %>%
  pivot_wider(names_from = sex,
              values_from = mx) %>% 
  group_nest(educ, year) %>% 
  mutate(data = map(data, ~ 
                      sen_arriaga_sym(mx1 = .x$Males,
                                      mx =  .x$Females,
                                      age = .x$age,
                                      closeout = T))) %>% 
  unnest(data) %>% 
  group_by(educ, year) %>% 
  mutate(age = 35:100) %>% 
  ungroup() %>% 
  rename(sensitivity = data)

take arriaga sym and multiply by age specific delta by education group e0

any subset with age specific mortality and calculate tge e1071::diff rates take sym and then sumproduct of sym and delta the difference
mx1 with lower e0 and mx2 will be females
total non st
total st_bind

# ----------------------------------------------------------------------- #
# cause, age, year, and education specific diff of  Females - Males 
# multiplied by the results from sen_arriaga_instantaneous. This object
# contains pairwise edu-specific decompositions.
mxc_decomp <- mxc_single |>
  filter(sex   != "Total",
         cause != "All",
         educ  != "Total") |>
  pivot_wider(names_from  = sex,
              values_from = mx) |> 
  mutate(mxc_diff = Females - Males) |> 
  full_join(averages, 
            by = join_by(educ, age, year)) |>
  mutate(result = mxc_diff * sensitivity)
# ----------------------------------------------------------------------- #
# kitagawa part: ex by groups
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
# mostly data cleaning similar to step 1.
# TR: I needed to make many changes here; we want prevalence
# at age 35-39, specifically, and causes play no role.
struct_kit <- data_5_prepped |> 
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
  # scale populaton
  mutate(prev = pop / sum(pop)) |> 
  ungroup() |> 
  select(-pop) |> 
  group_by(educ, year)  |> 
  pivot_wider(names_from   = sex, 
              values_from  = prev, 
              names_prefix = "st_") |> 
  mutate(st_diff = st_Females - st_Males,
         st_mean = (st_Females + st_Males) / 2) 
# ----------------------------------------------------------------------- #
# combine for Kitagawa
kit <- left_join(struct_kit, 
                 e35_kit, 
                 by = join_by(educ, year)) |> 
  mutate(e35_component = st_mean * e35_diff,
         st_component  = e35_avg * st_diff)

# sum st. component

# gaps
gaps <- kit |> 
  group_by(year) |> 
  summarize(e35_Males   = sum(st_Males * e35_Males),
            e35_Females = sum(st_Females * e35_Females)) |> 
  mutate(gap = e35_Females - e35_Males)

e35_total_compare |> 
  pivot_wider(names_from = sex, values_from = ex, names_prefix = "e35_") |> 
  mutate(gap = e35_Females - e35_Males)

e35_kit |> 
  select(educ, year, e35_Females, e35_Males) |> 
  mutate(gap = e35_Females - e35_Males)
# ----------------------------------------------------------------------- #
# now take educ-specific arriaga full decompositions and weight them
# according to Kitagawa e35 component
decomp_total <- kit |> 
  select(educ, year, e35_component) |> 
  right_join(mxc_decomp, by = join_by(educ, year)) |> 
  group_by(educ, year) |> 
  mutate(result_rescaled = result / sum(result) * e35_component)
# ----------------------------------------------------------------------- #
# separate figures for 2 periods
before <-
  decomp_total |>
  filter(year == "2016-2019") |>
  mutate(cause = case_when(
    !cause %in% c(
      "External",
      "Circulatory",
      "Digestive",
      "Neoplasms",
      "Respiratory"
    ) ~ "Other",
    TRUE ~ cause
  )) |>
  group_by(educ, year, cause, age) |>
  summarise(valuersc = sum(result_rescaled), 
            value = sum(result),.groups = "drop") |>
  mutate(cause = factor(cause),
         cause = reorder(cause, -valuersc)) |>
  mutate(
    educ = as.factor(educ),
    educ = fct_relevel(educ, "Higher", after = Inf)) |> 
  ggplot(aes(x = age, y = value, fill = educ)) +
  geom_density(stat = "identity",
               alpha = 0.8,
               color = "white") +
  facet_grid(vars(cause),vars(educ)) +
  theme_bw() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black")) +
  labs(fill = "Education groups")
# ----------------------------------------------------------------------- #
after <-
  decomp_total |>
  filter(year == "2016-2019") |>
  mutate(cause = case_when(
    !cause %in% c(
      "External",
      "Circulatory",
      "Digestive",
      "Neoplasms",
      "Respiratory"
    ) ~ "Other",
    TRUE ~ cause
  )) |>
  group_by(educ, year, cause, age) |>
  summarise(valuersc = sum(result_rescaled), 
            value = sum(result),.groups = "drop") |>
  mutate(cause = factor(cause),
         cause = reorder(cause, -valuersc)) |>
  mutate(
    educ = as.factor(educ),
    educ = fct_relevel(educ, "Higher", after = Inf)) |> 
  ggplot(aes(x = age, y = valuersc, fill = educ)) +
  geom_density(stat = "identity",
               alpha = 0.8,
               color = "white") +
  facet_grid(vars(cause),vars(educ)) +
  theme_bw() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black")) +
  labs(fill = "Education groups")

st_gaps <- kit |> 
  group_by(year) |> 
  summarize(cc_str = sum(st_component))

tot_gps <- gaps %>% 
  filter(year == "2016-2019") %>% 
  dplyr::select(year, Females = e35_Females, Males = e35_Males) %>% 
  pivot_longer(-year,
               names_to = "sex",
               values_to = "e35") %>% 
  mutate(educ = "Total")

ed_gps <- e35_kit |> 
  select(educ, year, Females = e35_Females, Males = e35_Males)  |>      
  pivot_longer(Females:Males, names_to = "sex", values_to = "e35") |> 
  filter(year == "2016-2019")

tot_gps <- gaps %>% 
  filter(year == "2016-2019") %>% 
  dplyr::select(year, Females = e35_Females, Males = e35_Males) %>% 
  pivot_longer(-year,
               names_to = "sex",
               values_to = "e35") %>% 
  mutate(educ = "Total")


mort_gaps <- decomp_total |> 
  group_by(year) |> 
  summarize(cc_mort = sum(result_rescaled))


dec_gaps <- full_join(mort_gaps, st_gaps, by = join_by(year)) |> 
  mutate(dec_gap = cc_str + cc_mort)

ed_gps %>% 
  bind_rows(tot_gps) %>%
  mutate(educ = as.factor(educ),
         educ = relevel(educ, "Total"),
         educ = relevel(educ, "Primary"),
         educ = relevel(educ, "Secondary"),
         educ = relevel(educ, "Higher")) |> 
  ggplot(aes(y = educ, x = e35, fill = sex)) +
  geom_col(position = position_dodge(), color = "white") +
  theme_bw() + 
  coord_flip()+
  scale_x_continuous(breaks = pretty_breaks(n = 12)) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", color = "black",size = 14),
    axis.text.y = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_blank(),
    axis.text = element_text(color = "black",size = 12)) +
  labs(x = "", y = "")

mxc_single |> 
  filter(cause == "All",
         year == "2016-2019",
         educ == "Total") |> 
  pivot_wider(names_from = sex, values_from = mx) |> 
  mutate(`M vs F` = arriaga(Males, Females),
         `F vs M` = -arriaga(Females, Males)) |> 
  select(age, `M vs F`, `F vs M`) |> 
  pivot_longer(-age, names_to = "Arriaga\ndirection",values_to = "result") |> 
  ggplot(aes(x = age, y = result, color = `Arriaga\ndirection`)) +
  geom_line(linewidth = 2) +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))


composed <-
  decomp_total |>
  filter(year == "2016-2019") |>
  mutate(cause = case_when(
    !cause %in% c(
      "External",
      "Circulatory",
      "Digestive",
      "Neoplasms",
      "Respiratory"
    ) ~ "Other",
    TRUE ~ cause
  )) |>
  group_by(educ, year, cause, age) |>
  summarise(valuersc = sum(result_rescaled), 
            value = sum(result),.groups = "drop") |>
  mutate(cause = factor(cause),
         cause = reorder(cause, -valuersc)) |>
  mutate(
    educ = as.factor(educ),
    educ = fct_relevel(educ, "Higher", after = Inf)) |> 
  ggplot(aes(x = age, y = valuersc, fill = educ)) +
  geom_density(stat = "identity",
               alpha = 0.8,
               color = "white",
               position = "stack",
               linewidth = .25) +
  facet_grid(rows = vars(cause)) +
  theme_bw() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    axis.text = element_text(color = "black")) +
  labs(fill = "Education groups")
# --------------------------------------------------------- #
# combine figures into 1
ggarrange(
  before,
  after,
  composed,
  ncol = 3,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom",
  labels = c("education-specific decomposition", "total decomposition","total decomp, composed"),
  vjust = 1
)
# ----------------------------------------------------------------------- #
# combine figures into 1
# ggarrange(
#   before,
#   after,
#   ncol = 2,
#   nrow = 1,
#   common.legend = TRUE,
#   legend = "bottom",
#   labels = c("subgroup-specific", "total composition"),
#   vjust = 1
# )
# ----------------------------------------------------------------------- #
# margins are interesting:


st_bind <- kit |> 
  ungroup() |> 
  summarize(margin = sum(st_component)) |> 
  mutate(educ = "Educ. Composition", .before = 1)

# educational component
# the one with causes and component
# change the figure with stack
# say in caption that we decompose the stationary
# sum st component for education
# figure with aggregate and weighted we decompose the purple

decomp_total |> 
  filter(year == "2016-2019") |> 
  group_by(educ) |> 
  summarize(margin = sum(result_rescaled)) |>
  full_join(dplyr::select(st_gaps[1, ], margin = cc_str) %>% 
              mutate(educ = "Educ. component")) %>%
  mutate(educ = fct_relevel(educ, "Higher", after = Inf)) |> 
  ggplot(aes(y = educ, 
             x = margin,
             fill = educ)) +
  geom_col() +
  guides(color= "none") +
  theme_minimal() +
  xlab("Contribution to sex-gap")+
  scale_x_continuous(breaks = pretty_breaks())+
  xlab("Contribution to sex-gap") + 
    ylab("Education") +
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text =element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black"))


decomp_total |> 
  filter(year == "2016-2019") |> 
  group_by(cause) |> 
  summarize(margin = sum(result_rescaled)) |> 
  filter(cause != "Covid-19") |> 
  full_join(dplyr::select(st_gaps[1, ], margin = cc_str) %>% 
              mutate(cause = "Educ. component")) %>% 
  mutate(margin_sign = if_else(sign(margin) == 1, "#eb4034","#3483eb")) |>
  mutate(margin_sign = if_else(cause == "Educ. component", "#BF2C34", margin_sign)) |>
  ggplot(aes(y = reorder(cause, margin), 
             x = margin, 
             color = margin_sign,
             fill = margin_sign)) +
  geom_col() +
  guides(color= "none") +
  theme_minimal() +
  scale_color_identity() +
  scale_fill_identity()+
  xlab("contribution to sex-gap")
# ----------------------------------------------------------------------- #
# check gaps
mort_gaps <- decomp_total |> 
  group_by(year) |> 
  summarize(cc_mort = sum(result_rescaled))

st_gaps <- kit |> 
  group_by(year) |> 
  summarize(cc_str = sum(st_component))

dec_gaps <- full_join(mort_gaps, st_gaps, by = join_by(year)) |> 
  mutate(dec_gap = cc_str + cc_mort)

# exact match to gaps
dec_gaps
gaps

non_stationary_gap <- 
  e35_total_compare |> 
  pivot_wider(names_from = sex, values_from = ex, names_prefix = "e35_") |> 
  mutate(gap = e35_Females - e35_Males)

e35_kit |> 
  select(educ, year, e35_Females, e35_Males) |> 
  mutate(gap = e35_Females - e35_Males) 
# ----------------------------------------------------------------------- #
# TODO:
# 1: include Total from stationary 'gaps' in this plot on the right or left.
ed_gps <- e35_kit |> 
  select(educ, year, Females = e35_Females, Males = e35_Males)  |>      
  pivot_longer(Females:Males, names_to = "sex", values_to = "e35") |> 
  filter(year == "2016-2019")

tot_gps <- gaps %>% 
  filter(year == "2016-2019") %>% 
  dplyr::select(year, Females = e35_Females, Males = e35_Males) %>% 
  pivot_longer(-year,
               names_to = "sex",
               values_to = "e35") %>% 
  mutate(educ = "Total")

ed_gps %>% 
  full_join(tot_gps) %>%
  ggplot(aes(x = educ, y = e35, fill = sex)) +
  geom_col(position = position_dodge(), color = "white") +
  theme_bw() + 
  coord_flip()+
  scale_y_continuous(breaks = pretty_breaks(n = 12)) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_blank(),
    axis.text = element_text(color = "black")) +
  labs(fill = "Education groups") + 
  labs(x = "Education level", y = "Age") + 
  ggtitle("Life expectancy at the are of 35 by education group")
# ----------------------------------------------------------------------- #
# 2: make a plot comparing e35 stationary and non-stationary
e35_total_compare %>% 
  filter(year == "2016-2019") %>%
  mutate(type = "Non-stationary") %>%
  rename(e35 = ex) %>% 
  full_join(tot_gps) %>%
  dplyr::select(-educ) %>% 
  mutate(type = ifelse(is.na(type), "Stationary", type)) %>% 
  ggplot(aes(x = sex, y = e35, fill = type)) +
  geom_col(position = position_dodge(), color = "white") +
  theme_bw() + 
  coord_flip()+
  scale_y_continuous(breaks = pretty_breaks(n = 12)) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_blank(),
    axis.text = element_text(color = "black")) +
  labs(fill = "Education groups") + 
  labs(x = "Education level", y = "Age") + 
  ggtitle("Life expectancy at the age of 35 in stationary and non-stationary populations")
# ----------------------------------------------------------------------- #
# 3: make a plot of education-specific gaps and the stationary gap
education <- e35_kit |> 
  select(educ, year, e35_diff)  |>
  filter(year == "2016-2019") %>% 
  mutate(type = "By education") %>% 
  rename(gap = e35_diff)

tot_gps <- gaps %>% 
  filter(year == "2016-2019") %>% 
  dplyr::select(year, gap) %>%
  mutate(educ = "Total") %>% 
  mutate(type = "Stationary")

orig_gap <- e35_total_compare |> 
  pivot_wider(names_from = sex, values_from = ex, names_prefix = "e35_") |> 
  mutate(gap = e35_Females - e35_Males) %>% 
  filter(year == "2016-2019") %>%
  dplyr::select(year, gap) %>% 
  mutate(type = "Non-Stationary") %>% 
  mutate(educ = "Total")
  
education %>% 
  full_join(tot_gps) %>% 
  full_join(orig_gap) %>% 
  ggplot(aes(x = type, y = gap, fill = educ)) + 
  geom_col(position = position_dodge(), color = "white") + 
  theme_minimal() + 
  coord_flip()+
  scale_y_continuous(breaks = pretty_breaks(n = 12)) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_blank(),
    axis.text = element_text(color = "black")) +
  labs(fill = "Education groupse") +
  labs(x = "e(35) difference type", 
       y = "Difference in years") + 
  ggtitle("The Female-Male diference in life expectancy at the age of 35 by education and population type")
# ----------------------------------------------------------------------- #
# 4: copy all plots into a Google presentation, link to be shared by email.
# ----------------------------------------------------------------------- #