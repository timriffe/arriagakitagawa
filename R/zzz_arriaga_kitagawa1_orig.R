# ----------------------------------------------------------------------- #
library(tidyverse)
library(coddecomp)
library(ggridges)
library(scales)
library(ggpubr)
library(mgcv)
# ----------------------------------------------------------------------- #
# data from Sergi
load("data_esp_1621.RData")
# ----------------------------------------------------------------------- #
# I will use ths data further for predict() with GAM model results
new_data <- expand_grid(
  age5  = seq(35, 100, 1),
  causa = unique(data_esp_1621$causa),
  year  = c("2016-2019", "2020-2021"),
  educ  = unique(data_esp_1621$educ),
  sex   = unique(data_esp_1621$sex),
  pop   = 1) %>%
  # these categories are redundant
  filter(year  != "Total") %>%
  # for model
  mutate(year = as.factor(year)) %>%
  # there are 2 secondary education group. I combine them into 1
  mutate(educ = ifelse(str_detect(educ, "Secundaria"), "Secundaria", educ)) %>%
  distinct()
# ----------------------------------------------------------------------- #
# Cleaning and rearranging original data original data 
tst <- data_esp_1621 %>%
  as_tibble() %>%
  # redundant category
  filter(year  != "Total") %>% 
  # turn age and year to numeric 
  mutate(age5 = parse_number(as.character(age5)),
         year = as.numeric(year),
         # group years into 2 periods pre and post covid
         year = ifelse(between(year, 2016, 2019), "2016-2019", "2020-2021"),
         year = as.factor(year),
         # combine secondary education into 1 
         educ = ifelse(str_detect(educ, "Secundaria"), "Secundaria", educ)) %>%
  group_by(sex, educ, year, causa, age5) %>% 
  # add up data to match new periods and education groups
  summarise(deaths  = sum(deaths),
            pop     = sum(pop),
            .groups = "drop") %>%
  # nest by sex education groups and causes 
  # unfortunately all together do not work, I tried
  # the confidence intervals in this case are too narrow
  # and variability is very low
  # so we fit separately for sex, educ group and cause
  group_nest(sex, educ, causa) %>%
  # add data for predict()
  nest_join(new_data, by = c("causa", "sex", "educ")) %>%
  # GAM model with log pop as offset and p-spline on age
  mutate(model = map(
    data,
    ~ gam(
      deaths ~ s(age5, bs = "ps") + year + offset(log(pop)),
      data   = .,
      family = quasipoisson
    )
  )) %>%
  # here predict() happens
  mutate(pred = map2(.x = model, 
                     .y = new_data, ~ .y %>%
                       mutate(new = predict(.x, 
                                            newdata = .,
                                            type    = "response"))
  )) %>%
  # remove intermediate columns and unnest
  dplyr::select(-c(data:model)) %>%
  unnest(pred) %>% 
  dplyr::select(-pop)
# ----------------------------------------------------------------------- #
# calculate life expectancy
Lt <- tst %>%
  group_by(sex, year, educ, age5) %>%
  summarise(mx = sum(new)) %>%
  group_by(sex, year, educ) %>%
  mutate(n  = c(diff(as.numeric(unique(age5))), 999),
         ax = 0.5) %>%
  nest() %>%
  mutate(data = map(data, ~ .x       %>%
                      mutate(qx = n * mx / (1 + (n - ax) * mx),
                             qx = ifelse(n == 999, 1, qx),
                             px = 1 - qx,
                             lx = head(c(1, cumprod(px)), -1) * 100000,
                             dx = lx - c(tail(lx, -1), 0),
                             Lx = c(head(n * c(tail(lx, -1), 0), -1), tail(lx, 1) / tail(mx, 1)),
                             Tx = rev(cumsum(rev(Lx))),
                             ex = Tx / lx) %>%
                      ungroup())) %>%
  unnest(data) %>%
  ungroup()
# ----------------------------------------------------------------------- #
# Average male and female overall mortality for each age, year and educ
avgs <- tst %>%
  filter(causa == "All",
         sex   != "Total") %>%
  group_by(educ, year, age5) %>% 
  summarise(avg = mean(new), .groups = "drop")
# ----------------------------------------------------------------------- #
# feed the averages to the corresponding function
pt1 <- avgs %>% 
  group_nest(educ, year) %>% 
  mutate(data = map(data, ~ .x %>%
                      pull(avg) %>% 
                      sen_arriaga_instantaneous() %>% 
                      as_tibble() %>% 
                      mutate(age5 = unique(avgs$age5)))) %>% 
  unnest(data)
# ----------------------------------------------------------------------- #
# cause, age, year, and education specific diff of Males - Females 
# multiplied by the results from sen_arriaga_instantaneous
pt2 <- tst %>% 
  filter(sex   != "Total", 
         causa != "All") %>%
  pivot_wider(names_from  = sex,
              values_from = new) %>% 
  mutate(dif = Males - Females) %>% 
  full_join(pt1) %>%
  mutate(result = dif * value) %>% 
  dplyr::select(-c(Females:value))
# ----------------------------------------------------------------------- #
# now make a matrix with causes in columns and ages in rows
# by groups
zz <- pt2 %>%
  pivot_wider(names_from  = "causa",
              values_from = "result") %>% 
  dplyr::select(-age5) %>%
  group_nest(educ, year) %>%
  mutate(data = map(data, ~ .x %>% 
                      as.matrix()))
# ----------------------------------------------------------------------- #
# kitagawa part: ex by groups
kitag <- Lt %>%
  filter(sex  != "Total",
         age5 == min(age5)) %>% 
  dplyr::select(sex, educ, year, ex) 
# ----------------------------------------------------------------------- #
# population prevalence data from original source
# mostly data cleaning similar to step 1
pop_k <- data_esp_1621 %>% 
  as_tibble() %>%
  filter(year  != "Total",
         sex  != "Total") %>% 
  mutate(age5 = parse_number(as.character(age5)),
         year = as.numeric(year),
         year = ifelse(between(year, 2016, 2019), "2016-2019", "2020-2021"),
         year = as.factor(year),
         educ = ifelse(str_detect(educ, "Secundaria"), "Secundaria", educ)) %>%
  group_by(sex, educ, year, causa, age5) %>% 
  summarise(deaths  = sum(deaths),
            pop     = sum(pop),
            .groups = "drop") %>% 
  group_by(sex, educ, year) %>% 
  # population by groups
  summarise(pop = sum(pop)) %>% 
  group_by(sex) %>%
  # scale populaton
  mutate(pop = pop / sum(pop)) %>% 
  ungroup()
# ----------------------------------------------------------------------- #
# population prevalence by groups
pop_new <- pop_k %>% 
  group_by(educ, year) %>% 
  summarise(pop = mean(pop), .groups = "drop")
# ----------------------------------------------------------------------- #
# calculate the kitagawa rate component for ex
pt3 <- kitag %>% 
  pivot_wider(names_from = sex,
              values_from = ex) %>% 
  mutate(df = Males - Females) %>% 
  full_join(pop_new) %>% 
  mutate(res = df * pop) %>% 
  dplyr::select(-c(Females:pop))
# ----------------------------------------------------------------------- #
# decomposition result
final <- zz %>% 
  full_join(pt3) %>% 
  # divide each matrix by its sum and go back to tibble
  mutate(data = map(data, ~ (.x / sum(.x)) %>% 
                      as_tibble() %>% 
                      mutate(age = unique(pt2$age5)) %>% 
                      pivot_longer(-c(age),
                                   names_to  = "cause",
                                   values_to = "value"))) %>% 
  unnest(data) %>% 
  # multiply the resulting data by corresponding kitagawa component
  mutate(value = value * res) %>% 
  dplyr::select(-res)
# ----------------------------------------------------------------------- #
# separate figures for 2 periods
a <-
  final %>%
  filter(year == "2016-2019") %>%
  filter(educ != "Total") %>%
  mutate(cause = str_remove_all(cause, "^[:digit:]+. ")) %>%
  mutate(cause = case_when(
    !cause %in% c(
      "External",
      "Circulatory",
      "Digestive",
      "Neoplasms",
      "Respiratory"
    ) ~ "Other",
    TRUE ~ cause
  )) %>%
  group_by(educ, year, cause, age) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(cause = factor(
    cause,
    levels = c(
      "Neoplasms",
      "Circulatory",
      "External",
      "Respiratory",
      "Digestive",
      "Other"
    )
  )) %>%
  mutate(
    educ = case_when(
      educ == "Primaria"  ~ "Primary",
      educ == "Secundaria"  ~ "Secondary",
      educ == "Superior"  ~ "Higher"
    )
  ) %>%
  ggplot(aes(x = age, y = value, fill = educ)) +
  geom_density(stat = "identity",
               alpha = 0.8,
               color = "white") +
  facet_wrap( ~ cause) +
  theme_bw() +
  coord_flip() +
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
    axis.text = element_text(color = "black")
  ) +
  labs(fill = "Education groups")
# ----------------------------------------------------------------------- #
b <- final %>%
  filter(year == "2020-2021") %>%
  filter(educ != "Total") %>%
  mutate(cause = str_remove_all(cause, "^[:digit:]+. ")) %>%
  mutate(cause = case_when(
    !cause %in% c(
      "External",
      "Circulatory",
      "Neoplasms",
      "Respiratory",
      "Covid-19"
    ) ~ "Other",
    TRUE ~ cause
  )) %>%
  group_by(educ, year, cause, age) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  mutate(cause = factor(
    cause,
    levels = c(
      "Neoplasms",
      "Circulatory",
      "Covid-19",
      "External",
      "Respiratory",
      "Other"
    )
  )) %>%
  mutate(
    educ = case_when(
      educ == "Primaria"  ~ "Primary",
      educ == "Secundaria"  ~ "Secondary",
      educ == "Superior"  ~ "Higher"
    )
  ) %>%
  ggplot(aes(x = age, y = value, fill = educ)) +
  geom_density(stat = "identity",
               alpha = 0.8,
               color = "white") +
  facet_wrap( ~ cause) +
  theme_bw() +
  coord_flip() +
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
    axis.text = element_text(color = "black")
  ) +
  labs(fill = "Education groups")
# ----------------------------------------------------------------------- #
# combine figures into 1
ggarrange(
  a,
  b,
  ncol = 2,
  nrow = 1,
  common.legend = TRUE,
  legend = "bottom",
  labels = c("2016-2019", "2020-2021"),
  hjust = -1.7,
  vjust = 1
)
# ----------------------------------------------------------------------- #