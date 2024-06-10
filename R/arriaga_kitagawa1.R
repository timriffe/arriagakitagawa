# ----------------------------------------------------------------------- #
mx
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
  mutate(age = as.integer(age),
         n  = if_else(age == 100, 999, 1),
         ax = .5,
         qx = n * mx / (1 + (n - ax) * mx),
         qx = ifelse(n == 999, 1, qx),
         lx = DemoTools::lt_id_q_l(qx, radix = 1),
         dx = lx * qx,
         Lx = lx - (1 - ax) * dx,
         Tx = rev(cumsum(rev(Lx))),
         ex = Tx / lx) |> 
  ungroup()

# we compare this later with our weighted-average e35 values just to see.
e35_total_compare <- Lt |>
  filter(educ == "Total", age == 35, sex != "Total") |> 
  select(sex, year, ex)

# ----------------------------------------------------------------------- #
# Average male and female overall mortality for each age, year and educ
mean_mx <- mxc_single |>
  filter(cause == "All",
         sex   != "Total",
         educ  != "Total") |>
  group_by(educ, year, age) |> 
  summarise(mean_mx = mean(mx), .groups = "drop")

# ----------------------------------------------------------------------- #
# feed the averages to the corresponding sensitivity function
sensitivity_at_mean <- mean_mx |> 
  group_nest(educ, year) |> 
  mutate(data = map(data, ~ .x |>
                      pull(mean_mx) |> 
                      sen_arriaga_instantaneous() |> 
                      as_tibble() |> 
                      mutate(age = 35:100))) |> 
  unnest(data) |> 
  rename(sensitivity = value)
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
  full_join(sensitivity_at_mean, 
            by = join_by(educ, age, year)) |>
  mutate(result = mxc_diff * sensitivity) 
# ----------------------------------------------------------------------- #
# kitagawa part: ex by groups
e35_kit <- Lt |>
  filter(sex  != "Total",
         educ !="Total",
         age == 35) |> 
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
struct_kit
# ----------------------------------------------------------------------- #
# combine for Kitagawa
kit <- left_join(struct_kit, 
                 e35_kit, 
                 by = join_by(educ, year)) |> 
  mutate(e35_component = st_mean * e35_diff,
         st_component  = e35_avg * st_diff)

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
  facet_grid(vars(cause),vars(educ),switch = "y") +
  theme_bw() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size =14),
    strip.text.y.left = element_text(angle = 0),
    legend.text = element_text(face = "bold", color = "black", size =12),
    legend.title = element_text(face = "bold", color = "black", size =12),
    axis.text = element_text(color = "black", size =12)) +
  labs(fill = "Education groups")
before
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
  facet_grid(vars(cause),vars(educ),switch = "y") +
  theme_bw() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size =14),
    strip.text.y.left = element_text(angle = 0),
    legend.text = element_text(face = "bold", color = "black", size =12),
    legend.title = element_text(face = "bold", color = "black", size =12),
    axis.text = element_text(color = "black", size =12)) +
  labs(fill = "Education groups")
after
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
  facet_grid(rows = vars(cause), switch = "y") +
  theme_bw() +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    strip.background = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size =14),
    strip.text.y.left = element_text(angle = 0),
    legend.text = element_text(face = "bold", color = "black", size =12),
    legend.title = element_text(face = "bold", color = "black", size =12),
    axis.text = element_text(color = "black", size =12)) +
  labs(fill = "Education groups")
composed
# --------------------------------------------------------- #
# combine figures into 1
# ggarrange(
#   before,
#   after,
#   composed,
#   ncol = 3,
#   nrow = 1,
#   common.legend = TRUE,
#   legend = "bottom",
#   labels = c("education-specific decomposition", "total decomposition","total decomp, composed"),
#   vjust = 1
# )
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
#margins are interesting:
st_bind <-
  kit |> 
  ungroup() |> 
  summarize(margin = sum(st_component)) |> 
  mutate(educ = "Educ. Composition", .before = 1)

decomp_total |> 
  filter(year == "2016-2019") |> 
  group_by(educ) |> 
  summarize(margin = sum(result_rescaled)) |> 
  bind_rows(st_bind) |> 
  mutate(color_comp = if_else(educ == "Educ. Composition","#a13b8e","#eb4034")) |> 
  mutate(educ = fct_relevel(educ, "Higher", after = Inf)) |> 
  ggplot(aes(x = educ, 
             y = margin,
             fill = color_comp)) +
  geom_col() +
    scale_fill_identity()+
  guides(color= "none") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  ylab("contribution to sex-gap") +
  xlab("") +
  guides(fill = "none")
st_tib <- tibble(cause = "Educ. Composition", margin = st_bind$margin)
decomp_total |> 
  filter(year == "2016-2019") |> 
  group_by(cause) |> 
  summarize(margin = sum(result_rescaled)) |> 
  bind_rows(st_tib) |> 
  filter(cause != "Covid-19") |> 
  mutate(margin_sign = if_else(sign(margin) == 1, "#eb4034","#3483eb"),
         margin_sign = if_else(cause == "Educ. Composition","#a13b8e",margin_sign)) |> 
  ggplot(aes(y = reorder(cause, margin), 
             x = margin, 
             color = margin_sign,
             fill = margin_sign)) +
  geom_col() +
  guides(color= "none") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  scale_color_identity() +
  scale_fill_identity()+
  xlab("contribution to sex-gap") +
  ylab("")
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
# ----------------------------------------------------------------------- #
# 2: make a plot comparing e35 stationary and non-stationary
e35_total_compare %>% 
  filter(year == "2016-2019") %>%
  mutate(type = "Aggregate") %>%
  rename(e35 = ex) %>% 
  full_join(tot_gps,by = join_by(sex, year, e35)) %>%
  dplyr::select(-educ) %>% 
  mutate(type = ifelse(is.na(type), "Weighted by age 35 education", type)) %>% 
  ggplot(aes(x = sex, y = e35, fill = type)) +
  geom_col(position = position_dodge(), color = "white") +
  theme_bw() + 
  scale_y_continuous(breaks = pretty_breaks(n = 12)) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", color = "black", size = 14),
    axis.text.y = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black", size = 12),
    legend.title = element_blank(),
    axis.text = element_text(color = "black", size =12)) +
  labs(fill = "Education groups") + 
  labs(x = "", y = "e35") + 
  scale_fill_manual(values = c("#c99c4d","#c253b9"))
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

# orig_gap <- e35_total_compare |> 
#   pivot_wider(names_from = sex, values_from = ex, names_prefix = "e35_") |> 
#   mutate(gap = e35_Females - e35_Males) %>% 
#   filter(year == "2016-2019") %>%
#   dplyr::select(year, gap) %>% 
#   mutate(type = "Non-Stationary") %>% 
#   mutate(educ = "Total") 
  
education %>% 
  bind_rows(tot_gps) %>% 
  mutate(
    educ = as.factor(educ),
    educ = relevel(educ, "Total"),
    educ = relevel(educ, "Primary"),
    educ = relevel(educ, "Secondary"),
    educ = relevel(educ, "Higher")) |> 
  # full_join(orig_gap) %>% 
  ggplot(aes(x = educ, y = gap, fill = educ)) + 
  geom_col(position = position_dodge(), color = "white") + 
  theme_bw() + 
  scale_y_continuous(breaks = pretty_breaks(n = 12)) +
  scale_color_viridis_b() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(face = "bold", color = "black", size = 14),
    axis.text.y = element_text(face = "bold", color = "black"),
    strip.text = element_text(face = "bold", color = "black", size = 14),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_blank(),
    axis.text = element_text(color = "black",size=12)) +
  labs(fill = "Education groups") +
  labs(x="")+
  guides(fill = "none")

# ----------------------------------------------------------------------- #
# 4: copy all plots into a Google presentation, link to be shared by email.
# ----------------------------------------------------------------------- #

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

arriaga

kit |> 
  ungroup() |> 
  filter(year == "2016-2019") |> 
  select(educ, Females = st_Females, Males = st_Males) |> 
  pivot_longer(-educ, names_to = "sex", values_to = "Educ. Composition") |> 
  mutate(
    educ = as.factor(educ),
    educ = relevel(educ,"Secondary"),
    educ = relevel(educ,"Primary")) |> 
  ggplot(aes(x=sex,y=`Educ. Composition`,fill=educ)) +
  geom_col() +
  theme_minimal() +
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size = 14)) +
  xlab("")
