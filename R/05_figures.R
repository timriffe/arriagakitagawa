source("R/00_initial_data_preparation.R")
source("R/01_smoothing_and_ungroupping.R")
source("R/02_LT_and_average_quantities.R")
source("R/03_decomposition_results.R")
source("R/04_gaps.R")

# table
mxc_decomp %>%
  group_by(educ, year) %>%
  summarise(decomp = sum(result)) %>% 
  full_join(struct_kit) %>%
  full_join(e35_kit) %>%
  dplyr::select(-c(st_mean, st_mean)) %>% 
  relocate(decomp, .after = e35_diff) %>%
  relocate(e35_avg, .after = e35_Males) %>%
  set_names(c("Education", 
              "Period", 
              "C(x) Females",
              "C(x) Males",
              "C(x) Difference",
              "e(35) Females", 
              "e(35) Males", 
              "e(35) Average",
              "e(35) Difference",
              "Decomposition")) %>% 
  mutate(across(contains("C(x)"), ~ .x %>% 
                  round(3))) %>% 
  mutate(across(contains("e(35)"), ~ .x %>% 
                  round(2))) %>% 
  mutate(Decomposition = round(Decomposition, 2)) %>%
  xtable()

# causes 
decomp_total |> 
filter(year == "2016-2019") |> 
  group_by(cause) |> 
  summarize(margin = sum(result_rescaled)) |> 
  filter(cause != "Covid-19") |> 
  full_join(dplyr::select(st_gaps[1, ], margin = cc_str) %>% 
              mutate(cause = "Educ. component")) %>% 
  mutate(margin_sign = if_else(sign(margin) == 1, "#F8766D","#C77CFF")) |>
  mutate(margin_sign = if_else(cause == "Educ. component", "#00BFC4", margin_sign)) %>% 
  ggplot(aes(y = reorder(cause, margin), 
             x = margin, 
             color = margin_sign,
             fill = margin_sign)) +
  geom_col() +
  guides(color= "none") +
  theme_minimal() +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(breaks = pretty_breaks())+
  xlab("Contribution to sex-gap") + 
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text =element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black"))


# causes post covid
decomp_total |> 
  filter(year != "2016-2019") |> 
  group_by(cause) |> 
  summarize(margin = sum(result_rescaled)) |> 
  filter(cause != "Covid-19") |> 
  full_join(dplyr::select(st_gaps[1, ], margin = cc_str) %>% 
              mutate(cause = "Educ. component")) %>% 
  mutate(margin_sign = if_else(sign(margin) == 1, "#F8766D","#C77CFF")) |>
  mutate(margin_sign = if_else(cause == "Educ. component", "#00BFC4", margin_sign)) %>% 
  ggplot(aes(y = reorder(cause, margin), 
             x = margin, 
             color = margin_sign,
             fill = margin_sign)) +
  geom_col() +
  guides(color= "none") +
  theme_minimal() +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(breaks = pretty_breaks())+
  xlab("Contribution to sex-gap") + 
  theme(
    legend.position = "none",
    axis.title = element_blank(),
    strip.background = element_blank(),
    strip.text =element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black"))

# education | Stationary, non-stationary
education %>% 
  full_join(tot_gps) %>% 
  full_join(orig_gap) %>%
  mutate(labels = educ) %>% 
  mutate(educ = ifelse(type == "Stationary", "Total (Stationary)", educ)) %>% 
  mutate(educ = ifelse(type == "Non-Stationary", "Total (Non-Stationary)", educ)) %>%
  mutate(educ = factor(educ, levels = c("Total (Stationary)", "Total (Non-Stationary)", 
                                        "Higher", "Secondary", "Primary"))) %>% 
  mutate(type = ifelse(str_detect(type, "By education"), "Education", "e(35) variant")) %>% 
  mutate(type = factor(type, levels = c("Education", "e(35) variant"))) %>%
  filter(year == "2016-2019") %>% 
  ggplot(aes(x = educ, y = gap, fill = educ)) + 
  geom_col(position = position_dodge(), color = "white") + 
  theme_minimal() + 
  coord_flip()+
  facet_grid(type ~., switch = "y",
             scales = "free_y",
             space = "free_y")  +
  scale_y_continuous(breaks = pretty_breaks()) +
  # scale_fill_discrete(labels = c("Higher", "Primary", "Secondary", "Total", "Total")) + 
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    strip.placement = "outside",
    strip.text =  element_text(face = "bold", color = "black", size = 10),
    axis.text = element_text(face = "bold", color = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", color = "black"),
    legend.text = element_text(face = "bold", color = "black")) +
  labs(fill = "Education groups",) + 
  labs(y = "Difference (years)") 

# post covid
  education %>% 
  full_join(tot_gps) %>% 
  full_join(orig_gap) %>%
  mutate(labels = educ) %>% 
  mutate(educ = ifelse(type == "Stationary", "Total (Stationary)", educ)) %>% 
  mutate(educ = ifelse(type == "Non-Stationary", "Total (Non-Stationary)", educ)) %>%
  mutate(educ = factor(educ, levels = c("Total (Stationary)", "Total (Non-Stationary)", 
                                        "Higher", "Secondary", "Primary"))) %>% 
  mutate(type = ifelse(str_detect(type, "By education"), "Education", "e(35) variant")) %>% 
  mutate(type = factor(type, levels = c("Education", "e(35) variant"))) %>%
  filter(year != "2016-2019") %>% 
  ggplot(aes(x = educ, y = gap, fill = educ)) + 
  geom_col(position = position_dodge(), color = "white") + 
  theme_minimal() + 
  coord_flip()+
  facet_grid(type ~., switch = "y",
             scales = "free_y",
             space = "free_y")  +
  scale_y_continuous(breaks = pretty_breaks()) +
  # scale_fill_discrete(labels = c("Higher", "Primary", "Secondary", "Total", "Total")) + 
  theme(
    legend.position = "none",
    legend.direction = "horizontal",
    strip.placement = "outside",
    strip.text =  element_text(face = "bold", color = "black", size = 12),
    axis.text = element_text(face = "bold", color = "black", size = 12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", color = "black", size = 12),
    legend.text = element_text(face = "bold", color = "black")) +
  labs(fill = "Education groups",) + 
  labs(y = "Difference (years)") 

  # composed figure
  fig <- decomp_total |>
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
              value = sum(result), .groups = "drop") |>
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
    facet_grid(rows = vars(cause), switch = "y", labeller = label_value) +
    theme_minimal() +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_fill_manual(values = c("#e86bf3", "#00b0f6", "#00bf7c")) +
    coord_flip() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      # axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.x = element_text(face = "bold", color = "black"),
      strip.text.y.left = element_text(face = "bold", color = "black", angle = 0, hjust = 0, size = 12),  # Set angle to 0 and adjust hjust
      legend.text = element_text(face = "bold", color = "black", size = 12),
      legend.title = element_text(face = "bold", color = "black", size = 14),
      axis.text = element_text(color = "black")
    ) +
    labs(fill = "Education groups")
  
  
  ggsave("fig4.pdf", fig, width = 10, height = 12)
  
  
  # composed figure post covid
  decomp_total |>
    filter(year != "2016-2019") |>
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
              value = sum(result), .groups = "drop") |>
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
    facet_grid(rows = vars(cause), switch = "y", labeller = label_value) +
    theme_minimal() +
    scale_x_continuous(breaks = pretty_breaks()) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_fill_manual(values = c("#e86bf3", "#00b0f6", "#00bf7c")) +
    coord_flip() +
    theme(
      legend.position = "bottom",
      axis.title = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      # axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.x = element_text(face = "bold", color = "black"),
      strip.text.y.left = element_text(face = "bold", color = "black", angle = 0, hjust = 0, size = 12),  # Set angle to 0 and adjust hjust
      legend.text = element_text(face = "bold", color = "black", size = 10),
      legend.title = element_text(face = "bold", color = "black", size = 14),
      axis.text = element_text(color = "black")
    ) +
    labs(fill = "Education groups")

# education gaps
  decomp_total |>
    filter(year == "2016-2019") |>
    group_by(educ) |> 
    summarize(margin = sum(result_rescaled)) |>
    full_join(dplyr::select(st_gaps[1, ], margin = cc_str) |> 
                mutate(educ = "Educ. component")) |>
    mutate(educ = factor(educ, levels = c("Educ. component", 
                                          "Higher", "Secondary", 
                                          "Primary"))) |>  
    ggplot(aes(y = educ, 
               x = margin,
               fill = educ)) +
    geom_col() +
    guides(color= "none") +
    theme_minimal() +
    xlab("Contribution to sex-gap")+
    scale_x_continuous(breaks = pretty_breaks())+
    scale_fill_manual(values = c("#00BFC4", "#00bf7c", "#00b0f6", "#e86bf3"))+
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


# education gaps post covid
  decomp_total |>
    filter(year != "2016-2019") |>
    group_by(educ) |> 
    summarize(margin = sum(result_rescaled)) |>
    full_join(dplyr::select(st_gaps[1, ], margin = cc_str) |> 
                mutate(educ = "Educ. component")) |>
    mutate(educ = factor(educ, levels = c("Educ. component", 
                                          "Higher", "Secondary", 
                                          "Primary"))) |>  
    ggplot(aes(y = educ, 
               x = margin,
               fill = educ)) +
    geom_col() +
    guides(color= "none") +
    theme_minimal() +
    xlab("Contribution to sex-gap")+
    scale_x_continuous(breaks = pretty_breaks())+
    scale_fill_manual(values = c("#00BFC4", "#00bf7c", "#00b0f6", "#e86bf3"))+
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
