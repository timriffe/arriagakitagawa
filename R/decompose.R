
library(coddecomp)
library(tidyverse)

mxc <- read_csv("data/mxc.csv") |> 
  filter(year == 2019)

head(mxc)

