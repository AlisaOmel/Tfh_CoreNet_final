#!/usr/bin/env Rscript
library(SLIDE)
library(foreach)
library(tidyverse)
yaml_path = 'cm_tfh_slide_seacells.yaml' 
input_params = yaml::read_yaml(yaml_path)
#checkDataParams(input_params)
optimizeSLIDE(input_params, sink_file = TRUE)

#SLIDEcv(yaml_path)
