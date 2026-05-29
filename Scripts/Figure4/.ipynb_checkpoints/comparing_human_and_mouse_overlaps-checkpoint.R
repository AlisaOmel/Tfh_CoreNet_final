library(dplyr)

human_translate <- read.csv("outputs/upset_membership_df_humantranslate.csv")
mouse_translate <- read.csv("outputs/upset_membership_df_mousetranslate.csv")

# Genes in CoreNet.DI ∩ Crotty.2020 for each species
human_hits <- human_translate$Genes[human_translate$CoreNet.DI == 1 &
                                      human_translate$Crotty.2020 == 1]

mouse_hits <- mouse_translate$Genes[mouse_translate$CoreNet.DI == 1 &
                                      mouse_translate$Crotty.2020 == 1]

cat("Human hits (", length(human_hits), "):\n", sep = "")
print(sort(human_hits))

cat("\nMouse hits (", length(mouse_hits), "):\n", sep = "")
print(sort(mouse_hits))