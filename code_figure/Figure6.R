library(readxl)
library(tidyverse)
library(ggpubr)
library(patchwork)

df <- read_excel("results/荧光染色统计.xlsx")
df_clean <- df %>% 
  fill(切片, .direction = "down")

head(df_clean)

my_theme <- theme_classic() +
  theme(
    text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none" 
  )

my_colors <- c("LS" = "#E69F00", "SS" = "#56B4E9")

p1 <- ggplot(df_clean, aes(x = 分组, y = `OLIG2+占比`, fill = 分组)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) + 
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +     
  scale_fill_manual(values = my_colors) +
  labs(title = "Ki67+ OLIG2+ Proliferation", 
       y = "Proportion (OLIG2+ / Ki67+)", 
       x = "Group") +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("LS", "SS")),
                     label = "p.signif") + 
  my_theme

p2 <- ggplot(df_clean, aes(x = 分组, y = `CD11b+占比`, fill = 分组)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 2.5, alpha = 0.6) +
  scale_fill_manual(values = my_colors) +
  labs(title = "Ki67+ CD11b+ Proliferation", 
       y = "Proportion (CD11b+ / Ki67+)", 
       x = "Group") +
  stat_compare_means(method = "t.test", 
                     comparisons = list(c("LS", "SS")),
                     label = "p.signif") + 
  my_theme

combined_plot <- p1 + p2
print(combined_plot)

ggsave("figures/Fig6B.pdf", combined_plot, width = 10, height = 5)
