library(tidyverse)
library(reshape2)
library(ggplot2)



# Load Data
cols <- c("S1", "R1", "S2", "R2", "S3", "R3", "S4", "R4", "S5", "R5", "Class")
df <- read.csv("data/poker-hand-training-true.data", header = FALSE, col.names = cols)

# class distribution
p1 <- ggplot(df, aes(x = as.factor(Class))) +
  geom_bar(fill = "steelblue") +
  scale_y_log10() +
  labs(title = "Figure 1: class distribution of training data (shown in log scale)",
       x = "Poker Hand Class (0-9)", y = "Count (Log Scale)")
ggsave("plot/Figure1.png", p1, dpi = 300)

# rank distriution across
ranks <- df %>% select(starts_with("R")) %>% melt()
p2 <- ggplot(ranks, aes(x = value)) +
  geom_histogram(binwidth = 1, fill = "coral", color = "white") +
  labs(title = "Figure 2: Distribution of Card Ranks (1-13)",
       subtitle = "Uniform distribution indicates reasonable deck",
       x = "Rank", y = "Frequency") +
  theme_minimal()
ggsave("plot/Figure2.png", p2, dpi = 300)

# suit correlation
cormat <- cor(df %>% select(starts_with("S")))
melted_cormat <- melt(cormat)
p3 <- ggplot(
  melted_cormat,
  aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(title = "Figure 3: Correlation among suits") +
  theme_minimal()
ggsave("plot/Figure3.png", p3, dpi = 300)

# ranks distribution
df$RankSum <- rowSums(df %>% select(starts_with("R")))
p4 <- ggplot(df, aes(x = as.factor(Class), y = RankSum)) +
  geom_boxplot(fill = "lightgreen") +
  labs(title = "Figure 4: Sum of ranks in each class",
       x = "Class", y = "Sum of 5 Card Ranks") +
  theme_minimal()
ggsave("plot/Figure4.png", p4, dpi = 300)
