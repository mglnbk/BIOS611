library(xgboost)
library(caret)
library(tidyverse)
library(reshape2)

engineer_features <- function(df) {
  ranks <- df %>%
    select(matches("^R")) %>%
    as.matrix()
  suits <- df %>%
    select(matches("^S")) %>%
    as.matrix()

  # sort: (e.g., 10, 2, 5 -> 2, 5, 10), make it easy to see order
  ranks_sorted <- t(apply(ranks, 1, sort))
  colnames(ranks_sorted) <- paste0("Rank_Sorted_", 1:5)

  # feature engineer
  # mannually extract valid information
  feat_list <- apply(df, 1, function(row) {
    r_vals <- as.numeric(row[c(2, 4, 6, 8, 10)])
    s_vals <- as.numeric(row[c(1, 3, 5, 7, 9)])

    r_counts <- sort(table(r_vals), decreasing = TRUE)
    s_counts <- sort(table(s_vals), decreasing = TRUE)

    # Feature: max count of a single rank
    max_rank_count <- as.numeric(r_counts[1])
    # Feature: second max count
    sec_rank_count <- if(length(r_counts) > 1) as.numeric(r_counts[2]) else 0

    # Feature: max suit count
    max_suit_count <- as.numeric(s_counts[1])

    # Feature: how many unique ranks
    n_unique_ranks <- length(r_counts)

    # Feature: straight
    sorted_r <- sort(r_vals)
    is_straight_raw <- (sorted_r[5] - sorted_r[1] == 4) && (n_unique_ranks == 5)

    c(max_rank_count = max_rank_count,
      sec_rank_count = sec_rank_count,
      max_suit_count = max_suit_count,
      n_unique_ranks = n_unique_ranks,
      rank_span = sorted_r[5] - sorted_r[1])
  })

  feat_df <- as.data.frame(t(feat_list))

  # combine old and new
  final_df <- cbind(as.data.frame(ranks_sorted), feat_df)
  return(final_df)
}

# Data loading
cat("Data Loading...\n")
cols <- c("S1", "R1", "S2", "R2", "S3", "R3", "S4", "R4", "S5", "R5", "Class")
train_raw <- read.csv(
  "data/poker-hand-training-true.data",
  header = FALSE,
  col.names = cols
)
test_raw  <- read.csv(
  "data/poker-hand-testing.data",
  header = FALSE,
  col.names = cols
)

# Feature Engineering
cat("Begin feature engineering, this step may take relatively long time...\n")
cat("Engineering features for Training data...\n")
X_train_eng <- engineer_features(train_raw)
y_train <- train_raw$Class

cat("Engineering features for Test data...\n")
X_test_eng <- engineer_features(test_raw)
y_test <- test_raw$Class

X_train_mat <- as.matrix(X_train_eng)
X_test_mat  <- as.matrix(X_test_eng)

# use weight
cls_count <- table(y_train)
cls_weight <- max(cls_count) / (cls_count + 1) # prevent rare is missing care
sample_weight <- cls_weight[as.character(y_train)]

dtrain <- xgb.DMatrix(data = X_train_mat,
                      label = y_train,
                      weight = sample_weight)
dtest  <- xgb.DMatrix(data = X_test_mat,
                      label = y_test)

# XGboost parameters
params <- list(
  booster = "gbtree",
  objective = "multi:softprob",
  num_class = 10,
  eta = 0.1,
  max_depth = 8,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  gamma = 0.1,
  eval_metric = "merror"
)

set.seed(2025611)
cat("Training XGBoost Model...\n")
model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 1000,
  watchlist = list(train = dtrain, test = dtest),
  early_stopping_rounds = 50,
  print_every_n = 50
)

# preditc
cat("Predicting....\n")
pred_prob <- predict(model, dtest)
pred_mat  <- matrix(pred_prob, ncol = 10, byrow = TRUE)
pred_lab  <- max.col(pred_mat) - 1 # 1-10 2 0-9

# confusion Matrix
conf_matrix <- confusionMatrix(
  factor(pred_lab, levels=0:9), 
  factor(y_test, levels=0:9)
)

# glimpse at the result
cat("glimpse at the accuracy, though not comprehensive metric...\n")
print(conf_matrix$overall["Accuracy"])

# save model
timestamp <- format(Sys.time(), "%y%m%d%H%M%S")
filename <- paste0("saved_model/poker_prediction_model_", timestamp, ".model")
xgb.save(model, filename)
xgb.save(model, "saved_model/latest_model.model")

# plotting
cat("Plotting...\n")
mat <- as.matrix(conf_matrix$table)
df_plot <- melt(mat)
colnames(df_plot) <- c("reference", "prediction", "freq")

# log scale for fill helps visualize rare classes better
cm_plot <- ggplot(df_plot, aes(x = reference, y = prediction, fill = freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = freq), color = "white", size = 3) +
  scale_fill_gradient(low = "steelblue1", high = "darkblue", trans = "log1p") +
  theme_minimal(base_size = 14) +
  labs(title = "Confusion matrix",
       x = "True class", y = "Predicted class")


# save plot
ggsave(
  "plot/confusion_matrix_engineered.png",
  plot = cm_plot,
  width = 12,
  height = 12,
  units = "in",
  dpi = 300
)
