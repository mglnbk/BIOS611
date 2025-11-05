library(tidyverse)
library(data.table)
library(cluster)

#### Do the task one ####
# n represents dimension and clusters (#dimension = #cluster)
# side_length represents how large the cube is 
# k represents how many pts generated in one cluster
generate_hypercube_clusters <- function(n, k, side_length, noise_sd = 1.0) {
  # centers are n dims * n dims
  centers <- diag(side_length, nrow = n, ncol = n) # n x n; row i is e_i * L
  
  data_list <- vector("list", n) # create null list of n size
  
  n_cluster_number <- n
  for (kth in 1:n_cluster_number) { # this n is cluster number
    center_idx <- matrix(centers[kth, ], nrow = 1, ncol = n)
    pts_in_kth_cluster <- matrix(rep(center_idx, k), nrow = k, ncol = n) + 
      matrix(rnorm(k*n, mean = 0, sd = noise_sd), nrow = k, ncol = n)
    data_list[[kth]] <- data.table(pts_in_kth_cluster) %>%
      mutate(true_cluster = kth)
  }
  
  df <- rbindlist(data_list)
  colnames(df)[1:n] <- paste0("dim", seq_len(n))
  df
}

generate_hypercube_clusters(n = 100, k = 5, side_length = 20)

dims_to_test <- c(6, 5, 4, 3, 2)
side_lengths_range <- 10:1  # 10 down to 1
number_pts_in_cluster <- 100
noise_sd_val <- 1.0


kmeans_robust <- function(x, k) {
  stats::kmeans(x, centers = k, nstart = 20, iter.max = 50)
}

task1_results <- expand_grid(
  n_dim = dims_to_test,
  side_length = side_lengths_range
) %>%
  mutate(
    # for each (n_dim, side_length), generate data and estimate optimal K via Gap Statistic
    est_k = map2_int(n_dim, side_length, function(n_dim, L) {
      dat <- generate_hypercube_clusters(
        n = n_dim,
        k = number_pts_in_cluster,
        side_length = L,
        noise_sd = noise_sd_val
      )
      # keep only the coordinate columns
      coords <- dat %>%
        select(starts_with("dim")) %>%
        as.matrix()
      
      # run Gap Statistic for k = 1,..., n_dim+2 just to be safe
      # ndim = n_clusters
      gap_obj <- clusGap(
        coords,
        FUNcluster = kmeans_robust,
        K.max = n_dim + 2,
        B = 20
      )
      
      best_k <- maxSE(
        f = gap_obj$Tab[, "gap"],
        SE.f = gap_obj$Tab[, "SE.sim"],
        method = "firstSEmax"
      )
      best_k
    })
  ) %>% mutate(true_k = n_dim)


tail(task1_results)
sum(task1_results$est_k == task1_results$true_k)/nrow(task1_results)

# Save .rds for plotting script
if (!dir.exists("results")) dir.create("results", recursive = TRUE)
saveRDS(task1_results, file = "results/task1_gap_results.rds")

p1 <- ggplot(task1_results,
         aes(x = side_length,
             y = est_k,
             group = n_dim,
             color = factor(n_dim))) +
    geom_line() +
    geom_point() +
    geom_hline(aes(yintercept = true_k),
               linetype = "dashed", color = "black") +
    scale_x_reverse(breaks = unique(task1_results$side_length)) +
    facet_wrap(~ n_dim, scales = "free_y") +
    labs(
      title = "Task 1",
      subtitle = "Dashed = # clusters / # corners / # dims",
      x = "side_length for each cube",
      y = "Estimated # clusters",
      color = "Dimension"
    ) +
    theme_minimal(base_size = 14)

ggsave("results/task1_gap_plot.png", p1, width = 8, height = 5, dpi = 300)

# Task 2. generate and classify on sphere

generate_shell_clusters <- function(
    n_shells,
    k_per_shell,
    max_radius,
    noise_sd = 0.1
  ) {
  
  inner_radius <- max_radius / (n_shells * 2) # a reasonable inner shell
  radii <- seq(inner_radius, max_radius, length.out = n_shells)
  
  all_pts <- vector("list", n_shells)
  
  for (s in seq_len(n_shells)) {
    r0 <- radii[s]
    # use polar index in 3D controlled by theta and phi
    # sample angles uniformly on shell
    theta <- runif(k_per_shell, 0, 2 * pi)
    phi <- acos(runif(k_per_shell, -1, 1)) # polar angle
    
    # radius with Gaussian noise (shell thickness)
    r     <- r0 + rnorm(k_per_shell, 0, noise_sd)
    
    x <- r*sin(phi)*cos(theta)
    y <- r*sin(phi)*sin(theta)
    z <- r*cos(phi)
    
    all_pts[[s]] <- data.table(
      x = x, y = y, z = z,
      true_shell = s
    )
  }
  
  rbindlist(all_pts)
}


library(tidyverse)
library(cluster)
library(plotly)

#---------------------------
# Build adjacency matrix given data matrix X and threshold d_threshold
#---------------------------
build_adjacency <- function(X, d_threshold) {
  # X: N x p
  # output: N x N matrix A of 0/1
  # We'll do brute force O(N^2). N here is ~400 so fine.
  n <- nrow(X)
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    dists <- sqrt(rowSums((t(t(X) - X[i, ]))^2))
    close_idx <- which(dists < d_threshold & dists > 0)
    A[i, close_idx] <- 1
    A[close_idx, i] <- 1
  }
  A
}

#---------------------------
# Spectral clustering core
#---------------------------
spectral_cluster_once <- function(X, k, d_threshold = 0.6) {
  
  A <- build_adjacency(X, d_threshold = d_threshold)
  
  degs <- rowSums(A)
  D <- diag(degs)
  L <- D - A
  
  # L_sym = D^{-1/2} L D^{-1/2}
  inv_sqrt_degs <- ifelse(degs > 0, 1 / sqrt(degs), 0)
  D_inv_sqrt <- diag(inv_sqrt_degs)
  L_sym <- D_inv_sqrt %*% L %*% D_inv_sqrt
  
  eig <- eigen(L_sym, symmetric = TRUE)
  ord <- order(eig$values, decreasing = FALSE)
  # select first k smallest component
  U <- eig$vectors[, ord[seq_len(k)], drop = FALSE]
  
  row_norms <- sqrt(rowSums(U^2))
  row_norms[row_norms == 0] <- 1
  
  U_norm <- U/row_norms
  
  res <- kmeans_robust(U_norm, k = k)
  res
}

#---------------------------
# Wrapper for clusGap()
#---------------------------



preview_shells <- function() {
  demo_dat <- generate_shell_clusters(
    n_shells    = 4,
    k_per_shell = 100,
    max_radius  = 6,
    noise_sd    = 0.1
  )
  
  plt <- plot_ly(
    demo_dat,
    x = ~x, y = ~y, z = ~z,
    color = ~factor(true_shell),
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 2)
  )
  
  if (!dir.exists("figs")) dir.create("figs", recursive = TRUE)
  htmlwidgets::saveWidget(plt, file = "figs/sample_shells_3d.html", selfcontained = TRUE)
}

preview_shells()  # create visualization

#---------------------------
# Simulation settings (Task 2)
#---------------------------
max_radius_values <- 10:0
n_shells_val <- 4
k_per_shell_val <- 100
noise_sd_shell <- 0.1

d_threshold_fixed <- 0.6
spectral_wrapper_for_clusGap <- function(X, k) {
  # X: matrix/data.frame
  # k: assumed number of clusters
  spectral_cluster_once(as.matrix(X), k = k, d_threshold = d_threshold_fixed)
}


task2_results <- tibble(
  max_radius = max_radius_values
) %>%
  mutate(
    est_k = map_int(max_radius, function(Rmax) {
      
      dat_shell <- generate_shell_clusters(
        n_shells = n_shells_val,
        k_per_shell = k_per_shell_val,
        max_radius = Rmax,
        noise_sd = noise_sd_shell
      )
      
      coords <- dat_shell %>%
        select(x, y, z) %>%
        as.matrix()
      
      gap_obj <- clusGap(
        coords,
        FUNcluster = spectral_wrapper_for_clusGap,
        K.max = n_shells_val + 2,
        B = 20
      )
      
      best_k <- maxSE(
        f = gap_obj$Tab[, "gap"],
        SE.f = gap_obj$Tab[, "SE.sim"],
        method = "firstSEmax"
      )
      best_k
    })
  ) %>%
  mutate(true_k = n_shells_val)

if (!dir.exists("results")) dir.create("results", recursive = TRUE)
saveRDS(task2_results, file = "results/task2_gap_results.rds")

p2 <- ggplot(task2_results,
       aes(x = max_radius,
           y = est_k)) +
  geom_line() +
  geom_point() +
  geom_hline(aes(yintercept = true_k),
             linetype = "dashed", color = "black") +
  scale_x_reverse(breaks = unique(task2_results$max_radius)) +
  labs(
    title = "Task 2",
    subtitle = "Dashed # true clusters",
    x = "max_radius",
    y = "Estimated #clusters"
  ) +
  theme_minimal(base_size = 14)


ggsave("results/task2_gap_plot.png", p2, width = 8, height = 5, dpi = 300)
