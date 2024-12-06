# 필요한 패키지 설치 및 로드
required_packages <- c("igraph", "Matrix", "MASS", "splines", "ggplot2", "reshape2", "glmnet")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(igraph)
library(Matrix)
library(MASS)       # 다변량 정규분포 샘플링을 위해
library(splines)    # B-spline 기저 함수 생성을 위해
library(ggplot2)    # 시각화를 위해
library(reshape2)   # 데이터 재구성을 위해
library(glmnet)     # 페널티 기반 회귀를 위해

# 1. 데이터 생성 및 네트워크 시각화

# 시간에 따라 변하는 함수 f1(t), f2(t), f3(t) 정의
f1 <- function(t) {
  if (t >= 0 && t <= 0.342) {
    return(5 * (t - 0.5)^2 - 0.125)
  } else if (t > 0.342 && t <= 0.658) {
    return(0)
  } else if (t > 0.658 && t <= 1) {
    return(-5 * (t - 0.5)^2 + 0.125)
  } else {
    stop("t는 [0, 1] 범위 내에 있어야 합니다.")
  }
}

f2 <- function(t) {
  if (t >= 0 && t <= 0.3) {
    return(-3 * t + 0.9)
  } else if (t > 0.3 && t <= 0.7) {
    return(0)
  } else if (t > 0.7 && t <= 1) {
    return(3 * t - 2.1)
  } else {
    stop("t는 [0, 1] 범위 내에 있어야 합니다.")
  }
}

f3 <- function(t) {
  if (t >= 0.3 && t <= 0.7) {
    return(-22.5 * (t - 0.5)^2 + 0.9)
  } else {
    return(0)
  }
}

# 집중 행렬 생성 함수 정의 (대각 우위 보장)
generate_precision_matrix <- function(block_sizes, fk, rho = 0.08, epsilon = 1e-4) {
  p <- sum(block_sizes)  # 전체 노드 수
  A <- matrix(0, nrow = p, ncol = p)  # 0으로 초기화
  
  # 블록의 시작과 끝 인덱스 계산
  block_ends <- cumsum(block_sizes)
  block_starts <- c(1, head(block_ends, -1) + 1)
  
  for (k in 1:length(block_sizes)) {
    start <- block_starts[k]
    end <- block_ends[k]
    size <- block_sizes[k]
    
    # 베르누이 랜덤 변수 생성
    U <- matrix(rbinom(size * size, 1, rho), nrow = size, ncol = size)
    
    # 대각 성분은 0으로 설정 (자기 자신과의 연결은 없음)
    diag(U) <- 0
    
    # 비대각 성분에 fk[k] * U 적용
    A_block <- fk[k] * U
    
    # 행렬 A에 블록 삽입
    A[start:end, start:end] <- A_block
  }
  
  # 대칭 행렬로 만들기 (교환 가능성 유지)
  A <- (A + t(A)) / 2
  
  # 대각 성분 설정 (대각 우위 보장)
  diag(A) <- abs(rowSums(A)) + 1  # 각 행의 절대값 합 + 1
  
  # 양의 정부호성 확인 및 보정
  eigenvalues <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  
  if (min(eigenvalues) <= 0) {
    A <- A + (abs(min(eigenvalues)) + epsilon) * diag(p)
    cat("양의 정부호성 보정을 수행했습니다.\n")
  }
  
  return(A)
}

# 블록별 노드 수 설정
block_sizes <- c(50, 30, 20)  # 블록1: 50노드, 블록2: 30노드, 블록3: 20노드

# 시간 단계 설정
time_steps <- 10

# 샘플 크기 설정
n_samples <- 200

# 전체 노드 수
p <- sum(block_sizes)

# 생성된 데이터를 저장할 배열 초기화 (n_samples x p x time_steps)
data_array <- array(0, dim = c(n_samples, p, time_steps))

# 각 시간 단계에서의 공분산 행렬을 저장할 리스트 초기화
covariance_matrices <- vector("list", time_steps)

# 각 시간 단계에서의 집중 행렬을 저장할 리스트 초기화
precision_matrices <- vector("list", time_steps)

# 시간 점들 계산 (0부터 1까지 균등한 간격으로)
time_points <- seq(0, 1, length.out = time_steps)

# 그래프 생성 및 시각화 반복을 위한 그래픽 파라미터 설정
par(mfrow = c(2, ceiling(time_steps / 2)), mar = c(1, 1, 2, 1))  # 그래프 간 여백 조정

for (t_idx in 1:time_steps) {
  t_real <- time_points[t_idx]
  
  cat("\n시간 단계:", t_idx, "(t_real =", round(t_real, 3), ")\n")
  
  # 각 함수 f1(t), f2(t), f3(t) 계산
  fk_current <- c(f1(t_real), f2(t_real), f3(t_real))
  
  cat("fk_current:", fk_current, "\n")
  
  # 집중 행렬 생성
  A_t <- generate_precision_matrix(block_sizes, fk_current, rho = 0.08)
  
  # A_t의 차원 확인
  cat("A_t dimensions:", dim(A_t), "\n")
  
  # 공분산 행렬 계산 시도
  Sigma_t <- tryCatch({
    solve(A_t)
  }, error = function(e) {
    cat("공분산 행렬을 계산하는 중 에러 발생:", e$message, "\n")
    # 정밀 행렬에 약간의 노이즈 추가하여 양의 정부호성 보장
    A_t_adj <- A_t + 1e-4 * diag(p)
    solve(A_t_adj)
  })
  
  # Sigma_t의 차원 확인
  cat("Sigma_t dimensions:", dim(Sigma_t), "\n")
  
  # 공분산 행렬 대칭화 (수치 오차 방지)
  Sigma_t <- (Sigma_t + t(Sigma_t)) / 2
  
  # 양의 정부호성 확인 및 보정
  eigenvalues <- eigen(Sigma_t, symmetric = TRUE, only.values = TRUE)$values
  min_eigenvalue <- min(eigenvalues)
  
  if (min_eigenvalue <= 0) {
    Sigma_t <- Sigma_t + (abs(min_eigenvalue) + 1e-4) * diag(p)
    cat("양의 정부호성 보정을 수행했습니다.\n")
  }
  
  # 공분산 행렬 및 집중 행렬 저장
  covariance_matrices[[t_idx]] <- Sigma_t
  precision_matrices[[t_idx]] <- A_t
  
  # 다변량 정규분포로부터 데이터 샘플링
  data_samples <- mvrnorm(n = n_samples, mu = rep(0, p), Sigma = Sigma_t)
  
  # 데이터 배열에 저장
  data_array[, , t_idx] <- data_samples
  
  # 집중 행렬을 인접 행렬로 변환 (비대각 성분이 0이 아닌 경우에만 연결)
  adj_matrix_t <- as.matrix(A_t)
  adj_matrix_t[adj_matrix_t != 0] <- 1  # 연결 여부만 표시
  diag(adj_matrix_t) <- 0              # 자기 자신과의 연결 제거
  
  # 그래프 생성
  g_t <- graph_from_adjacency_matrix(adj_matrix_t, mode = "undirected", diag = FALSE)
  
  # 정점 수와 간선 수 확인
  cat("Time_step", t_idx, "has", vcount(g_t), "vertices and", ecount(g_t), "edges.\n")
  
  # 정점 이름을 명시적으로 설정
  V(g_t)$name <- as.character(1:p)
  
  # 블록 정보 추가 (정점 순서대로 할당)
  V(g_t)$block <- rep(1:length(block_sizes), times = block_sizes)
  
  # 블록별 색상 설정
  block_colors <- c("red", "green", "blue")
  
  # 블록 색상이 올바르게 설정되었는지 확인
  if (length(V(g_t)$block) != vcount(g_t)) {
    cat("Error: Length of block assignments (", length(V(g_t)$block), 
        ") does not match number of vertices (", vcount(g_t), ").\n", sep = "")
  } else {
    V(g_t)$color <- block_colors[V(g_t)$block]
  }
  
  # 네트워크 플롯
  plot(g_t,
       vertex.size = 3,
       vertex.label = NA,
       main = paste("시간 단계", t_idx, "\n(t =", round(t_real, 2), ")"),
       layout = layout_with_fr(g_t))
}

# 데이터 구조 확인
cat("\n데이터 배열의 차원:", dim(data_array), "\n")  # n_samples x p x time_steps

# 데이터 배열의 예시 출력
# 첫 번째 샘플, 첫 번째 노드, 시간 단계별 값
cat("첫 번째 샘플, 첫 번째 노드의 시간 단계별 값:\n")
print(data_array[1, 1, ])

# 2. 부분 상관 계수 계산

# 부분 상관 계수를 저장할 리스트 초기화
partial_correlations <- vector("list", time_steps)

for (t_idx in 1:time_steps) {
  Sigma_t <- covariance_matrices[[t_idx]]
  Theta_t <- precision_matrices[[t_idx]]
  
  cat("\nTheta_t for time_step", t_idx, ":\n")
  
  # Theta_t가 NULL인지 확인
  if (is.null(Theta_t)) {
    cat("Theta_t for time_step", t_idx, "is NULL. Skipping.\n")
    next
  }
  
  # Theta_t의 크기 확인
  cat("Theta_t dimensions:", dim(Theta_t), "\n")
  
  if (!all(dim(Theta_t) == c(p, p))) {
    cat("Theta_t for time_step", t_idx, "does not match expected dimensions. Skipping.\n")
    next
  }
  
  cat("Processing time_step", t_idx, "\n")
  
  # 부분 상관 계수 계산
  rho_matrix_t <- matrix(0, nrow = p, ncol = p)
  
  for (i in 1:(p-1)) {
    for (j in (i+1):p) {
      rho_ij <- -Theta_t[i, j] / sqrt(Theta_t[i, i] * Theta_t[j, j])
      rho_matrix_t[i, j] <- rho_ij
      rho_matrix_t[j, i] <- rho_ij
    }
  }
  
  partial_correlations[[t_idx]] <- rho_matrix_t
}

# 부분 상관 계수 존재 여부 확인
cat("\npartial_correlations 리스트의 NULL 여부:\n")
print(sapply(partial_correlations, is.null))  # Should all be FALSE

# 3. B-spline 기저 함수 생성

# B-spline 기저 함수 생성 (계수 중첩을 반영하도록 수정)
num_knots <- 5  # 내부 노드 수를 5로 설정
degree <- 3     # 3차 B-spline

# 시간 점들 계산 (0부터 1까지 균등한 간격으로)
time_points <- seq(0, 1, length.out = time_steps)

# 노드 설정 (경계 노드 포함)
knots <- c(rep(time_points[1], degree), 
           quantile(time_points, probs = seq(0, 1, length.out = num_knots + 2))[-c(1, num_knots + 2)], 
           rep(time_points[length(time_points)], degree))

# splineDesign 함수 호출 시 outer.ok = TRUE 설정
bs_basis <- splineDesign(knots = knots, x = time_points, ord = degree + 1, outer.ok = TRUE)

# 기저 함수의 차원 확인
cat("\nbs_basis dimensions:", dim(bs_basis), "\n")  # time_steps x (num_knots + degree)

# B-spline 기저 함수 데이터를 데이터 프레임으로 변환
spline_data <- data.frame(Time = time_points)
for (h in 1:ncol(bs_basis)) {
  spline_data[[paste0("B", h)]] <- bs_basis[, h]
}

# 데이터 프레임을 길게 변환 (melt)
if (exists("spline_data")) {
  spline_data_long <- melt(spline_data, id.vars = "Time", variable.name = "Spline", value.name = "Value")
} else {
  stop("Error: 'spline_data' 객체가 생성되지 않았습니다.")
}

# `spline_data_long` 객체 확인
if (!exists("spline_data_long")) {
  stop("Error: 'spline_data_long' 객체가 생성되지 않았습니다.")
}

cat("\nspline_data_long의 첫 몇 행:\n")
print(head(spline_data_long))

# 4. 그룹 SCAD 패널티 및 최적화 알고리즘 구현

# 선택된 변수 쌍 선택 (예: 5개의 변수 쌍)
selected_pairs <- list(
  c(1, 2),
  c(1, 3),
  c(2, 4),
  c(5, 6),
  c(7, 8)
)

# 그룹 구조 정의: 각 스플라인 기저 함수는 degree+1개의 인접한 계수와 연관됨
group_size <- degree + 1
groups <- list()
num_basis <- ncol(bs_basis)
for (h in 1:num_basis) {
  # 각 기저 함수에 대해 그룹 정의 (중첩)
  # h번째 기저 함수는 h, h+1, ..., h+degree 계수와 관련됨
  group_indices <- h:(h + degree)
  # 계수가 범위를 벗어나지 않도록 조정
  group_indices <- group_indices[group_indices <= num_basis]
  groups[[h]] <- group_indices
}

# 예시: 그룹 구조 확인
cat("\n그룹 구조의 첫 5개 그룹:\n")
print(groups[1:5])

# SCAD 페널티 함수 정의
scad <- function(beta, lambda, gamma = 3.7) {
  abs_beta <- abs(beta)
  penalty <- ifelse(abs_beta <= lambda, lambda * abs_beta,
                    ifelse(abs_beta <= gamma * lambda, 
                           (-abs_beta^2 + 2 * gamma * lambda * abs_beta - lambda^2) / (2 * (gamma - 1)),
                           (gamma + 1) * lambda^2 / 2))
  return(penalty)
}

# SCAD 페널티 함수 정의 (그룹 단위)
scad_group <- function(beta_group, lambda, gamma = 3.7) {
  norm_beta <- sqrt(sum(beta_group^2))
  if (norm_beta <= lambda) {
    return(lambda * norm_beta)
  } else if (norm_beta <= gamma * lambda) {
    return((-norm_beta^2 + 2 * gamma * lambda * norm_beta - lambda^2) / (2 * (gamma - 1)))
  } else {
    return((gamma + 1) * lambda^2 / 2)
  }
}

# 전체 손실 함수 계산
calculate_piecewise_loss <- function(beta, sigma, time_points, y, w, b_spline_basis, groups, lambda, gamma = 3.7) {
  n <- dim(y)[1]       # 샘플 수
  p <- dim(y)[2]       # 변수 수
  time_steps <- length(time_points)
  Jn <- dim(b_spline_basis)[2]  # B-spline 계수의 수
  
  # 잔차 제곱합 계산
  rss <- 0
  for (k in 1:n) {
    for (i in 1:p) {
      for (u in 1:time_steps) {
        tku <- time_points[u]
        prediction <- 0
        for (j in 1:p) {
          if (i != j) {
            Bh <- b_spline_basis[u, ]
            beta_ij <- beta[i, j, ]
            prediction <- prediction +
              sum(beta_ij * Bh) *
              sqrt(sigma[j, j, u] / sigma[i, i, u]) * y[k, j, u]
          }
        }
        residual <- y[k, i, u] - prediction
        rss <- rss + w[k, i, u] * residual^2
      }
    }
  }
  
  # SCAD 페널티 계산 (그룹 단위)
  penalty <- 0
  for (g in 1:length(groups)) {
    group_indices <- groups[[g]]
    beta_group <- beta[group_indices]
    penalty <- penalty + scad_group(beta_group, lambda, gamma)
  }
  
  # 최종 손실 함수 값
  loss <- (1 / (2 * n * time_steps)) * rss + penalty
  return(loss)
}

# Proximal Gradient Method을 이용한 최적화 함수 정의
optimize_with_proximal_gradient <- function(y, X, groups, lambda, gamma = 3.7, max_iter = 1000, tol = 1e-6) {
  J <- ncol(X)
  beta <- rep(0, J)  # 초기값
  step_size <- 1 / (max(eigen(t(X) %*% X)$values) + 1)  # 스텝 사이즈 설정
  
  for (iter in 1:max_iter) {
    # 그라디언트 계산
    residual <- y - X %*% beta
    grad <- -t(X) %*% residual
    
    # 그라디언트 스텝
    beta_new <- beta - step_size * grad
    
    # 프로시멀 연산 (그룹 SCAD 패널티)
    for (g in 1:length(groups)) {
      idx <- groups[[g]]
      beta_subset <- beta_new[idx]
      
      # SCAD 프로시멀 함수 적용
      norm_beta <- sqrt(sum(beta_subset^2))
      if (norm_beta == 0) {
        beta_new[idx] <- 0
      } else {
        # 그룹 SCAD 프로시멀 연산
        shrinkage <- max(0, 1 - scad(norm_beta, lambda, gamma) / norm_beta)
        beta_new[idx] <- shrinkage * beta_subset
      }
    }
    
    # 수렴 여부 확인
    if (sum(abs(beta_new - beta)) < tol) {
      cat("Converged at iteration", iter, "\n")
      break
    }
    
    beta <- beta_new
    if (iter == max_iter) {
      cat("Reached maximum iterations\n")
    }
  }
  
  return(beta)
}

# 선택된 변수 쌍에 대해 B-spline 모델을 피팅하고 계수 저장
spline_coefficients <- list()

for (pair in selected_pairs) {
  i <- pair[1]
  j <- pair[2]
  
  cat("\n선택된 변수 쌍:", i, "-", j, "\n")
  
  # 각 시간 단계에서의 부분 상관 계수 추출
  rho_ij_t <- sapply(partial_correlations, function(rho_mat) {
    if (!is.null(rho_mat)) rho_mat[i, j] else NA
  })
  
  if (all(is.na(rho_ij_t))) {
    cat("Warning: 부분 상관 계수에 NA 값만 포함되어 있습니다. 해당 변수 쌍을 건너뜁니다.\n")
    spline_coefficients[[paste(i, j, sep = "_")]] <- rep(NA, ncol(bs_basis))
    next
  }
  
  # 결측치 제거
  valid_idx <- which(!is.na(rho_ij_t))
  y_fit <- rho_ij_t[valid_idx]
  X_fit <- bs_basis[valid_idx, ]
  
  # SCAD 페널티 적용을 위한 최적화
  lambda <- 0.1  # 페널티 매개변수 (조정 가능)
  gamma <- 3.7    # SCAD 매개변수
  
  # 프로시멀 그래디언트 알고리즘을 사용한 최적화
  beta_est <- optimize_with_proximal_gradient(y = y_fit, X = X_fit, groups = groups, 
                                              lambda = lambda, gamma = gamma, 
                                              max_iter = 1000, tol = 1e-6)
  
  spline_coefficients[[paste(i, j, sep = "_")]] <- beta_est
}

# 선택된 변수 쌍의 부분 상관 계수 및 스플라인 기저 함수 데이터를 데이터 프레임으로 변환
plot_data <- data.frame(
  Time = rep(time_points, length(selected_pairs)),
  Pair = rep(sapply(selected_pairs, function(x) paste(x[1], x[2], sep = "_")), each = time_steps),
  PartialCorrelation = unlist(lapply(selected_pairs, function(pair) {
    i <- pair[1]
    j <- pair[2]
    sapply(partial_correlations, function(rho_mat) if (!is.null(rho_mat)) rho_mat[i, j] else NA)
  })),
  SplineValue = unlist(lapply(selected_pairs, function(pair) {
    pair_name <- paste(pair[1], pair[2], sep = "_")
    if (!is.null(spline_coefficients[[pair_name]]) && !any(is.na(spline_coefficients[[pair_name]]))) {
      beta_est <- spline_coefficients[[pair_name]]
      # 계수 중첩을 반영하여 스플라인 값을 계산
      fitted_values <- bs_basis %*% beta_est
      return(fitted_values)
    } else {
      return(rep(NA, time_steps))
    }
  }))
)

# `plot_data` 상태 확인
cat("\nplot_data 상태:\n")
print(head(plot_data))

# `plot_data` NA 값 확인
if (any(is.na(plot_data$SplineValue))) {
  cat("Warning: SplineValue에 NA 값이 포함되어 있습니다. 시각화 시 일부 데이터가 제외될 수 있습니다.\n")
}

# `Pair` 변수를 factor로 변환하여 이산형 변수로 인식하게 함
plot_data$Pair <- as.factor(plot_data$Pair)

# ggplot2를 사용하여 부분 상관 계수 시각화
ggplot(data = plot_data, aes(x = Time, y = PartialCorrelation, color = Pair)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("blue", "green", "red", "purple", "orange"),
    name = "Variable Pairs"
  ) +
  labs(title = "시간에 따른 부분 상관 계수",
       x = "시간 (t)",
       y = "부분 상관 계수") +
  theme_minimal()

# spline_coefficients의 첫 몇 개 확인
print(spline_coefficients[1:2])

# ggplot2를 사용하여 부분 상관 계수와 B-spline 기저 함수 시각화
ggplot() +
  geom_line(data = spline_data_long, aes(x = Time, y = Value, group = Spline), 
            color = "grey", size = 0.5, alpha = 0.3) +
  geom_line(data = plot_data, aes(x = Time, y = SplineValue, color = Pair), 
            size = 1) +
  geom_point(data = plot_data, aes(x = Time, y = PartialCorrelation, color = Pair), 
             size = 1, alpha = 0.6) +
  scale_color_manual(
    values = c("blue", "green", "red", "purple", "orange"),
    name = "Variable Pairs"
  ) +
  labs(title = "시간에 따른 부분 상관 계수 및 B-spline 기저 함수",
       x = "시간 (t)",
       y = "부분 상관 계수 / 스플라인 값") +
  theme_minimal()
#회색 스플라인: 스플라인 기저 함수의 중첩성을 시각적으로 확인할 수 있습니다.
#희소성: 많은 변수 쌍이 회색 스플라인만 표시되고, 주요 변수 쌍만 색상별로 강조되었다면 희소성이 잘 유도된 것입니다

