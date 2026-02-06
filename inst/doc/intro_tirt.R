## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tirt)

## ----sim_irt------------------------------------------------------------------
set.seed(123)
# Define item structure: 1 block of 15 items using 2PL
sim_data <- sim_irt(
  n_people = 500, 
  item_structure = list(
    list(model = "2PL", n_items = 15)
  ),
  theta_mean = 0,
  theta_sd = 1
)

# Access responses and true parameters
head(sim_data$resp)
head(sim_data$true_params)

## ----sim_trt------------------------------------------------------------------
# Define 2 testlets with 5 items each
trt_data <- sim_trt(
  n_people = 500,
  item_structure = list(
    list(model = "2PLT", n_items = 5, testlet_id = "T1"),
    list(model = "2PLT", n_items = 5, testlet_id = "T2")
  )
)

head(trt_data$resp)

## ----binary_irt---------------------------------------------------------------
# Estimate 2PL model
# Note: max_iter set low for vignette speed; use default (100) in practice
fit_2pl <- binary_irt(
  data = sim_data$resp, 
  model = "2PL", 
  method = "EM",
  control = list(max_iter = 20, verbose = FALSE) 
)

# View estimated item parameters
print(head(fit_2pl$item_params))

# View model fit indices
print(fit_2pl$model_fit)

## ----trt_binary---------------------------------------------------------------
# Items 1-5 are Testlet 1, Items 6-10 are Testlet 2
testlet_structure <- list(
  c(1, 2, 3, 4, 5),   # Indices for Testlet 1
  c(6, 7, 8, 9, 10)   # Indices for Testlet 2
)

fit_trt <- trt_binary(
  data = trt_data$resp,
  group = testlet_structure,
  model = "2PLT",
  method = "EM",
  control = list(max_iter = 15, verbose = FALSE)
)

# The output includes item parameters
head(fit_trt$item_params)

# And person parameters including the testlet specific effects (Gamma)
head(fit_trt$person_params)

## ----sim----------------------------------------------------------------------
cat("\n\n--- Generating Simulated Data ---\n")
set.seed(2025)

N <- 1000
# 20 Items Total:
# 1-10: Binary (2PL) - Independent
# 11-15: Polytomous (GRM, 3 cats) - Independent
# 16-17: Binary (2PLT) - Testlet 1
# 18-20: Polytomous (GPCT, 3 cats) - Testlet 2

theta <- rnorm(N, 0, 1)
gamma_1 <- rnorm(N, 0, 0.5) # Testlet 1 effect
gamma_2 <- rnorm(N, 0, 0.6) # Testlet 2 effect

resp_matrix <- matrix(NA, N, 20)
item_names <- paste0("Item_", 1:20)
colnames(resp_matrix) <- item_names

# Define Item Parameters (True Values)
a_true <- runif(20, 0.8, 1.5)
b_true <- seq(-1.5, 1.5, length.out = 20)

# 1. Simulate Binary Independent (1-10)
for(j in 1:10) {
  p <- 1 / (1 + exp(-a_true[j] * (theta - b_true[j])))
  resp_matrix[,j] <- rbinom(N, 1, p)
}

# 2. Simulate Poly Independent (11-15) - GR
for(j in 11:15) {
  b_k <- sort(c(b_true[j] - 0.7, b_true[j] + 0.7))
  p1 <- 1 / (1 + exp(-a_true[j] * (theta - b_k[1])))
  p2 <- 1 / (1 + exp(-a_true[j] * (theta - b_k[2])))
  # GRM Probabilities
  prob_0 <- 1 - p1
  prob_1 <- p1 - p2
  prob_2 <- p2
  # Sample
  r <- runif(N)
  resp_matrix[,j] <- ifelse(r < prob_0, 0, ifelse(r < prob_0 + prob_1, 1, 2))
}

# 3. Simulate Binary Testlet 1 (16-17)
for(j in 16:17) {
  eff_theta <- theta + gamma_1
  p <- 1 / (1 + exp(-a_true[j] * (eff_theta - b_true[j])))
  resp_matrix[,j] <- rbinom(N, 1, p)
}

# 4. Simulate Poly Testlet 2 (18-20) - GPCMT
for(j in 18:20) {
  eff_theta <- theta + gamma_2
  # GPCM steps (centered at b_true)
  b_vec <- c(b_true[j] - 0.5, b_true[j] + 0.5) 
  
  numer <- matrix(0, N, 3)
  numer[,1] <- 0
  numer[,2] <- a_true[j] * (eff_theta - b_vec[1])
  numer[,3] <- a_true[j] * (eff_theta - b_vec[1]) + a_true[j] * (eff_theta - b_vec[2])
  
  max_n <- apply(numer, 1, max)
  exps <- exp(numer - max_n)
  probs <- exps / rowSums(exps)
  
  r <- runif(N)
  resp_matrix[,j] <- apply(cbind(r, probs), 1, function(x) {
    if(x[1] < x[2]) return(0)
    if(x[1] < x[2]+x[3]) return(1)
    return(2)
  })
}

sim_data <- as.data.frame(resp_matrix)
#view the data
print(head(sim_data))

## -----------------------------------------------------------------------------
# Create Item Specification
spec <- data.frame(
  item = item_names,
  model = c(rep("2PL", 10), rep("GRM", 5), rep("2PLT", 2), rep("GPCMT", 3)),
  testlet = c(rep(NA, 15), rep("testlet1", 2), rep("testlet2", 3)),
  stringsAsFactors = FALSE
)
#view the spec
print(spec)

## -----------------------------------------------------------------------------
# Run the Universal Function
results <- irt_trt(sim_data, spec, method = "EM", control = list(max_iter=50, verbose=TRUE))

# Save Results
item_par=results$item_params
fit=results$model_fit
person_par=results$person_params

## ----simulation---------------------------------------------------------------
# --- 1. Simulation Helper Functions ---

# Function to simulate Binary 2PL
sim_2pl <- function(theta, a, b) {
  D <- 1.7
  prob <- 1 / (1 + exp(-D * a * (theta - b)))
  resp <- ifelse(runif(length(theta)) < prob, 1, 0)
  return(resp)
}

# Function to simulate Polytomous (GPCM structure to match fixed_item logic)
sim_poly <- function(theta, a, steps) {
  n_cat <- length(steps) + 1
  n_stud <- length(theta)
  
  # Calculate numerators for categories 0, 1, 2...
  # P(k) proportional to exp( sum(a*(theta - step_j)) )
  
  probs <- matrix(0, nrow = n_stud, ncol = n_cat)
  probs[, 1] <- 1 # Reference category (unnormalized exp is 1 or e^0)
  
  current_sum <- 0
  for (k in 2:n_cat) {
    # Step parameters correspond to boundaries
    current_sum <- current_sum + a * (theta - steps[k-1])
    probs[, k] <- exp(current_sum)
  }
  
  # Normalize
  prob_sum <- rowSums(probs)
  probs_final <- probs / prob_sum
  
  # Random draw
  resp <- numeric(n_stud)
  for (i in 1:n_stud) {
    resp[i] <- sample(0:(n_cat-1), 1, prob = probs_final[i, ])
  }
  return(resp)
}

# --- 2. Generate Parameters for 50 Items ---

n_students <- 5497
true_theta <- rnorm(n_students, mean = 0.4, sd = 1.7)

all_items_meta <- list()
response_matrix <- matrix(NA, nrow = n_students, ncol = 50)
colnames(response_matrix) <- paste0("Item_", 1:50)

# --- Set 1: Known Dichotomous (Items 1-10) ---
# Model: 2PL
for (i in 1:10) {
  a <- runif(1, 0.8, 1.4)
  b <- rnorm(1, 0, 1)
  resp <- sim_2pl(true_theta, a, b)
  response_matrix[, i] <- resp
  all_items_meta[[i]] <- list(item = colnames(response_matrix)[i], type = "Known", model = "2PL", a = a, b = b)
}

# --- Set 2: Known Polytomous (Items 11-18) ---
# Model: Poly (3 categories: 0, 1, 2). Requires 2 steps (d1, d2).
for (i in 11:18) {
  a <- runif(1, 0.7, 1.2)
  d1 <- rnorm(1, -0.5, 0.3)
  d2 <- rnorm(1, 0.5, 0.3) # Steps usually ordered
  resp <- sim_poly(true_theta, a, c(d1, d2))
  response_matrix[, i] <- resp
  all_items_meta[[i]] <- list(item = colnames(response_matrix)[i], type = "Known", model = "Poly", a = a, d1 = d1, d2 = d2)
}

# --- Set 3: Unknown Dichotomous (Items 19-40) ---
# 22 Items
for (i in 19:40) {
  a <- runif(1, 0.8, 1.4)
  b <- rnorm(1, 0, 1)
  resp <- sim_2pl(true_theta, a, b)
  response_matrix[, i] <- resp
  all_items_meta[[i]] <- list(item = colnames(response_matrix)[i], type = "Unknown", model = "2PL", a = a, b = b)
}

# --- Set 4: Unknown Polytomous (Items 41-50) ---
# 10 Items (Varying categories: some 3 cat, some 4 cat)
for (i in 41:50) {
  a <- runif(1, 0.7, 1.2)
  # Mix of 3 categories (2 steps) and 4 categories (3 steps)
  if (i %% 2 == 0) {
    steps <- c(-0.8, 0, 0.8) # 4 cats
  } else {
    steps <- c(-0.5, 0.5)    # 3 cats
  }
  resp <- sim_poly(true_theta, a, steps)
  response_matrix[, i] <- resp
  all_items_meta[[i]] <- list(item = colnames(response_matrix)[i], type = "Unknown", model = "Poly", a = a, steps = steps)
}
#check
response_df <- as.data.frame(response_matrix)
head(response_df[, c(1:3, 11:13, 41:43)]) # Peak at data

## -----------------------------------------------------------------------------
# Create empty dataframe structure
known_params_df <- data.frame(
  item = character(),
  model = character(),
  a = numeric(),
  b = numeric(),
  d1 = numeric(),
  d2 = numeric(),
  stringsAsFactors = FALSE
)

# Populate with Known Items (1-18)
for (i in 1:18) {
  meta <- all_items_meta[[i]]
  
  if (meta$model == "2PL") {
    # Add row for Binary
    known_params_df <- rbind(known_params_df, data.frame(
      item = meta$item,
      model = "2PL",
      a = meta$a,
      b = meta$b,
      d1 = NA, d2 = NA
    ))
  } else {
    # Add row for Poly
    # Note: Using "GPCM" model string as our function uses GPCM math
    known_params_df <- rbind(known_params_df, data.frame(
      item = meta$item,
      model = "GPCM", 
      a = meta$a,
      b = NA,
      d1 = meta$d1,
      d2 = meta$d2
    ))
  }
}

# Display input to verify
print("Known Item Parameters Input:")
print(known_params_df)


## -----------------------------------------------------------------------------
# Run the calibration
# We set model_default = "2PL" (this applies to unknown binary items).
# The function automatically detects poly items in the unknown set based on response categories.
results <- fixed_item(
  response_df = response_df, 
  item_params_df = known_params_df,
  control = list(max_iter = 100)
)

# Check that Items 1-18 are "Fixed" and match inputs.
# Check that Items 19-50 are "Estimated".
# Note how Poly items have step_1, step_2, etc.
item=results$item_params
person=results$person_params
fit=results$model_fit

## ----fix_person---------------------------------------------------------------
# --- Example: With Package Data ---
data("ela1", package = "tirt")

# Select Item Responses (Cols 1-30)
df_real <- ela1[, 1:30]
fixed_theta <- ela1$THETA
fixed_cov <- ela1$COVARIATE
real_res <- fix_person(df = df_real,
                            theta = fixed_theta,
                            model = "2PL",
                            covariate = fixed_cov)
head(real_res)

## ----equate_irt---------------------------------------------------------------
# Create dummy parameters for Form X (Base)
base_params <- data.frame(
  item = c("Item1", "Item2", "Item3", "Item4", "Item5"),
  model = "2PL",
  a = c(1.0, 1.2, 0.9, 1.1, 1.0),
  b = c(-1.0, -0.5, 0.0, 0.5, 1.0),
  stringsAsFactors = FALSE
)

# Create dummy parameters for Form Y (New)
# Suppose Form Y is harder and has different discrimination scaling
new_params <- data.frame(
  item = c("Item1", "Item2", "Item3", "Item6", "Item7"), # Items 1-3 are anchors
  model = "2PL",
  a = c(1.0, 1.2, 0.9, 1.1, 1.0) / 0.9,  # Scale shift
  b = (c(-1.0, -0.5, 0.0, 0.8, 1.2) * 0.9) + 0.2, # Location shift
  stringsAsFactors = FALSE
)

# Perform Equating
linked <- equate_irt(
  base_params = base_params,
  new_params = new_params,
  methods = c("Stocking-Lord", "Mean-Mean")
)

# View Linking Constants (A and B)
print(linked$linking_constants)

# View Transformed Parameters for Form Y (put on Form X scale)
print(linked$transformed_item_params)

