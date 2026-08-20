## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tirt)

## ----calibrate----------------------------------------------------------------
set.seed(2025)
sim <- sim_irt(
  n_people = 600,
  item_structure = list(list(model = "2PL", n_items = 12))
)

fit <- binary_irt(sim$resp, model = "2PL", method = "EM",
                  control = list(max_iter = 20, verbose = FALSE))
head(fit$item_params)

## ----information--------------------------------------------------------------
theta_grid <- seq(-3, 3, by = 0.5)

# Item information (items in rows, theta in columns)
info <- item_info(fit$item_params, theta = theta_grid)
round(info[1:3, ], 3)

# Test information function and conditional SEM
tif <- test_info(fit$item_params, theta = theta_grid)
tif

# Where does the test measure most precisely?
tif$theta[which.max(tif$test_info)]

## ----score_table--------------------------------------------------------------
# EAP conversion table (0 to 12 correct)
score_table(fit$item_params, method = "EAP")

# Maximum-likelihood conversion
score_table(fit$item_params, method = "MLE")

## ----person_fit---------------------------------------------------------------
pf <- person_fit(sim$resp, fit$item_params, fit$person_params)
head(pf)

# Number of examinees flagged as potentially misfitting
sum(pf$flag, na.rm = TRUE)

## ----item_fit-----------------------------------------------------------------
item_fit(sim$resp, fit$item_params)

## ----ld_stats-----------------------------------------------------------------
q3 <- ld_stats(sim$resp, fit$item_params)
round(q3[1:5, 1:5], 3)

# Largest absolute residual correlation
attr(q3, "max_abs_q3")

## ----dif----------------------------------------------------------------------
resp_dif <- sim$resp
grp <- rep(c("Reference", "Focal"), each = 300)
flip <- grp == "Focal" & resp_dif[[3]] == 1
resp_dif[[3]][flip] <- rbinom(sum(flip), 1, 0.55)

dif(resp_dif, group = grp)[, c("item", "MH_delta", "ETS_class", "LR_p")]

## ----reliability--------------------------------------------------------------
reliability(person_params = fit$person_params,
            data = sim$resp,
            item_params = fit$item_params)

## ----tcc----------------------------------------------------------------------
curves <- tcc(fit$item_params, theta = seq(-3, 3, by = 1))
curves$test_curve

## ----sim_mirt-----------------------------------------------------------------
# Two correlated dimensions, simple structure
mdat <- sim_mirt(
  n_people = 400,
  dimension = 2,
  Sigma = matrix(c(1, 0.4, 0.4, 1), 2, 2),
  item_structure = list(
    list(model = "M2PL", n_items = 6, dims = 1),
    list(model = "M2PL", n_items = 6, dims = 2)
  )
)
head(mdat$true_params)

## ----sim_tirt-----------------------------------------------------------------
# Independent items plus two testlets
tdat <- sim_tirt(
  n_people = 400,
  item_structure = list(
    list(model = "2PL",  n_items = 6),
    list(model = "2PLT", n_items = 4, testlet_id = "P1", testlet_var = 0.6),
    list(model = "GPCT", n_items = 3, categories = 3, testlet_id = "P2")
  )
)
tdat$true_item_params[, c("item_id", "model", "testlet")]

## ----mixture------------------------------------------------------------------
set.seed(11)
N <- 300; J <- 8
b1 <- seq(-1.5, 1.5, length.out = J); b2 <- rev(b1)
theta <- rnorm(N); cls <- rep(1:2, each = N / 2)
rmat <- matrix(0, N, J)
for (i in 1:N) {
  b <- if (cls[i] == 1) b1 else b2
  rmat[i, ] <- rbinom(J, 1, 1 / (1 + exp(-(theta[i] - b))))
}
mdf <- as.data.frame(rmat); names(mdf) <- paste0("I", 1:J)

mix <- mixture_irt(mdf, n_class = 2, model = "Rasch",
                   control = list(max_iter = 40, verbose = FALSE))
mix$class_params
mix$model_fit

