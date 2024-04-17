## code to prepare the test datasets goes here
library(tidyverse)

# Import dummy claims dataset from SynthETIC
test_claims_object <- SynthETIC::test_claims_object
# Convert claims dataset into triangular format
test_claims <- SynthETIC::claim_output(
    frequency_vector = test_claims_object$frequency_vector,
    payment_time_list = test_claims_object$payment_time_list,
    payment_size_list = test_claims_object$payment_size_list
)

# Convert claims dataset into required format for ADLP package:
#   - `origin` and `dev` columns as columns 1 and 2 respectively for AP and DP
test_claims_dataset <- test_claims %>%
    as.data.frame() %>%
    tibble::rownames_to_column(., "origin") %>%
    tidyr::pivot_longer(., cols = -c("origin"),
                 names_to = "dev",
                 values_to = "claims") %>%
    dplyr::mutate(origin = as.numeric(substring(origin, 3)),
           dev = as.numeric(substring(dev, 3)),
           calendar = origin + dev) %>%
    dplyr::select(
        origin, dev, calendar, claims
    ) %>%
    as.data.frame()

usethis::use_data(test_claims_dataset, overwrite = T)

## Code to create test_adlp_component
train_val <- train_val_split_method1(
    df = test_claims_dataset,
    tri.size = 40,
    val_ratio = 0.3,
    test = TRUE
)

train_data <- train_val$train
valid_data <- train_val$valid
insample_data <- rbind(train_data, valid_data)
base_model1 <- glm(formula = claims~factor(dev),
                family=gaussian(link = "identity"), data=train_data)

base_model1_full <- update(base_model1, data = insample_data)

dens_normal <- function(y, model, newdata){
    pred_model <- predict(model, newdata=newdata, type="response", se.fit=TRUE)
    mu <- pred_model$fit
    sigma <- pred_model$residual.scale
    return(dnorm(x=y, mean=mu, sd=sigma))
}
cdf_normal<-function(y, model, newdata){
    pred_model <- predict(model, newdata=newdata, type="response", se.fit=TRUE)
    mu <- pred_model$fit
    sigma <- pred_model$residual.scale
    return(pnorm(q=y, mean=mu, sd=sigma))
}
mu_normal<-function(model, newdata){
    mu <- predict(model, newdata=newdata, type="response")
    mu <- pmax(mu, 0)
    return(mu)
}
sim_normal<-function(model, newdata){
    pred_model <- predict(model, newdata=newdata, type="response", se.fit=TRUE)
    mu <- pred_model$fit
    sigma <- pred_model$residual.scale
    sim <- rnorm(length(mu), mean=mu, sd=sigma)
    sim <- pmax(sim, 0)
    return(sim)
}

test_adlp_component = adlp_component(
    model_train = base_model1,
    model_full = base_model1_full,
    calc_dens = dens_normal,
    calc_cdf = cdf_normal,
    calc_mu = mu_normal,
    sim_fun = sim_normal
)

usethis::use_data(test_adlp_component, overwrite = T)
