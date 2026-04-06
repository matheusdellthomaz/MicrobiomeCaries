n_cores <- detectCores() - 2
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)

setwd("C:/Users/Home/Documents")

# Data
data <- as.data.frame(fread("subset_db_rf.csv"))
data <- data %>%
  mutate(STATUS = factor(STATUS, levels = c("Health", "Caries")),
         response = factor(if_else(STATUS == "Caries", 1, 0), 
                           levels = c(0, 1)))
data$response <- as.factor(data$response)

asvs <- data %>% select(-STATUS, -response)
# asvs <- data %>% select(starts_with("ASV"))
responses <- data %>% select(response)

# Filter ASV f > 10%
asv_filter <- asvs[, colSums(asvs > 0) >= 0.1 * nrow(asvs)]
asvres <- cbind(asv_filter,responses)

#Feature importance and selection
set.seed(111)

folds <- vfold_cv(asvres, v = 5, strata = response)

glmnet_fs <- map_df(folds$splits, function(split) {
  
  train_fold <- analysis(split)
  
  X <- model.matrix(response ~ . -1, train_fold)
  y <- train_fold$response
  
  fit <- cv.glmnet(X, y, family = "binomial", alpha = 1)
  
  coef_sel <- coef(fit, s = "lambda.min")
  
  tibble(
    feature = rownames(coef_sel),
    selected = as.numeric(coef_sel[,1] != 0)
  ) %>%
    filter(feature != "(Intercept)")
})

# Frequency selection
glmnet_freq <- glmnet_fs %>%
  group_by(feature) %>%
  summarise(freq = mean(selected)) %>%
  arrange(desc(freq))

# Threshold (> 60%)
features_glmnet_cv <- glmnet_freq %>%
  filter(freq >= 0.6)

# Random Forest Selection
rf_fs <- map_df(folds$splits, function(split) {
  
  train_fold <- analysis(split)
  
  rf <- ranger(
    response ~ ., 
    data = train_fold, 
    importance = "permutation"
  )
  
  tibble(
    feature = names(rf$variable.importance),
    importance = rf$variable.importance
  )
})

rf_importance_cv <- rf_fs %>%
  group_by(feature) %>%
  summarise(mean_importance = mean(importance)) %>%
  arrange(desc(mean_importance))

# Threshold
features_rf_cv <- rf_importance_cv %>%
  filter(mean_importance > 5e-4)

# Boruta
boruta_fs <- map_df(folds$splits, function(split) {
  
  train_fold <- analysis(split)
  
  boruta_model <- Boruta(response ~ ., data = train_fold, doTrace = 0)
  stats <- attStats(boruta_model)
  
  tibble(
    feature = rownames(stats),
    decision = stats$decision
  )
})

boruta_summary <- boruta_fs %>%
  group_by(feature) %>%
  summarise(
    prop_confirmed = mean(decision == "Confirmed")
  ) %>%
  arrange(desc(prop_confirmed))

# Threshold 
features_boruta_cv <- boruta_summary %>%
  filter(prop_confirmed >= 0.6)

# Recreating steps from Combinatory analysis of features for plotting and statistics
# Select combination or subset or combination for input in test_combinations function
asvres <- asvres %>% select("ASV765", "ASV393","ASV1445","ASV273","ASV100","response")

# Train/test split
set.seed(123)
split <- initial_split(asvres, 
                       strata = response, 
                       prop= 0.75)
train  <- training(split)
test  <- testing(split)

# Recipe
final_rec <- recipe(response ~ ., data = train) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_zv(all_numeric_predictors())

# Select workflows and fits for each combination
final_wf <- results$workflows$size_5$ASV669_ASV2294_ASV273_ASV765_ASV652glmnet

final_wff <- results$fits$size_5$ASV669_ASV2294_ASV273_ASV765_ASV652glmnet

# Crossvalidation
set.seed(456)

folds <- vfold_cv(train, v = 10, strata = response)

cv_res <- fit_resamples(final_wf, folds, control = control_resamples(verbose=T, save_pred = T))
cv_preds <- collect_predictions(cv_res)
cv_metrics <- calculate_class_metrics(cv_preds) 

# Test Metrics
test_preds <- predict(final_wff, test, type = "prob") %>%
  bind_cols(predict(final_wff, test, type = "class"), test)
test_metrics <- calculate_class_metrics(test_preds)

cv_preds$response <- as.factor(cv_preds$response)
test_preds$response <- as.factor(test_preds$response)

# ROC_AUC metrics CV
roc_curve <- roc_curve(cv_preds, response, .pred_1, event_level = "second")

# AUC
auc_value <- roc_auc(cv_preds, response, .pred_1, event_level = "second")

# Plot
ggplot(roc_curve, aes(x = 1 - specificity, y = sensitivity)) +
  geom_path(color = "blue", size = 1) +
  geom_abline(lty = 2, color = "gray") +
  coord_equal() +
  labs(title = paste("CV10F GLMNet ROC Curve - AUC =", round(auc_value$.estimate, 3)),
       x = "1 - Specificity",
       y = "Sensibility") +
  theme_minimal()

# Confusion Matrix
conf_mat <- test_preds %>% 
  conf_mat(truth = response, estimate = .pred_class)
conf_mat_plot <- conf_mat %>% 
  autoplot(type = "heatmap")
conf_mat_plot

#PCA

pca_data <- asvres
pca_data <- pca_data %>%
  mutate(across(starts_with("ASV"), ~ as.numeric(scale(.))))

pca_result <- prcomp(pca_data %>% select(starts_with("ASV")), 
                     scale. = F, center = F)
summary(pca_result)

autoplot(
  pca_result,
  data = pca_data,
  colour = "response",
  loadings = FALSE,
  frame = TRUE
) +
  scale_color_manual(values = c(
    "0" = "blue",
    "1" = "red"
  )) +
  labs(
    title = "ASVs PCA",
    color = "Status"
  ) +
  theme_minimal()

scores <- as.data.frame(pca_result$x)
scores$response <- pca_data$response

# Loadings
loadings <- as.data.frame(pca_result$rotation[, 1:2]) %>%
  rownames_to_column("ASV")

loadings %>% 
  arrange(desc(abs(PC1))) %>% 
  head()

# % Contribution
contrib <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Loadings contribution
ggplot(loadings, aes(x = PC1, y = PC2, label = ASV)) +
  geom_segment(aes(xend = 0, yend = 0), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "blue") +
  geom_text(size = 5, vjust = 1.5, color = "darkred") +
  geom_point(size = 3, color = "red") +
  xlab(paste0("PC1 (", round(contrib[1], 1), "%)")) +  
  ylab(paste0("PC2 (", round(contrib[2], 1), "%)")) +  
  ggtitle("ASVs contribution for PCA") +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")

# Scree plot
plot(pca_result, type = "l", main = "Scree Plot")

# Correlation between ASVs abundance and Status (0=health, 1=Carie)
# GLM logistic
glm <- pca_data %>%
  select(starts_with("ASV")) %>%
  map_df(~ tidy(glm(pca_data$response ~ .x, family = binomial)),
         .id = "ASV")

# Pearson point-bisserial
y <- as.numeric(pca_data$response) - 1  
correlation <- pca_data %>%
  select(starts_with("ASV")) %>%
  map(~ cor.test(.x, y))

cor_coef <- map_dbl(correlation, "estimate")
p_value  <- map_dbl(correlation, "p.value")

# ANOVA
asv_anova <- pca_data %>%
  select(response, starts_with("ASV")) %>%
  pivot_longer(-response, names_to = "ASV", values_to = "Abundance") %>%
  group_by(ASV) %>%
  summarise(p_valor = wilcox.test(Abundance ~ response)$p.value) %>%
  mutate(p_ajustado = p.adjust(p_valor, method = "fdr"))

#PLS-DA

X <- asvres %>% select(-response)
Y <- asvres$response  %>% as.factor()

PLSDA_model <- mixOmics::plsda(X = X, 
                               Y = Y, 
                               ncomp = 2, 
                               scale = TRUE)

# Loadings
PLSDA_loadings <- PLSDA_model$loadings$X

# scores
PLSDA_scores <- data.frame(PLSDA_model$variates$X) %>%
  rownames_to_column("sample_ID") %>%
  mutate(response = Y)

# variance
var_explained <- PLSDA_model$prop_expl_var$X * 100

# Score plot
PLSDA_plot <- ggplot(PLSDA_scores, aes(x = comp1, y = comp2, color = response)) +
  geom_point(alpha = 0.6, size = 3) +
  stat_ellipse(level = 0.95, linetype = 2, linewidth = 0.8) +
  geom_point(data = PLSDA_scores %>% 
               group_by(response) %>% 
               summarise(comp1 = mean(comp1), comp2 = mean(comp2)),
             aes(x = comp1, y = comp2, color = response), 
             size = 5, shape = 8) +
  labs(title = "PLS-DA - ASVs",
       x = paste("Component 1 (", round(var_explained[1], 1), "%)", sep = ""),
       y = paste("Component 2 (", round(var_explained[2], 1), "%)", sep = ""),
       color = "Group") +
  scale_color_manual(values = c("0" = "#6a3d9a", "1" = "#ff7f00"),
                     labels = c("0" = "Healthy", "1" = "Caries")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "right") +
  coord_fixed()
print(PLSDA_plot)

# Loadings DF
loadings_df <- data.frame(
  ASV = rownames(PLSDA_loadings),
  comp1 = PLSDA_loadings[, 1],
  comp2 = PLSDA_loadings[, 2]
)

# loadings
loadings_plot <- ggplot(loadings_df, aes(x = comp1, y = comp2)) +
  geom_segment(aes(x = 0, y = 0, xend = comp1, yend = comp2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "red", alpha = 0.8) +
  geom_text(aes(label = ASV), vjust = -0.5, size = 4, fontface = "bold") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(title = "ASVs Loadings - Direcionalidad",
       x = paste("Loading Component 1 (", round(var_explained[1], 1), "%)", sep = ""),
       y = paste("Loading Component 2 (", round(var_explained[2], 1), "%)", sep = "")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  xlim(-1, 1) + ylim(-1, 1) +
  coord_fixed()

print(loadings_plot)

library(ggrepel)

combined_plot <- ggplot() +
  # Scores
  geom_point(data = PLSDA_scores, 
             aes(x = comp1, y = comp2, color = response), 
             alpha = 0.6, size = 3) +
  # ASVs Loadings
  geom_segment(data = loadings_df,
               aes(x = 0, y = 0, xend = comp1 * 2.5, yend = comp2 * 2.5),  
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "darkred", size = 1.2) +  
  
  geom_text_repel(data = loadings_df,
                  aes(x = comp1 * 2.6, y = comp2 * 2.6, label = ASV),
                  color = "darkred", 
                  size = 5,
                  min.segment.length = 0,  
                  segment.color = "grey50",  
                  segment.size = 0.5,  
                  box.padding = 0.5,  
                  point.padding = 0.5,  
                  force = 1,  
                  max.overlaps = Inf) +  
  # Centroids
  geom_point(data = PLSDA_scores %>% 
               group_by(response) %>% 
               summarise(comp1 = mean(comp1), comp2 = mean(comp2)),
             aes(x = comp1, y = comp2, color = response), 
             size = 6, shape = 8) + 
  labs(title = "PLS-DA: ASVs Scores and Loadings for RF combination",
       x = paste("Component 1 (", round(var_explained[1], 1), "%)", sep = ""),
       y = paste("Component 2 (", round(var_explained[2], 1), "%)", sep = ""),
       color = "Group") +
  scale_color_manual(values = c("0" = "#6a3d9a", "1" = "#ff7f00"),
                     labels = c("0" = "Healthy", "1" = "Caries")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 12))

print(combined_plot)

# VIP
vip_scores <- mixOmics::vip(PLSDA_model)
print(vip_scores)

# Plotar VIP scores
vip_df <- data.frame(ASV = rownames(vip_scores), VIP = vip_scores[, 1])
vip_plot <- ggplot(vip_df, aes(x = reorder(ASV, VIP), y = VIP)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "VIP Scores - ASVs Importance",
       x = "ASV", y = "VIP Score") +
  theme_minimal()

print(vip_plot)

