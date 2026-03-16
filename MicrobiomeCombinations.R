stopImplicitCluster()
plan(sequential)
gc()
n_cores <- detectCores() - 2
cl <- makePSOCKcluster(n_cores)
registerDoParallel(cl)
calculate_class_metrics <- function(pred_df) {
  
  # Metrics
  metrics <- pred_df %>%
    yardstick::metrics(truth = response, estimate = .pred_class) %>%
    bind_rows(
      yardstick::sens(pred_df, truth = response, estimate = .pred_class, event_level = "second"),
      yardstick::spec(pred_df, truth = response, estimate = .pred_class, event_level = "second"),
      yardstick::ppv(pred_df, truth = response, estimate = .pred_class, event_level = "second"),
      yardstick::npv(pred_df, truth = response, estimate = .pred_class, event_level = "second"),
      yardstick::f_meas(pred_df, truth = response, estimate = .pred_class, event_level = "second"),
      yardstick::roc_auc(pred_df, truth = response, .pred_1, event_level = "second")
    )
  
  return(metrics)
}

test_combinations <- function(data, time_predictors, fixed_predictors = NULL, response,
                              external_dfs = list(), combo_sizes = c(2, 3),
                              include_fixed_predictors = TRUE,
                              include_ratio_predictors = FALSE, 
                              include_diff_predictors = FALSE,
                              algorithms = c("rf", "xgb", "svm", "glmnet"),
                              return_fits = FALSE,
                              return_workflows = FALSE) {
  
  # Lists for storing results
  all_results <- list()
  all_best_params <- list()
  all_tabs <- list()
  all_fits <- list()
  all_workflows <- list()
  
  
  # Control settings
  ctrl <- control_resamples(save_pred = TRUE, save_workflow = TRUE)
  race_ctrl <- control_race(save_pred = TRUE)
  
  # Verify algorithms
  valid_algs <- c("rf", "xgb", "svm", "glmnet")
  algorithms <- intersect(algorithms, valid_algs)
  
  # Possible derived predictors
  create_derived_predictors <- function(df, predictors, 
                                        include_ratio, include_diff) {
    new_df <- df
    ratio_names <- c()
    diff_names <- c()
    
    if (include_ratio && length(predictors) >= 2) {
      ratios <- combn(predictors, 2, simplify = FALSE) %>%
        map(function(pair) {
          name <- paste0(pair[2], "_div_", pair[1])
          new_df <<- new_df %>%
            mutate(!!name := ifelse(!!sym(pair[1]) == 0, NA, !!sym(pair[2]) / !!sym(pair[1])))
          name
        })
      ratio_names <- unlist(ratios)
    }
    
    if (include_diff && length(predictors) >= 2) {
      diffs <- combn(predictors, 2, simplify = FALSE) %>%
        map(function(pair) {
          name <- paste0(pair[2], "_sub_", pair[1])
          new_df <<- new_df %>%
            mutate(!!name := !!sym(pair[2]) - !!sym(pair[1]))
          name
        })
      diff_names <- unlist(diffs)
    }
    
    list(data = new_df, ratio_names = ratio_names, diff_names = diff_names)
  }
  
  # Algorithms specs for tuning
  model_specs <- list(
    rf = rand_forest(
      mtry = tune(),
      min_n = tune(),
      trees = 1000
    ) %>% 
      set_engine("ranger", importance = "impurity") %>% 
      set_mode("classification"),
    
    xgb = boost_tree(
      mtry = tune(),
      trees = tune(),
      min_n = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      loss_reduction = tune(),
      sample_size = tune()
    ) %>% 
      set_engine("xgboost", nthread = 10) %>% 
      set_mode("classification"),
    
    svm = svm_rbf(
      cost = tune(),
      rbf_sigma = tune()
    ) %>% 
      set_engine("kernlab") %>% 
      set_mode("classification"),
    
    glmnet = logistic_reg(
      penalty = tune(),
      mixture = tune()
    ) %>%
      set_engine("glmnet") %>% 
      set_mode("classification")
  )
  
  # Grid search ranges
  grids <- list(
    rf = grid_regular(
      mtry(range = c(10, 50)),
      min_n(range = c(2, 20)),
      levels = 5
    ),
    
    xgb = grid_space_filling(
      mtry(range = c(10, 50)),
      trees(range = c(50, 1000)),
      min_n(range = c(2, 20)),
      tree_depth(range = c(3, 10)),
      learn_rate(range = c(-3, -1)),
      loss_reduction(range = c(-5, -1)),
      sample_size(range = c(0,1)),
      size = 20
    ),
    
    svm = grid_space_filling(
      cost(range = c(-5, 5)),
      rbf_sigma(range = c(-5, 0)),
      size = 15
    ),
    
    glmnet = grid_regular(
      penalty(range = c(-5, 0)),
      mixture(range = c(0, 1)),
      levels = 10
    )
  )
  
  # Combination sizes
  for (n in combo_sizes) {
    size_results <- list()
    size_best_params <- list()
    size_tabs <- list()
    size_fits <- list()
    size_workflows <- list()
    
    # Combinations
    time_combos <- combn(time_predictors, n, simplify = FALSE)
    
    for (combo in time_combos) {
      # Combinations names
      combo_name <- paste(combo, collapse = "_")
      suffix <- ""
      
      # Suffix for derived predictors
      if (include_fixed_predictors && length(fixed_predictors) > 0) {
        suffix <- paste0(suffix, "_fixed(", paste(fixed_predictors, collapse = "+"), ")")
      }
      
      # Derived predictors
      derived <- create_derived_predictors(
        data, combo, include_ratio_predictors, include_diff_predictors
      )
      main_data <- derived$data
      all_predictors <- c(combo, 
                          if(include_fixed_predictors) fixed_predictors,
                          derived$ratio_names, 
                          derived$diff_names)
      
      # Suffix inputation
      if (include_ratio_predictors && length(derived$ratio_names) > 0) {
        suffix <- paste0(suffix, "_ratio")
      }
      if (include_diff_predictors && length(derived$diff_names) > 0) {
        suffix <- paste0(suffix, "_diff")
      }
      
      combo_base_name <- paste0(combo_name, suffix)
      cat("\n--- Processing combination:", combo_base_name, "---\n")
      
      # Select data
      df_subset <- main_data %>% 
        select(all_of(c(all_predictors, response))) %>%
        na.omit()
      
      # Train/test split
      set.seed(123)
      split <- initial_split(df_subset, strata = !!sym(response), prop = 0.75)
      train <- training(split)
      test <- testing(split)
      
      # Recipe
      rec <- recipe(as.formula(paste(response, "~ .")), data = df_subset) %>%
        step_normalize(all_numeric_predictors()) %>%
        step_zv(all_numeric_predictors())
      
      # Algorithm processing
      for (alg in algorithms) {
        alg_name <- paste0(combo_base_name,alg)
        cat("  Algorithm:", alg, "\n")
        
        tryCatch({
          # Workflow
          wf <- workflow() %>% 
            add_recipe(rec) %>% 
            add_model(model_specs[[alg]])
          
          # Crossvalidation
          set.seed(456)
          folds <- vfold_cv(train, strata = !!sym(response))
          
          # Tuning          
          tune_res <- tune_grid(
            # microbiome_wf,
            wf,
            resamples = folds,
            grid = grids[[alg]],
            metrics = metric_set(roc_auc, accuracy, sens, spec),
            control = control_grid(verbose = TRUE, allow_par = T)
          )
          
          
          # Best parameters based on ROC AUC
          best_params <- select_best(tune_res, metric = "roc_auc")
          final_wf <- finalize_workflow(wf, best_params)
          
          # Fit
          final_fit <- fit(final_wf, train)
          
          # Save WF 
          if (return_fits) {
            size_fits[[alg_name]] <- final_fit
          }
          if (return_workflows) {
            size_workflows[[alg_name]] <- final_wf
          }
          
          # CV metrics
          cv_res <- fit_resamples(final_wf, folds, control = control_resamples(verbose=T, save_pred = T))
          cv_preds <- collect_predictions(cv_res)
          cv_metrics <- calculate_class_metrics(cv_preds) %>% 
            mutate(val = "CV10F", Combination = alg_name, Algoritmo = alg)
          
          # Test metrics
          test_preds <- predict(final_fit, test, type = "prob") %>%
            bind_cols(predict(final_fit, test, type = "class"), test)
          test_metrics <- calculate_class_metrics(test_preds) %>% 
            mutate(val = "Test", Combination = alg_name, Algoritmo = alg)
          
          # Possible external validation
          external_metrics <- list()
          for (df_name in names(external_dfs)) {
            df <- external_dfs[[df_name]]
            
            # Derived predictors for external data
            derived_ext <- create_derived_predictors(
              df, combo, include_ratio_predictors, include_diff_predictors
            )
            df_ext <- derived_ext$data
            
            #External validation metrics
            missing_vars <- setdiff(all_predictors, colnames(df_ext))
            
            preds <- predict(final_fit, df_ext) %>% 
              bind_cols(df_ext %>% select(all_of(response)))
            metrics <- calculate_class_metrics(preds) %>% 
              mutate(val = df_name, Combination = alg_name, Algoritmo = alg)
            
            external_metrics[[df_name]] <- metrics
          }
          
          # Store results
          size_results[[alg_name]] <- list(
            cv_metrics = cv_metrics,
            test_metrics = test_metrics,
            external_metrics = external_metrics
          )
          size_best_params[[alg_name]] <- best_params
          all_metrics <- c(list(cv_metrics, test_metrics), external_metrics)
          size_tabs[[alg_name]] <- bind_rows(all_metrics)
          
        }, error = function(e) {
          message("Erro no algoritmo ", alg, " para combinação ", combo_base_name, ": ", e$message)
        })
      }
    }
    
    # Store results for each size
    all_results[[paste0("size_", n)]] <- size_results
    all_best_params[[paste0("size_", n)]] <- size_best_params
    all_tabs[[paste0("size_", n)]] <- bind_rows(size_tabs)
    if (return_fits) {
      all_fits[[paste0("size_", n)]] <- size_fits
    }
    if (return_workflows) {
      all_workflows[[paste0("size_", n)]] <- size_workflows
    }
  }
  if (return_fits && return_workflows) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      fits = all_fits,
      workflows = all_workflows
    ))
  } else if (return_fits) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      fits = all_fits
    ))
  } else if (return_workflows) {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs,
      workflows = all_workflows
    ))
  } else {
    return(list(
      results = all_results,
      best_params_list = all_best_params,
      all_tabs = all_tabs
    ))
  }
}

results <- test_combinations(
  data = asvres,
  time_predictors = c("ASV669", "ASV2294","ASV273","ASV765","ASV652"),
  fixed_predictors = NULL,
  response = "response",
  external_dfs = list(),
  combo_sizes = c(5),
  return_fits = T,
  return_workflows = T,
  include_fixed_predictors = F,
  include_ratio_predictors = F,
  include_diff_predictors = F,
  algorithms = c("glmnet"))

final_results <- bind_rows(results$all_tabs)
fwrite(final_results, "results_completembio13c.csv")



final_resultsall <- rbind(final_results,final_resultsrf,final_resultsSVM,final_resultsXGB)
final_results2 <- final_results %>% filter(.metric == "roc_auc")%>% pivot_wider(id_cols = Combination, names_from = val, values_from = .estimate)
final_results2$cvtest <- (final_results2$CV10F + final_results2$Test)/2



top_combinations_both <- final_results2 %>%
  arrange(desc(cvtest)) %>%
  slice_head(n = 5) %>%
  select(Combination)
top_combinations_cv <- final_results %>%
  filter(.metric == "roc_auc", val == "CV10F") %>%
  arrange(desc(.estimate)) %>%
  slice_head(n = 5) %>%
  select(Combination)
top_combinations_test <- final_results %>%
  filter(.metric == "roc_auc", val == "Test") %>%
  arrange(desc(.estimate)) %>%
  slice_head(n = 5) %>%
  select(Combination)
top_combinations <- rbind(top_combinations_both,top_combinations_cv,top_combinations_test)

top_pairs <- final_results %>%
  filter(Combination %in% top_combinations$Combination) %>%
  filter(val %in% c("CV10F", "Test")) %>%
  arrange(Combination, val)

top_roc <- top_pairs %>% filter(.metric == "roc_auc")