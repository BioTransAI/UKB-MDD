# =============================================================================
# UK Biobank proteomics and major depressive disorder
# Analysis workflow accompanying the manuscript
# =============================================================================
# Purpose
#   Describe the two-stage pathway-based prediction workflow and its downstream
#   phenotype, genetic, mediation, and druggability analyses.
#
# Reading this script
#   Sections follow the workflow described in read.md. The implemented analysis
#   logic, model settings, statistical formulas, and thresholds are retained.
#   Descriptive <PLACEHOLDERS> replace study-specific locations and file names.
#   Replace them with the corresponding local inputs/outputs before execution.
#   The same placeholder denotes the same resource wherever it appears.
#
#   Section 10 is a comment-only outline because its implementation was not
#   included. Section 11 retains the mediation code as a commented template,
#   with protein pairs to be specified. GREP commands are terminal examples.
#   See read.md for the workflow overview and input descriptions.
#
# Prepared inputs
#   Supply train_data_final and test_data_final before running this workflow.
#   Each is a data frame with label (0 = control; 1 = MDD) in the first column,
#   protein features in the remaining columns, and participant IDs in row names.
#   Supply confounder_metadata with prepared covariates and participant IDs in
#   row names for analyses that use covariate adjustment.
# =============================================================================
library(data.table)
library(parallel)
library(caret)
library(mlr3verse)
library(BioM2)

# ----------------------------------------------------------
# 1. Define model-fitting functions
# ----------------------------------------------------------

# Fit a base learner or obtain internal cross-validation predictions.
fit_base_model <- function(train_data, test_data = NULL, prediction_mode = c("classification", "regression", "probability"), classifier, param_list = NULL, inner_folds = 10) {
  prediction_mode <- match.arg(prediction_mode)

  if (!is.null(test_data)) {
    if (colnames(train_data)[1] != "label") stop("The first column of 'train_data' must be 'label'!")
    if (colnames(test_data)[1] != "label") stop("The first column of 'test_data' must be 'label'!")

    if (prediction_mode == 'probability') {
      if (is.character(classifier)) {
        classifier_name <- paste0('classif.', classifier)
        model <- lrn(classifier_name, predict_type = "prob")
      } else {
        model <- classifier
      }

      train_data[, 1] <- as.factor(train_data[, 1])
      test_data[, 1] <- as.factor(test_data[, 1])
      task_train <- as_task_classif(train_data, target = 'label')
      task_test <- as_task_classif(test_data, target = 'label')

      if (!is.null(param_list)) {
        auto_tuner_instance <- auto_tuner(
          tuner = tnr("grid_search", resolution = 10, batch_size = 5),
          learner = model,
          search_space = param_list,
          resampling = rsmp("cv", folds = 5),
          measure = msr("classif.auc")
        )
        auto_tuner_instance$train(task_train)
        model$param_set$values <- auto_tuner_instance$tuning_result$learner_param_vals[[1]]
        model$train(task_train)
      } else {
        sink(nullfile())
        model$train(task_train)
        sink()
      }
      return(model$predict(task_test)$prob[, 2])

    } else if (prediction_mode == 'regression') {
      classifier_name <- paste0('regr.', classifier)
      train_data[, 1] <- as.numeric(train_data[, 1])
      test_data[, 1] <- as.numeric(test_data[, 1])
      task_train <- as_task_regr(train_data, target = 'label')
      task_test <- as_task_regr(test_data, target = 'label')
      model <- lrn(classifier_name)

      if (!is.null(param_list)) {
        auto_tuner_instance <- auto_tuner(
          tuner = tnr("grid_search", resolution = 5, batch_size = 5),
          learner = model,
          search_space = param_list,
          resampling = rsmp("cv", folds = 5),
          measure = msr("regr.mae")
        )
        auto_tuner_instance$train(task_train)
        model$param_set$values <- auto_tuner_instance$tuning_result$learner_param_vals[[1]]
        model$train(task_train)
      } else {
        model$train(task_train)
      }
      return(model$predict(task_test)$response)
    }

  } else {
    if (colnames(train_data)[1] != "label") stop("The first column of 'train_data' must be 'label'!")

    if (prediction_mode == 'probability') {
      classifier_name <- paste0('classif.', classifier)
      train_data[, 1] <- as.factor(train_data[, 1])
      task_train <- as_task_classif(train_data, target = 'label')
      model <- lrn(classifier_name, predict_type = "prob")

      sink(nullfile())
      resampling_result <- resample(task_train, model, rsmp("cv", folds = inner_folds))$prediction()
      sink()

      results_df <- as.data.frame(as.data.table(resampling_result))[, c(1, 5)]
      results_df <- results_df[order(results_df$row_ids), ][, 2]
      return(results_df)
    }
  }
}

# Fit cross-validated learners and average their test-set predictions.
fit_model_with_cv <- function(train_data, test_data, classifier, inner_folds = 10, param_list = NULL) {
  if (is.character(classifier)) {
    classifier_name <- paste0('classif.', classifier)
    model <- lrn(classifier_name, predict_type = "prob")
  } else {
    model <- classifier
  }

  num_train_rows <- nrow(train_data)
  num_test_rows <- nrow(test_data)

  tryCatch({
    sink(nullfile())
    train_data[, 1] <- as.factor(train_data[, 1])
    test_data[, 1] <- as.factor(test_data[, 1])

    task_train <- as_task_classif(train_data, target = 'label')
    task_test <- as_task_classif(test_data, target = 'label')

    if (!is.null(param_list)) {
      auto_tuner_instance <- auto_tuner(
        tuner = tnr("grid_search", resolution = 10, batch_size = 5),
        learner = model,
        search_space = param_list,
        resampling = rsmp("cv", folds = 5),
        measure = msr("classif.auc")
      )
      auto_tuner_instance$train(task_train)
      model$param_set$values <- auto_tuner_instance$tuning_result$learner_param_vals[[1]]
    }

    resampling_result <- resample(task_train, model, rsmp("cv", folds = inner_folds), store_models = TRUE, store_backends = FALSE)
    predictions_df <- as.data.frame(as.data.table(resampling_result$prediction()))[, c(1, 5)]

    result_list <- list(
      pred_train = predictions_df[order(predictions_df$row_ids), ][, 2],
      pred_test = rowMeans(sapply(1:inner_folds, function(x) resampling_result$learners[[x]]$predict(task_test)$prob[, 2]))
    )
    sink()
    return(result_list)

  }, error = function(e) {
    sink()
    return(list(
      pred_train = rep(0, num_train_rows),
      pred_test = rep(0, num_test_rows)
    ))
  })
}

# Fit one classification learner and return training/test predictions.
fit_simple_model <- function(train_data, test_data, classifier, param_list = NULL) {
  if (is.character(classifier)) {
    classifier_name <- paste0('classif.', classifier)
    model <- lrn(classifier_name, predict_type = "prob")
  } else {
    model <- classifier
  }

  num_train_rows <- nrow(train_data)
  num_test_rows <- nrow(test_data)

  tryCatch({
    train_data[, 1] <- as.factor(train_data[, 1])
    test_data[, 1] <- as.factor(test_data[, 1])

    task_train <- as_task_classif(train_data, target = 'label')
    task_test <- as_task_classif(test_data, target = 'label')

    sink(nullfile())
    if (!is.null(param_list)) {
      auto_tuner_instance <- auto_tuner(
        tuner = tnr("grid_search", resolution = 10, batch_size = 5),
        learner = model,
        search_space = param_list,
        resampling = rsmp("cv", folds = 5),
        measure = msr("classif.auc")
      )
      auto_tuner_instance$train(task_train)
      model$param_set$values <- auto_tuner_instance$tuning_result$learner_param_vals[[1]]
    }

    model$train(task_train)
    sink()

    return(list(
      pred_train = model$predict(task_train)$prob[, 2],
      pred_test = model$predict(task_test)$prob[, 2]
    ))

  }, error = function(e) {
    sink()
    return(list(
      pred_train = rep(0, num_train_rows),
      pred_test = rep(0, num_test_rows)
    ))
  })
}

# ----------------------------------------------------------
# 2. Define feature-selection functions
# ----------------------------------------------------------

# Select reconstructed pathway features for the second-stage model.
step2_feature_selection <- function(selection_method = NULL, feature_data = NULL, target_labels = NULL, target_cutoff = NULL, prediction_mode = NULL, classifier_name = NULL) {
  if (prediction_mode %in% c('probability', 'classification')) {

    if (selection_method == 'cor') {
      print(paste0('     Using << correlation >>, selected cutoff => ', target_cutoff))
      correlation_values <- sapply(seq_along(feature_data), function(x) cor(feature_data[[x]], target_labels, method = 'pearson'))
      top_indices <- order(correlation_values, decreasing = TRUE)[seq_len(target_cutoff)]
      positive_indices <- which(correlation_values > 0)
      final_indices <- intersect(top_indices, positive_indices)
      print(paste0('     |> Final number of pathways >>> ', length(final_indices),
                   ' | Min correlation >>> ', round(min(correlation_values[final_indices]), digits = 3)))
      return(final_indices)

    } else if (selection_method == 'wilcox.test') {
      combined_data <- as.data.frame(cbind(label = target_labels, do.call(cbind, feature_data)))
      class_0_data <- combined_data[combined_data$label == unique(combined_data$label)[1], ]
      class_1_data <- combined_data[combined_data$label == unique(combined_data$label)[2], ]

      p_values <- unlist(mclapply(2:ncol(combined_data), function(x) wilcox.test(class_0_data[, x], class_1_data[, x])$p.value, mc.cores = 20))
      correlation_values <- stats::cor(combined_data$label, combined_data[, -1])

      ordered_p_indices <- order(p_values)
      negative_corr_indices <- which(correlation_values < 0)
      valid_indices <- setdiff(ordered_p_indices, negative_corr_indices)

      if (target_cutoff < length(which(correlation_values > 0))) {
        final_indices <- valid_indices[seq_len(target_cutoff)]
      } else {
        final_indices <- valid_indices
      }

      print(paste0('     |> Final number of pathways >>> ', length(final_indices),
                   ' | Max p-value >>> ', round(max(p_values[final_indices]), digits = 3)))
      return(final_indices)

    } else if (selection_method == 'new') {
      if (!is.null(target_labels)) {
        correlation_values <- sapply(seq_along(feature_data), function(x) cor(feature_data[[x]], target_labels, method = 'pearson'))
        combined_features <- do.call(cbind, feature_data)
        corr_matrix <- cor(combined_features)
        redundant_indices <- findCorrelation(corr_matrix, cutoff = target_cutoff)

        positive_indices <- which(correlation_values > 0)
        final_indices <- setdiff(positive_indices, redundant_indices)

        print(paste0('     |> Final number of pathways >>> ', length(final_indices),
                     ' | Min correlation >>> ', round(min(correlation_values[final_indices]), digits = 3)))
        return(final_indices)
      } else {
        combined_features <- do.call(rbind, feature_data)
        labels <- combined_features[, 1]
        correlation_values <- sapply(2:ncol(combined_features), function(x) cor(combined_features[, x], labels, method = 'pearson'))

        positive_indices <- which(correlation_values > 0)
        features_only <- combined_features[, -1]
        corr_matrix <- cor(features_only)
        redundant_indices <- findCorrelation(corr_matrix, cutoff = target_cutoff)

        final_indices <- setdiff(positive_indices, redundant_indices) + 1
        return(final_indices)
      }

    } else if (selection_method == 'RemoveLinear') {
      if (!is.null(target_labels)) {
        correlation_values <- sapply(seq_along(feature_data), function(x) cor(feature_data[[x]], target_labels, method = 'pearson'))
        combined_features <- do.call(cbind, feature_data)
        linear_combos <- findLinearCombos(combined_features)$remove

        positive_indices <- which(correlation_values > 0)
        final_indices <- setdiff(positive_indices, linear_combos)

        print(paste0('     |> Final number of pathways >>> ', length(final_indices),
                     ' | Min correlation >>> ', round(min(correlation_values[final_indices]), digits = 3)))
        return(final_indices)
      } else {
        combined_features <- do.call(rbind, feature_data)
        labels <- combined_features[, 1]
        correlation_values <- sapply(2:ncol(combined_features), function(x) cor(combined_features[, x], labels, method = 'pearson'))

        positive_indices <- which(correlation_values > 0)
        features_only <- combined_features[, -1]
        linear_combos <- findLinearCombos(features_only)$remove

        final_indices <- setdiff(positive_indices, linear_combos) + 1
        return(final_indices)
      }

    } else {
      if (!is.null(target_labels)) {
        correlation_values <- sapply(seq_along(feature_data), function(x) cor(feature_data[[x]], target_labels, method = 'pearson'))
        upper_bound <- 100
        valid_indices <- which(correlation_values > 0 & correlation_values < upper_bound)

        print(paste0('     |> Final number of pathways >>> ', length(valid_indices),
                     ' | Min correlation >>> ', round(min(correlation_values[valid_indices]), digits = 3)))
        return(valid_indices)
      } else {
        combined_features <- do.call(rbind, feature_data)
        labels <- combined_features[, 1]
        correlation_values <- sapply(2:ncol(combined_features), function(x) cor(combined_features[, x], labels, method = 'pearson'))
        return(which(correlation_values > 0) + 1)
      }
    }
  } else if (prediction_mode == 'regression') {
    # Regression feature selection is not specified in this workflow.
  }
}

# Select optional protein features that are not mapped to pathways.
add_unmapped_features <- function(train_df = NULL, test_df = NULL, unmapped_limit = NULL, selection_method = NULL, annotations = NULL, total_length = NULL, core_count = 30) {

  unmapped_train <- train_df[, setdiff(colnames(train_df), annotations$ID)]
  unmapped_test <- test_df[, setdiff(colnames(train_df), annotations$ID)]

  if (is.null(selection_method) || selection_method == 'wilcox.test') {
    class_0_data <- unmapped_train[unmapped_train$label == unique(unmapped_train$label)[1], ]
    class_1_data <- unmapped_train[unmapped_train$label == unique(unmapped_train$label)[2], ]

    p_values <- unlist(mclapply(2:ncol(unmapped_train), function(x) wilcox.test(class_0_data[, x], class_1_data[, x])$p.value, mc.cores = core_count))

    if (is.null(unmapped_limit)) unmapped_limit <- total_length - 1
    unmapped_limit <- min(unmapped_limit, length(p_values))

    selected_indices <- order(p_values)[seq_len(unmapped_limit)] + 1
    return(list('train' = unmapped_train[, selected_indices], 'test' = unmapped_test[, selected_indices]))

  } else if (selection_method == 'cor') {
    correlations <- unlist(mclapply(2:ncol(unmapped_train), function(x) cor(unmapped_train$label, unmapped_train[, x]), mc.cores = core_count))
    abs_correlations <- abs(correlations)

    if (is.null(unmapped_limit)) unmapped_limit <- total_length - 1
    unmapped_limit <- min(unmapped_limit, length(abs_correlations))

    selected_indices <- order(abs_correlations, decreasing = TRUE)[seq_len(unmapped_limit)] + 1
    return(list('train' = unmapped_train[, selected_indices], 'test' = unmapped_test[, selected_indices]))
  }
}

# ----------------------------------------------------------
# 3. Compare two-stage model configurations by cross-validation
# ----------------------------------------------------------
print(table(train_data_final$label))

# Load database and annotations
pathway_db <- readRDS('<GO_BIOLOGICAL_PROCESS_GENE_SETS>')
feature_annotations <- readRDS('<PROTEIN_GENE_ANNOTATIONS>')
feature_annotations$ID <- gsub('[\\.\\_\\-\\/\\-]', '', feature_annotations$ID)
colnames(train_data_final) <- gsub('[\\.\\_\\-\\/\\-]', '', colnames(train_data_final))

# Analysis settings
global_prediction_mode <- "probability"
num_cv_folds <- 5
adjust_confounders <- FALSE
pathway_size_max <- 2000
pathway_size_min <- 10
min_features_per_pathway <- 5
num_processing_cores <- 50
stage1_classifiers <- c('glmnet')
stage1_selection_method <- 'wilcox'
stage1_cutoffs <- c(0.05, 0.01, 0.001)

unmapped_method <- 'wilcox'
add_unmapped_flag <- 'None'
unmapped_counts <- c(0)

stage2_selection_method <- 'wilcox.test'
stage2_cutoffs <- c(0, 50, 100, 500)
use_inner_cv <- 'Yes'
inner_cv_folds <- 100
stage2_classifiers <- c('glmnet')
tuning_param_list <- NULL
tuning_param_list2 <- NULL

final_results_list <- list()

for (classifier_idx in seq_along(stage1_classifiers)) {
  current_classifier <- stage1_classifiers[[classifier_idx]]
  cutoff_evaluations <- list()

  for (cutoff_idx in seq_along(stage1_cutoffs)) {
    print(paste("Evaluating Classifier:", current_classifier))
    set.seed(666)
    current_cutoff <- stage1_cutoffs[[cutoff_idx]]
    fold_assignments <- createFolds(train_data_final$label, k = num_cv_folds)
    fold_results <- list()
    start_time <- Sys.time()

    is_positive_majority <- sum(train_data_final$label == 1) > sum(train_data_final$label == 0)

    for (fold_idx in seq_len(num_cv_folds)) {
      print('Step 1: Reading and Splitting Data')
      train_fold_data <- train_data_final[unlist(fold_assignments[-fold_idx]), ]
      test_fold_data <- train_data_final[unlist(fold_assignments[fold_idx]), ]

      genes_per_pathway <- sapply(seq_along(pathway_db), function(i) length(pathway_db[[i]]))
      valid_pathways_db <- pathway_db[genes_per_pathway > pathway_size_min & genes_per_pathway < pathway_size_max]
      print(paste0('     |> Total valid pathways => ', length(valid_pathways_db)))

      print('Step 2: Feature Selection (Stage 1)')
      if (stage1_selection_method == 'cor') {
        print(paste0('      Using << correlation >>, cutoff: ', current_cutoff))
        feature_correlations <- abs(cor(train_fold_data$label, train_fold_data))
        names(feature_correlations) <- colnames(train_fold_data)
        selected_features <- feature_correlations[feature_correlations > current_cutoff]
        selected_feature_names <- names(selected_features)
      } else {
        print(paste0('      Using << wilcox.test >>, cutoff: ', current_cutoff))
        train_class_0 <- train_fold_data[train_fold_data$label == 0, ]
        train_class_1 <- train_fold_data[train_fold_data$label == 1, ]
        feature_correlations <- unlist(mclapply(1:ncol(train_fold_data), function(x) wilcox.test(train_class_0[, x], train_class_1[, x])$p.value, mc.cores = 30))
        names(feature_correlations) <- colnames(train_fold_data)
        selected_features <- feature_correlations[feature_correlations < current_cutoff]
        selected_feature_names <- names(selected_features)
      }

      fold_annotations <- feature_annotations[feature_annotations$ID %in% colnames(train_fold_data), ]
      min_features_required <- min_features_per_pathway + 1

      pathway_feature_indices <- mclapply(seq_along(valid_pathways_db), function(x) {
        pathway_ids <- c('label', fold_annotations$ID[fold_annotations$entrezID %in% valid_pathways_db[[x]]])
        if (length(pathway_ids) > min_features_required) {
          filtered_ids <- pathway_ids[pathway_ids %in% selected_feature_names]
          if (length(filtered_ids) < min_features_required) {
            correlation_subset <- feature_correlations[pathway_ids]
            if (stage1_selection_method == 'cor') {
              return(names(correlation_subset)[order(correlation_subset, decreasing = TRUE)[1:min_features_required]])
            } else {
              return(names(correlation_subset)[order(correlation_subset, decreasing = FALSE)[1:min_features_required]])
            }
          } else {
            return(filtered_ids)
          }
        } else {
          return(pathway_ids)
        }
      }, mc.cores = 10)

      pathway_lengths <- sapply(seq_along(pathway_feature_indices), function(x) length(pathway_feature_indices[[x]]))

      print('Step 3: Merging Pathway Data')
      train_pathway_list <- mclapply(seq_along(pathway_feature_indices), function(x) train_fold_data[, pathway_feature_indices[[x]]], mc.cores = 10)
      test_pathway_list <- mclapply(seq_along(pathway_feature_indices), function(x) test_fold_data[, pathway_feature_indices[[x]]], mc.cores = 10)

      names(train_pathway_list) <- names(valid_pathways_db)
      names(test_pathway_list) <- names(valid_pathways_db)

      train_pathway_list <- train_pathway_list[pathway_lengths > min_features_per_pathway]
      test_pathway_list <- test_pathway_list[pathway_lengths > min_features_per_pathway]
      filtered_pathway_lengths <- sapply(seq_along(train_pathway_list), function(i) length(train_pathway_list[[i]]))

      print(paste0('     |> Total selected pathways => ', length(train_pathway_list)))
      print(paste0('     |> Min features per pathway => ', min(filtered_pathway_lengths) - 1, ' | Max => ', max(filtered_pathway_lengths) - 1))

      print('Step 4: Model Reconstruction')
      if (use_inner_cv == 'Yes') {
        print('     |> Utilizing Inner Cross-Validation')
        reconstructed_train <- mclapply(seq_along(train_pathway_list), function(i) fit_base_model(train_data = train_pathway_list[[i]], test_data = NULL, prediction_mode = global_prediction_mode, classifier = current_classifier, inner_folds = inner_cv_folds, param_list = tuning_param_list), mc.cores = num_processing_cores)
        reconstructed_test <- mclapply(seq_along(test_pathway_list), function(i) fit_base_model(train_data = train_pathway_list[[i]], test_data = test_pathway_list[[i]], prediction_mode = global_prediction_mode, classifier = current_classifier, param_list = tuning_param_list), mc.cores = num_processing_cores)
      } else {
        reconstructed_preds <- mclapply(seq_along(train_pathway_list), function(i) fit_simple_model(train_data = train_pathway_list[[i]], test_data = test_pathway_list[[i]], classifier = current_classifier, param_list = tuning_param_list), mc.cores = num_processing_cores)
        reconstructed_train <- lapply(seq_along(reconstructed_preds), function(x) reconstructed_preds[[x]]$pred_train)
        reconstructed_test <- lapply(seq_along(reconstructed_preds), function(x) reconstructed_preds[[x]]$pred_test)
      }
      print('     <<< Reconstruction Completed >>>     ')

      print('Step 5: Pathway Feature Selection (Stage 2)')
      # Setup Unmapped Features (if needed)
      if (add_unmapped_flag == 'Yes') {
        db_mapped_genes <- unique(unlist(valid_pathways_db))
        anno_mapped_genes <- unique(fold_annotations$entrezID)
        intersected_genes <- intersect(anno_mapped_genes, db_mapped_genes)
        mapped_ids <- fold_annotations$ID[fold_annotations$entrezID %in% intersected_genes]

        unmapped_train_data <- train_fold_data[, setdiff(colnames(train_fold_data), mapped_ids)]
        unmapped_test_data <- test_fold_data[, setdiff(colnames(train_fold_data), mapped_ids)]
        print(paste0('Unmapped features count: ', ncol(unmapped_train_data)))

        if (unmapped_method == 'cor') {
          unmapped_pvals <- 1 / abs(cor(unmapped_train_data$label, unmapped_train_data[, -1]))
        } else {
          unmapped_class_0 <- unmapped_train_data[unmapped_train_data$label == unique(unmapped_train_data$label)[1], ]
          unmapped_class_1 <- unmapped_train_data[unmapped_train_data$label == unique(unmapped_train_data$label)[2], ]
          unmapped_pvals <- unlist(mclapply(2:ncol(unmapped_train_data), function(x) wilcox.test(unmapped_class_0[, x], unmapped_class_1[, x])$p.value, mc.cores = num_processing_cores))
        }
      }

      if (adjust_confounders) {
        print('Applying Confounder Adjustment')
        adjusted_metadata <- confounder_metadata[rownames(test_fold_data), ]
        print(paste0('Confounders utilized: ', paste(colnames(adjusted_metadata), collapse = ", ")))
      }

      stage2_results_list <- list()
      for (s2_idx in seq_along(stage2_cutoffs)) {
        if (stage2_cutoffs[s2_idx] == 0) {
          selected_stage2_indices <- step2_feature_selection(selection_method = 'None', feature_data = reconstructed_train, target_labels = train_pathway_list[[1]]$label, target_cutoff = stage2_cutoffs[s2_idx], prediction_mode = 'probability', classifier_name = current_classifier)
        } else {
          selected_stage2_indices <- step2_feature_selection(selection_method = stage2_selection_method, feature_data = reconstructed_train, target_labels = train_pathway_list[[1]]$label, target_cutoff = stage2_cutoffs[s2_idx], prediction_mode = 'probability', classifier_name = current_classifier)
        }

        active_stage2_classifiers <- if (is.null(stage2_classifiers)) current_classifier else stage2_classifiers
        total_evaluations <- length(unmapped_counts) * length(active_stage2_classifiers)
        evaluation_record <- data.frame(unmapped_num = 1:total_evaluations, stage2_learner = 1:total_evaluations, AUC = 1:total_evaluations, ACC = 1:total_evaluations, PCC = 1:total_evaluations, BAC = 1:total_evaluations, PRAUC = 1:total_evaluations, cutoff = 1:total_evaluations, MCC = 1:total_evaluations)

        for (u_idx in seq_along(unmapped_counts)) {
          merged_new_train <- do.call(cbind, reconstructed_train[selected_stage2_indices])
          colnames(merged_new_train) <- names(train_pathway_list)[selected_stage2_indices]
          merged_new_train <- cbind(label = train_pathway_list[[1]]$label, merged_new_train)

          merged_new_test <- do.call(cbind, reconstructed_test[selected_stage2_indices])
          colnames(merged_new_test) <- names(train_pathway_list)[selected_stage2_indices]
          merged_new_test <- cbind(label = test_pathway_list[[1]]$label, merged_new_test)

          colnames(merged_new_test) <- gsub(':', '', colnames(merged_new_test))
          colnames(merged_new_train) <- gsub(':', '', colnames(merged_new_train))

          if (add_unmapped_flag == 'Yes' && unmapped_counts[u_idx] > 0) {
            actual_unmapped_count <- min(unmapped_counts[u_idx], length(unmapped_pvals))
            unmapped_target_ids <- order(unmapped_pvals)[seq_len(actual_unmapped_count)] + 1
            merged_new_train <- cbind(merged_new_train, unmapped_train_data[, unmapped_target_ids])
            merged_new_test <- cbind(merged_new_test, unmapped_test_data[, unmapped_target_ids])
          }

          for (c2_idx in seq_along(active_stage2_classifiers)) {
            classifier_stage2 <- active_stage2_classifiers[[c2_idx]]
            record_row <- (u_idx - 1) * length(active_stage2_classifiers) + c2_idx

            prediction_results <- fit_base_model(train_data = merged_new_train, test_data = merged_new_test, prediction_mode = 'probability', classifier = classifier_stage2, param_list = tuning_param_list2)

            if (adjust_confounders) {
              adjustment_data <- cbind(y = prediction_results, adjusted_metadata)
              prediction_results <- lm(y ~ ., data = adjustment_data, family = gaussian(link = "identity"))$residual
            }

            actual_labels <- test_pathway_list[[1]]$label
            binary_predictions <- ifelse(prediction_results > 0.5, 1, 0)

            acc_class1 <- sum(binary_predictions[actual_labels == 1] == 1) / sum(actual_labels == 1)
            acc_class0 <- sum(binary_predictions[actual_labels == 0] == 0) / sum(actual_labels == 0)

            evaluation_record[record_row, 6] <- (acc_class1 + acc_class0) / 2
            evaluation_record[record_row, 1] <- unmapped_counts[u_idx]
            evaluation_record[record_row, 2] <- ifelse(is.character(classifier_stage2), classifier_stage2, classifier_stage2$id)
            evaluation_record[record_row, 5] <- cor(actual_labels, prediction_results, method = 'pearson')
            evaluation_record[record_row, 9] <- ModelMetrics::mcc(predicted = prediction_results, actual = actual_labels, cutoff = 0.5)
            evaluation_record[record_row, 3] <- ROCR::performance(ROCR::prediction(prediction_results, actual_labels), 'auc')@y.values[[1]]

            actual_labels_factor <- as.factor(actual_labels)
            evaluation_record[record_row, 7] <- ifelse(adjust_confounders, 0, ifelse(is_positive_majority, mlr3measures::prauc(actual_labels_factor, 1 - prediction_results, '0'), mlr3measures::prauc(actual_labels_factor, prediction_results, '1')))
            binary_predictions_factor <- as.factor(binary_predictions)
            evaluation_record[record_row, 4] <- confusionMatrix(binary_predictions_factor, actual_labels_factor)$overall['Accuracy'][[1]]
            evaluation_record[record_row, 8] <- stage2_cutoffs[s2_idx]
          }
        }
        stage2_results_list[[s2_idx]] <- evaluation_record
      }
      fold_results[[fold_idx]] <- do.call(rbind, stage2_results_list)
    }

    combined_fold_results <- do.call(rbind, fold_results)
    end_time <- Sys.time()
    print(paste("Fold evaluation time:", end_time - start_time))

    aggregated_results <- aggregate(combined_fold_results[, c(3:7, 9)], by = list(unmapped_num = combined_fold_results$unmapped_num, stage2_cutoff = combined_fold_results$cutoff, stage2_learner = combined_fold_results$stage2_learner), mean)
    aggregated_results$stage1_cutoff <- stage1_cutoffs[cutoff_idx]
    print(aggregated_results[, c(1:4, 6:8)])

    cutoff_evaluations[[cutoff_idx]] <- aggregated_results
  }

  compiled_cutoff_evaluations <- do.call(rbind, cutoff_evaluations)
  compiled_cutoff_evaluations$stage1_learner <- ifelse(is.character(current_classifier), current_classifier, current_classifier$id)
  compiled_cutoff_evaluations <- compiled_cutoff_evaluations[, c('stage1_learner', 'stage2_learner', 'stage1_cutoff', 'stage2_cutoff', 'unmapped_num', 'AUC', 'BAC', 'PRAUC', 'PCC')]

  final_results_list[[classifier_idx]] <- compiled_cutoff_evaluations
}

final_compiled_results <- do.call(rbind, final_results_list)
final_compiled_results <- final_compiled_results[order(final_compiled_results$AUC, decreasing = TRUE), ]
print(paste0('Top Performing Learner: ', final_compiled_results[1, 1]))

saveRDS(final_compiled_results, '<CROSS_VALIDATION_RESULTS>')

# ----------------------------------------------------------
# 4. Apply the selected configuration to training and test data
# ----------------------------------------------------------
cv_results_df <- readRDS('<CROSS_VALIDATION_RESULTS>')
optimal_params <- cv_results_df[cv_results_df$stage1_learner == cv_results_df$stage2_learner, ]

pathway_db <- readRDS('<GO_BIOLOGICAL_PROCESS_GENE_SETS>')
feature_annotations <- readRDS('<PROTEIN_GENE_ANNOTATIONS>')
feature_annotations$ID <- gsub('[\\.\\_\\-\\/\\-]', '', feature_annotations$ID)

print(optimal_params[1, ])
colnames(train_data_final) <- gsub('[\\.\\_\\-\\/\\-]', '', colnames(train_data_final))
colnames(test_data_final) <- gsub('[\\.\\_\\-\\/\\-]', '', colnames(test_data_final))

best_classifier_name <- optimal_params[1, 1]
final_classifier <- ifelse(grepl("class", best_classifier_name), lrn(best_classifier_name, predict_type = 'prob', alpha = 0), best_classifier_name)

# 4.1 Evaluate the selected configuration by cross-validation in training data
set.seed(666)
biom2_train_results <- BioM2(
  TrainData = train_data_final,
  TestData = NULL,
  pathlistDB = pathway_db,
  FeatureAnno = feature_annotations,
  classifier = final_classifier,
  nfolds = 5,
  PathwaySizeUp = 2000,
  PathwaySizeDown = 10,
  MinfeatureNum_pathways = 5,
  Add_UnMapped = FALSE,
  Unmapped_num = 50,
  Inner_CV = TRUE,
  inner_folds = 100,
  Stage1_FeartureSelection_Method = 'wilcox.test',
  cutoff = optimal_params[1, 3],
  Stage2_FeartureSelection_Method = "wilcox.test_rank",
  cutoff2 = ifelse(optimal_params[1, 4] == 0, 10000, optimal_params[1, 4]),
  target = 'predict',
  cores = 30
)
saveRDS(biom2_train_results, '<TRAINING_VALIDATION_RESULTS>')

# 4.2 Fit on training data and evaluate the held-out test samples
print('Running Independent Validation...')
set.seed(666)
biom2_test_results <- BioM2(
  TrainData = train_data_final,
  TestData = test_data_final,
  pathlistDB = pathway_db,
  FeatureAnno = feature_annotations,
  classifier = final_classifier,
  nfolds = 5,
  PathwaySizeUp = 2000,
  PathwaySizeDown = 10,
  MinfeatureNum_pathways = 5,
  Add_UnMapped = FALSE,
  Unmapped_num = 0,
  Inner_CV = TRUE,
  inner_folds = 100,
  Stage1_FeartureSelection_Method = 'wilcox.test',
  cutoff = optimal_params[1, 3],
  Stage2_FeartureSelection_Method = "wilcox.test_rank",
  cutoff2 = ifelse(optimal_params[1, 4] == 0, 10000, optimal_params[1, 4]),
  target = 'predict',
  cores = 30
)
saveRDS(biom2_test_results, '<TEST_VALIDATION_RESULTS>')

# ----------------------------------------------------------
# 5. Extract and rank pathway coefficients
# ----------------------------------------------------------
library(BioM2)
library(GO.db)

lasso_model <- lrn('classif.glmnet')
stage2_train_matrix <- biom2_test_results$Stage2_train
colnames(stage2_train_matrix) <- gsub('GO:', 'GO', colnames(stage2_train_matrix))

lasso_task <- as_task_classif(stage2_train_matrix, target = 'label')
lasso_model$train(lasso_task)
fitted_lasso <- lasso_model$model

lasso_coefficients <- coef(fitted_lasso, s = 0.01)
lasso_weights_df <- as.data.frame(as.matrix(lasso_coefficients))
rownames(lasso_weights_df) <- gsub('GO', 'GO:', rownames(lasso_weights_df))
lasso_weights_df$weight <- abs(lasso_weights_df$s1)

# Remove intercept and sort by weight
lasso_weights_df <- lasso_weights_df[-1, ]
lasso_weights_df <- lasso_weights_df[order(lasso_weights_df$weight, decreasing = TRUE), ]

top_pathway_count <- 20
lasso_top_pathways <- rownames(lasso_weights_df)[1:top_pathway_count]
saveRDS(lasso_top_pathways, '<TOP_PATHWAY_IDS>')

# Map GO Terms
lasso_final_weights <- lasso_weights_df[abs(lasso_weights_df$weight) > 0, ]
lasso_final_weights$term <- AnnotationDbi::select(GO.db, keys = rownames(lasso_final_weights), columns = c("GOID", "TERM"), keytype = "GOID")[, 2]
saveRDS(lasso_final_weights, '<PATHWAY_COEFFICIENT_RESULTS>')

# ----------------------------------------------------------
# 6. Extract within-pathway protein coefficients
# ----------------------------------------------------------
# Use pathways with nonzero second-stage coefficients to assemble protein
# datasets, estimate within-pathway coefficients, and collect input protein IDs.
selected_pathways_df <- readRDS('<PATHWAY_COEFFICIENT_RESULTS>')
lasso_top_pathway_names <- rownames(selected_pathways_df)
pathway_db_subset <- pathway_db[lasso_top_pathway_names]

stage1_cutoff_val <- optimal_params[1, 3]
min_features_required <- min_features_per_pathway + 1

# Screen protein features using the configured Wilcoxon threshold.
train_class_0 <- train_data_final[train_data_final$label == 0, ]
train_class_1 <- train_data_final[train_data_final$label == 1, ]
protein_p_values <- unlist(mclapply(1:ncol(train_data_final), function(x) wilcox.test(train_class_0[, x], train_class_1[, x])$p.value, mc.cores = 30))

names(protein_p_values) <- colnames(train_data_final)
selected_proteins <- protein_p_values[protein_p_values < stage1_cutoff_val]
selected_protein_names <- names(selected_proteins)

fold_annotations <- feature_annotations[feature_annotations$ID %in% colnames(train_data_final), ]

protein_pathway_indices <- mclapply(seq_along(pathway_db_subset), function(x) {
  pathway_ids <- c('label', fold_annotations$ID[fold_annotations$entrezID %in% pathway_db_subset[[x]]])
  if (length(pathway_ids) > min_features_required) {
    filtered_ids <- pathway_ids[pathway_ids %in% selected_protein_names]
    if (length(filtered_ids) < min_features_required) {
      p_val_subset <- protein_p_values[pathway_ids]
      return(names(p_val_subset)[order(p_val_subset, decreasing = FALSE)[1:min_features_required]])
    } else {
      return(filtered_ids)
    }
  } else {
    return(pathway_ids)
  }
}, mc.cores = 10)

pathway_protein_train_list <- mclapply(seq_along(protein_pathway_indices), function(x) train_data_final[, protein_pathway_indices[[x]]], mc.cores = 10)
names(pathway_protein_train_list) <- names(pathway_db_subset)

pathway_protein_weights <- mclapply(seq_along(pathway_protein_train_list), function(x) {
  protein_lasso_model <- lrn('classif.glmnet')
  protein_matrix <- pathway_protein_train_list[[x]]

  protein_task <- as_task_classif(protein_matrix, target = 'label')
  protein_lasso_model$train(protein_task)
  fitted_protein_lasso <- protein_lasso_model$model

  protein_coefs <- coef(fitted_protein_lasso, s = 0.01)
  protein_weights_df <- as.data.frame(as.matrix(protein_coefs))
  protein_weights_df$weight <- protein_weights_df$s1
  return(protein_weights_df[-1, ])
}, mc.cores = 10)

names(pathway_protein_weights) <- names(pathway_protein_train_list)
extracted_protein_names <- lapply(seq_along(pathway_protein_train_list), function(x) colnames(pathway_protein_train_list[[x]])[-1])
unique_important_proteins <- unique(unlist(extracted_protein_names))

# ----------------------------------------------------------
# 7. Calculate MDD polygenic risk scores with PRSice-2
# ----------------------------------------------------------
target_dataset <- "<TARGET_GENOTYPE_PREFIX>"
base_gwas_eur <- "<PRS_BASE_GWAS_BETA>"
prs_output_dir <- "<PRS_OUTPUT_PREFIX>"

prsice_cmd_1 <- paste(
  "Rscript <PRSICE_R_SCRIPT>",
  "--dir .",
  "--prsice <PRSICE_EXECUTABLE>",
  "--snp SNP --A1 A1 --bp BP --chr CHR --pvalue P",
  "--thread 20 --stat BETA --binary-target T",
  "--base", base_gwas_eur,
  "--target", target_dataset,
  "--bar-levels 1,0.5,0.1,0.05,0.01,1e-3,1e-4,1e-5,5e-8",
  "--clump-kb 250",
  "--clump-r2 0.1",
  "--no-regress --fastscore",
  "--out", prs_output_dir
)
system(prsice_cmd_1)

base_gwas_2023 <- "<PRS_BASE_GWAS_OR>"

prsice_cmd_2 <- paste(
  "Rscript <PRSICE_R_SCRIPT>",
  "--dir .",
  "--prsice <PRSICE_EXECUTABLE>",
  "--snp SNP --A1 A1 --bp BP --chr CHR --pvalue P",
  "--thread 20 --stat OR --binary-target T",
  "--base", base_gwas_2023,
  "--target", target_dataset,
  "--bar-levels 1,0.5,0.1,0.05,0.01,1e-3,1e-4,1e-5,5e-8",
  "--clump-kb 250",
  "--clump-r2 0.1",
  "--no-regress --fastscore",
  "--out", prs_output_dir
)
system(prsice_cmd_2)

# ----------------------------------------------------------
# 8. Assess associations with phenotypes
# ----------------------------------------------------------
library(MASS)
library(dplyr)

phenotype_data <- readRDS('<PHENOTYPE_DATA>')
biom2_validation_data <- readRDS('<TEST_VALIDATION_RESULTS>')

# Extract risk predictions
if (length(biom2_validation_data) == 3) {
  risk_predictions <- do.call(rbind, biom2_validation_data$Prediction)
} else {
  risk_predictions <- biom2_validation_data$Prediction
}
rownames(risk_predictions) <- risk_predictions$sample
risk_df <- data.frame(id = risk_predictions$sample, risk = risk_predictions$prediction, row.names = risk_predictions$sample)

# Align confounders and phenotypes
aligned_confounders <- confounder_metadata[rownames(risk_df), ]
aligned_phenotypes <- phenotype_data[rownames(risk_df), ]
merged_covariates <- cbind(risk = risk_df$risk, aligned_confounders)

association_results <- mclapply(seq_along(aligned_phenotypes), function(idx) {
  merged_model_data <- cbind(y = aligned_phenotypes[, idx], merged_covariates)
  merged_model_data <- as.data.frame(na.omit(merged_model_data))

  # Remove columns with zero variance
  zero_variance_cols <- which(sapply(seq_len(ncol(merged_model_data)), function(col_idx) length(unique(merged_model_data[, col_idx]))) == 1)
  if (length(zero_variance_cols) > 0) {
    merged_model_data <- merged_model_data[, -zero_variance_cols]
  }

  if (!'y' %in% colnames(merged_model_data)) return(c(0, 1))

  tryCatch({
    if (class(merged_model_data$y) == 'numeric') {
      lm_model <- lm(y ~ ., data = merged_model_data)
      p_val <- summary(lm_model)$coefficients['risk', 4]
      beta_val <- summary(lm_model)$coefficients['risk', 1]
    } else if (class(merged_model_data$y) == 'character') {
      merged_model_data$y <- as.numeric(merged_model_data$y)
      glm_model <- glm(y ~ ., family = binomial(link = "logit"), data = merged_model_data)
      p_val <- summary(glm_model)$coefficients['risk', 4]
      beta_val <- summary(glm_model)$coefficients['risk', 1]
    } else if (class(merged_model_data$y) == 'factor') {
      if (length(unique(merged_model_data$y)) >= 30 || length(unique(merged_model_data$y)) <= 2) {
        merged_model_data$y <- as.numeric(as.character(merged_model_data$y))
        lm_model <- lm(y ~ ., data = merged_model_data)
        p_val <- summary(lm_model)$coefficients['risk', 4]
        beta_val <- summary(lm_model)$coefficients['risk', 1]
      } else {
        merged_model_data$y <- as.factor(as.character(merged_model_data$y))
        polr_model <- polr(y ~ ., data = merged_model_data, Hess = TRUE, model = TRUE, method = "logistic")
        p_val <- pnorm(abs(coef(summary(polr_model))["risk", "t value"]), lower.tail = FALSE) * 2
        beta_val <- summary(polr_model)$coefficients['risk', 1]
      }
    }
    return(c(beta_val, p_val))

  }, error = function(e) {
    merged_model_data$y <- as.numeric(as.character(merged_model_data$y))
    lm_fallback <- lm(y ~ ., data = merged_model_data)
    p_val <- summary(lm_fallback)$coefficients['risk', 4]
    beta_val <- summary(lm_fallback)$coefficients['risk', 1]
    return(c(beta_val, p_val))
  })
}, mc.cores = 10)

phenotype_betas <- unlist(lapply(association_results, function(x) x[1]))
phenotype_pvals <- unlist(lapply(association_results, function(x) x[2]))

phenotype_summary_df <- data.frame(
  phenotype = colnames(phenotype_data),
  beta = phenotype_betas,
  pvalue = phenotype_pvals,
  pvalue_adj = p.adjust(phenotype_pvals, method = "fdr")
)

# ----------------------------------------------------------
# 9. Perform forward Mendelian randomization: protein -> MDD
# ----------------------------------------------------------
library(TwoSampleMR)
library(ggplot2)
library(data.table)
library(openxlsx)
library(ieugwasr)

mdd_gwas_data <- fread('<MDD_GWAS_SUMMARY_STATISTICS>', data.table = FALSE)
mdd_gwas_data$phenotype <- 'depression'

# Load mapping and annotation files
annotation_files <- list.files('<SNP_ANNOTATION_DIRECTORY>', '<SNP_ANNOTATION_FILE_PATTERN>', full.names = TRUE)
snp_annotations <- do.call(rbind, lapply(annotation_files, function(x) fread(x, data.table = FALSE)))

pqtl_files <- list.files('<PQTL_DIRECTORY>', '.tar', full.names = TRUE)
split_filenames <- strsplit(pqtl_files, '_')
protein_names <- sapply(split_filenames, function(x) gsub('<PQTL_DIRECTORY>/', '', x[1]))

lasso_top_proteins <- readRDS('<SELECTED_PROTEIN_SET>')
target_proteins <- intersect(lasso_top_proteins, protein_names)
protein_directories <- list.dirs('<PQTL_DIRECTORY>/')

mr_results_forward <- list()

for (i in seq_along(target_proteins)) {
  start_time <- Sys.time()
  current_protein <- target_proteins[i]
  print(paste("Processing MR (Forward):", current_protein))

  protein_dir <- protein_directories[grepl(current_protein, protein_directories)]
  if (length(protein_dir) > 1) {
    dir_splits <- strsplit(protein_dir, '_')
    dir_prefixes <- unlist(lapply(dir_splits, function(x) x[1]))
    protein_dir <- protein_dir[which.min(nchar(dir_prefixes))]
  }

  if (length(protein_dir) == 0) {
    mr_results_forward[[current_protein]] <- NA
    next
  }

  protein_files <- list.files(protein_dir, '', full.names = TRUE)

  # Filter significant SNPs
  protein_snps_list <- lapply(protein_files, function(x) {
    raw_data <- fread(x, data.table = FALSE)
    significance_threshold <- -log10(5e-8)
    significant_snps <- raw_data[raw_data$LOG10P > significance_threshold, ]
    if (nrow(significant_snps) > 0) return(significant_snps) else return(NA)
  })

  protein_snps_df <- do.call(rbind, protein_snps_list)
  if (is.null(protein_snps_df) || ncol(protein_snps_df) == 1) {
    mr_results_forward[[current_protein]] <- NA
    next
  }

  protein_snps_df <- protein_snps_df[!is.na(protein_snps_df$CHROM), ]
  merged_snps_df <- merge(protein_snps_df, snp_annotations, by = 'ID')
  merged_snps_df$phenotype <- current_protein
  merged_snps_df$P <- 10^(-merged_snps_df$LOG10P)
  merged_snps_df$MAF <- ifelse(merged_snps_df$A1FREQ > 0.5, 1 - merged_snps_df$A1FREQ, merged_snps_df$A1FREQ)
  merged_snps_df$R2 <- 2 * merged_snps_df$MAF * (1 - merged_snps_df$MAF) * (merged_snps_df$BETA^2)
  merged_snps_df$F_sta <- merged_snps_df$R2 * (merged_snps_df$N - 2) / (1 - merged_snps_df$R2)

  # LD Clumping
  tryCatch({
    clumped_snps <- ld_clump(
      dplyr::tibble(rsid = merged_snps_df$rsid, pval = merged_snps_df$P, id = merged_snps_df$phenotype),
      clump_kb = 1000,
      clump_r2 = 0.01,
      plink_bin = genetics.binaRies::get_plink_binary(),
      bfile = "<LD_REFERENCE_PREFIX>"
    )
    merged_snps_df <- merged_snps_df[merged_snps_df$rsid %in% clumped_snps$rsid, ]
  }, error = function(e) {
    print('LD clumping failed, skipping.')
  })

  merged_snps_df <- merged_snps_df[merged_snps_df$F_sta > 10, ]
  if (nrow(merged_snps_df) == 0 || sum(grepl('rs', merged_snps_df$rsid)) == 0) {
    mr_results_forward[[current_protein]] <- NA
    next
  }

  exposure_data <- format_data(
    dat = merged_snps_df, type = "exposure", snps = merged_snps_df$rsid,
    header = TRUE, phenotype_col = "phenotype", snp_col = "rsid",
    beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1",
    other_allele_col = "ALLELE0", pval_col = "P", chr_col = "CHROM", pos_col = "POS38"
  )

  outcome_snps <- mdd_gwas_data[mdd_gwas_data$SNP %in% exposure_data$SNP, ]
  if (nrow(outcome_snps) == 0) {
    mr_results_forward[[current_protein]] <- NA
    next
  }

  outcome_data <- format_data(
    dat = outcome_snps, type = "outcome", snps = outcome_snps$SNP,
    header = TRUE, phenotype_col = "phenotype", snp_col = "SNP",
    beta_col = "logOR", se_col = "SE", effect_allele_col = "EA",
    other_allele_col = "NEA", pval_col = "P", chr_col = "Chromosome", pos_col = "Position"
  )

  harmonized_data <- harmonise_data(exposure_data, outcome_data)
  mr_results_forward[[current_protein]] <- list(
    mr_result = mr(harmonized_data),
    heterogeneity = mr_heterogeneity(harmonized_data),
    pleiotropy_test = mr_pleiotropy_test(harmonized_data)
  )
  print(paste("Completed in:", Sys.time() - start_time))
}

saveRDS(mr_results_forward, '<FORWARD_MR_RESULTS>')

# ----------------------------------------------------------
# 10. Reverse Mendelian randomization: MDD -> protein (outline)
# ----------------------------------------------------------
# The reverse-MR implementation is not included in this script.
# Intended direction: MDD-associated variants as instruments, with protein
# abundance as the outcome. This section is retained as a workflow outline.
# No analysis or result file is generated by this section.

# ----------------------------------------------------------
# 11. Mediation analysis (commented template)
# ----------------------------------------------------------
# Specify the upstream/downstream protein sets before enabling this template.
# The models describe upstream protein -> MDD status -> downstream protein,
# with the supplied covariates included in both models.
#
# library(mediation)
#
# upstream_proteins <- c('MR-Upstream-Protein-Placeholder')
# downstream_proteins <- c('MR-Downstream-Protein-Placeholder')
#
# merged_clinical_data <- rbind(train_data_final, test_data_final)
# merged_clinical_data <- cbind(merged_clinical_data, confounder_metadata[rownames(merged_clinical_data), ])
# covariate_names <- colnames(confounder_metadata)
#
# mediation_results <- list()
#
# for (up_idx in seq_along(upstream_proteins)) {
#   start_time <- Sys.time()
#   current_upstream <- upstream_proteins[up_idx]
#   print(paste("Processing Upstream Protein:", current_upstream))
#
#   downstream_results <- list()
#   for (down_idx in seq_along(downstream_proteins)) {
#     current_downstream <- downstream_proteins[down_idx]
#
#     # Mediator Model (Label ~ Upstream + Covariates)
#     mediator_df <- merged_clinical_data[, c('label', current_upstream, covariate_names)]
#     mediator_fit <- glm(label ~ ., family = binomial(link = "logit"), data = mediator_df)
#
#     # Outcome Model (Downstream ~ Label + Upstream + Covariates)
#     outcome_df <- merged_clinical_data[, c('label', current_upstream, current_downstream, covariate_names)]
#     colnames(outcome_df)[3] <- 'y'
#     outcome_fit <- lm(y ~ ., data = outcome_df)
#
#     # Mediation Analysis
#     med_out <- mediate(mediator_fit, outcome_fit, treat = current_upstream, mediator = "label", robustSE = TRUE, sims = 1000)
#     summary_med <- summary(med_out)
#
#     formatted_result <- data.frame(
#       effect_type = c("ACME (Indirect Effect)", "ADE (Direct Effect)", "Total Effect", "Proportion Mediated"),
#       estimate = c(summary_med$d0, summary_med$z0, summary_med$tau.coef, summary_med$n0),
#       ci_lower = c(summary_med$d0.ci["2.5%"], summary_med$z0.ci["2.5%"], summary_med$tau.ci["2.5%"], summary_med$n0.ci["2.5%"]),
#       ci_upper = c(summary_med$d0.ci["97.5%"], summary_med$z0.ci["97.5%"], summary_med$tau.ci["97.5%"], summary_med$n0.ci["97.5%"]),
#       p_value = c(summary_med$d0.p, summary_med$z0.p, summary_med$tau.p, summary_med$n0.p)
#     )
#
#     downstream_results[[current_downstream]] <- formatted_result
#   }
#   mediation_results[[current_upstream]] <- downstream_results
#   print(paste("Mediation completed in:", Sys.time() - start_time))
# }
#
# saveRDS(mediation_results, '<MEDIATION_RESULTS>')

# ----------------------------------------------------------
# 12. Assess enrichment for druggable targets
# ----------------------------------------------------------
library(tidyr)
library(ggplot2)
library(dplyr)
library(scales)
library(openxlsx)

drug_db <- read.xlsx('<DRUGGABILITY_ANNOTATIONS>', sheet = 1)
selected_proteins <- readRDS('<SELECTED_PROTEIN_SET>')
all_analyzed_proteins <- colnames(train_data_final)
non_selected_proteins <- setdiff(all_analyzed_proteins, selected_proteins)

tier_1_proteins <- drug_db$hgnc_names[drug_db$druggability_tier == 'Tier 1']
tier_2_proteins <- drug_db$hgnc_names[drug_db$druggability_tier == 'Tier 2']
tier_3_proteins <- drug_db$hgnc_names[drug_db$druggability_tier %in% c('Tier 3A', 'Tier 3B')]

calculate_contingency <- function(subgroup, selected, non_selected) {
  in_subgroup_and_selected <- length(intersect(subgroup, selected))
  not_in_subgroup_selected <- length(selected) - in_subgroup_and_selected
  in_subgroup_not_selected <- length(intersect(subgroup, non_selected))
  not_in_subgroup_not_selected <- length(non_selected) - in_subgroup_not_selected

  return(data.frame(
    a = in_subgroup_and_selected,
    b = not_in_subgroup_selected,
    c = in_subgroup_not_selected,
    d = not_in_subgroup_not_selected
  ))
}

tier_statistics <- bind_rows(
  calculate_contingency(tier_1_proteins, selected_proteins, non_selected_proteins) %>% mutate(tier = "Tier1"),
  calculate_contingency(tier_2_proteins, selected_proteins, non_selected_proteins) %>% mutate(tier = "Tier2"),
  calculate_contingency(tier_3_proteins, selected_proteins, non_selected_proteins) %>% mutate(tier = "Tier3")
)

fisher_test_results <- list()
for (i in 1:nrow(tier_statistics)) {
  contingency_matrix <- matrix(
    c(tier_statistics$a[i], tier_statistics$c[i],
      tier_statistics$b[i], tier_statistics$d[i]),
    nrow = 2, byrow = TRUE
  )

  test_res <- fisher.test(contingency_matrix)

  fisher_test_results[[i]] <- data.frame(
    tier = tier_statistics$tier[i],
    overlap_count = tier_statistics$a[i],
    odds_ratio = as.numeric(test_res$estimate),
    ci_lower = test_res$conf.int[1],
    ci_upper = test_res$conf.int[2],
    p_value = test_res$p.value
  )
}
final_fisher_results <- do.call(rbind, fisher_test_results)
print(final_fisher_results)

# Export the selected gene list for the separate GREP analysis.
write.table(selected_proteins, file = "<SELECTED_PROTEIN_GENE_LIST>", quote = FALSE, sep = "\n", col.names = FALSE, row.names = FALSE)

# -------------------------------------------------------------------------
# NOTE: GREP Python Analysis Instructions
# The following commands should be executed in the terminal, not in R:
# -------------------------------------------------------------------------
# cd <GREP_DIRECTORY>
# python grep.py --genelist <SELECTED_PROTEIN_GENE_LIST> \
#                --out <GREP_ATC_OUTPUT_PREFIX> --test ATC --output-drug-name \
#                --background <BACKGROUND_GENE_LIST>
#
# python grep.py --genelist <SELECTED_PROTEIN_GENE_LIST> \
# --out <GREP_ICD_OUTPUT_PREFIX> --test ICD --output-drug-name \
# --background <BACKGROUND_GENE_LIST>
# -------------------------------------------------------------------------
