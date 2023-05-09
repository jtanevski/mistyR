# MISTy runner
# Copyleft (ɔ) 2020 Jovan Tanevski [jovan.tanevski@uni-heidelberg.de]


#' @importFrom rlang !! :=
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("mistyR is able to run computationally intensive functions
  in parallel. Please consider specifying a future::plan(). For example by running
  future::plan(future::multisession) before calling mistyR functions.")
}

#' @importFrom dplyr %>%
#' @inherit run_misty examples
#' @export
dplyr::`%>%`

# allow using tidyselect where
utils::globalVariables("where")

#' Train MISTy models
#'
#' Trains multi-view models for all target markers, estimates the performance,
#' the contributions of the view specific models and the importance of predictor
#' markers for each target marker.
#'
#' If \code{bypass.intra} is set to \code{TRUE} all variable in the intraview
#' the intraview data will be treated as targets only. The baseline intraview
#' model in this case is a trivial model that predicts the average of each
#' target. If the intraview has only one variable this switch is automatically
#' set to \code{TRUE}.
#'
#' \emph{Default} model to train the view-specific views is a Random Forest model
#' based on \code{\link[ranger]{ranger}()} --
#' \code{run_misty(views, model.function = random_forest_model)}
#'
#' The following parameters are the default
#' configuration: \code{num.trees = 100}, \code{importance = "impurity"},
#' \code{num.threads = 1}, \code{seed = seed}.
#'
#' \emph{Gradient boosting} is an alternative to model each view using gradient
#' boosting. The algorithm is based on \code{\link[xgboost]{xgb.train}()} --
#'  \code{run_misty(views, model.function = gradient_boosting_model)}
#'
#' The following parameters are the default configuration: \code{booster = "gbtree"},
#' \code{rounds = 10}, \code{objective = "reg:squarederror"}. Set \code{booster}
#' to \code{"gblinear"} for linear boosting.
#'
#' \emph{Bagged MARS} is an alternative to model each view using bagged MARS,
#' (multivariate adaptive spline regression models) trained with
#' bootstrap aggregation samples. The algorithm is based on
#' \code{\link[earth]{earth}()} --
#' \code{run_misty(views, model.function = bagged_mars_model)}
#'
#' The following parameters are the default configuration: \code{degree = 2}.
#' Furthermore 50 base learners are used by default (pass \code{n.bags} as
#' parameter via \code{...} to change this value).
#'
#' \emph{MARS} is an alternative to model each view using
#' multivariate adaptive spline regression model. The algorithm is based on
#' \code{\link[earth]{earth}()} --
#' \code{run_misty(views, model.function = mars_model)}
#'
#' The following parameters are the default configuration: \code{degree = 2}.
#'
#' \emph{Linear model} is an alternative to model each view using a simple linear
#' model. The algorithm is based on \code{\link[stats]{lm}()} --
#' \code{run_misty(views, model.function = linear_model)}
#'
#' \emph{SVM} is an alternative to model each view using a support vector
#' machines. The algorithm is based on \code{\link[kernlab]{ksvm}()} --
#' \code{run_misty(views, model.function = svm_model)}
#'
#' The following parameters are the default configuration: \code{kernel = "vanilladot"}
#' (linear kernel), \code{C = 1}, \code{type = "eps-svr"}.
#'
#' \emph{MLP} is an alternative to model each view using a multi-layer perceptron.
#' The alogorithm is based on \code{\link[RSNNS]{mlp}()} --
#' \code{run_misty(views, model.function = mlp_model)}
#'
#' The following parameters are the default configuration: \code{size = c(10)}
#' (meaning we have 1 hidden layer with 10 units).
#'
#'
#' @param views view composition.
#' @param results.folder path to the top level folder to store raw results.
#' @param seed seed used for random sampling to ensure reproducibility.
#' @param target.subset subset of targets to train models for. If \code{NULL},
#'     models will be trained for markers in the intraview.
#' @param bypass.intra a \code{logical} indicating whether to train a baseline
#'     model using the intraview data (see Details).
#' @param cv.folds number of cross-validation folds to consider for estimating
#'     the performance of the multi-view models.
#' @param cv.strict a \code{logical} indicating whether to drop markers that have
#'     less than \code{cv.fold} unique values.
#' @param cached a \code{logical} indicating whether to cache the trained models
#'     and to reuse previously cached ones if they already exist for this sample.
#' @param append a \code{logical} indicating whether to append the performance
#'     and coefficient files in the \code{results.folder}. Consider setting to
#'     \code{TRUE} when rerunning a workflow with different \code{target.subset}
#'     parameters.
#' @param model.function a function which is used to model each view, default
#'     model is \code{random_forest_model}. Other models included in mistyR are
#'     \code{gradient_boosting_model}, \code{bagged_mars_model},
#'     \code{mars_model}, \code{linear_model},
#'     \code{svm_model}, \code{mlp_model}.
#' @param ... all additional parameters are passed to the chosen ML model for
#' training the view-specific models.
#'
#' @return Path to the results folder that can be passed to
#'     \code{\link{collect_results}()}.
#'
#' @seealso \code{\link{create_initial_view}()} for
#'     starting a view composition.
#'
#' @examples
#' # Create a view composition of an intraview and a paraview with radius 10 then
#' # run MISTy for a single sample.
#'
#' library(dplyr)
#'
#' # get the expression data
#' data("synthetic")
#' expr <- synthetic[[1]] %>% select(-c(row, col, type))
#' # get the coordinates for each cell
#' pos <- synthetic[[1]] %>% select(row, col)
#'
#' # compose
#' misty.views <- create_initial_view(expr) %>% add_paraview(pos, l = 10)
#'
#' # run with default parameters
#' run_misty(misty.views)
#' @export
run_misty <- function(views, results.folder = "results", seed = 42,
                      target.subset = NULL, bypass.intra = FALSE, cv.folds = 10,
                      cv.strict = FALSE, cached = FALSE, append = FALSE,
                      model.function = random_forest_model,
                      ...) {
  model.name <- as.character(rlang::enexpr(model.function))

  if (!exists(model.name, envir = globalenv())) {
    model.function <- utils::getFromNamespace(model.name, "mistyR")
  }

  on.exit(sweep_cache())

  clean.views <- views %>%
    select_markers(
      "intraview",
      where(~ (stats::sd(.) > 0))
    ) %>%
    select_markers(
      "intraview",
      where(~ (length(unique(.)) >= ifelse(cv.strict, cv.folds, 1)))
    )

  view.abbrev <- clean.views %>%
    rlist::list.remove(c("misty.uniqueid")) %>%
    purrr::map_chr(~ .x[["abbrev"]])


  header <- stringr::str_glue("target intercept {views} p.intercept {p.views}",
    views = paste0(view.abbrev, collapse = " "),
    p.views = paste0("p.", view.abbrev, collapse = " "),
    .sep = " "
  )

  expr <- clean.views[["intraview"]][["data"]]

  check.cv <- assertthat::see_if(nrow(expr) >= cv.folds,
    msg = "The data has less rows than the requested number of cv folds. Stopping."
  )

  if (!check.cv) {
    warning(attr(check.cv, "msg"))
    return(NULL)
  }

  check.var <- assertthat::see_if(ncol(expr) >= 1,
    msg = "The intraview doesn't contain any relevant information that can be modeled. Stopping."
  )

  if (!check.var) {
    warning(attr(check.var, "msg"))
    return(NULL)
  }

  normalized.results.folder <- R.utils::getAbsolutePath(results.folder)

  if (!dir.exists(normalized.results.folder)) {
    dir.create(normalized.results.folder, recursive = TRUE)
  }

  if (ncol(expr) == 1) bypass.intra <- TRUE

  coef.file <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "coefficients.txt"
  )
  coef.lock <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "coefficients.txt.lock"
  )
  on.exit(file.remove(coef.lock), add = TRUE)

  if (!append) {
    current.lock <- filelock::lock(coef.lock)
    write(header, file = coef.file)
    filelock::unlock(current.lock)
  }

  header <- "target intra.RMSE intra.R2 multi.RMSE multi.R2 p.RMSE p.R2"

  perf.file <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "performance.txt"
  )
  perf.lock <- paste0(
    normalized.results.folder, .Platform$file.sep,
    "performance.txt.lock"
  )
  on.exit(file.remove(perf.lock), add = TRUE)

  if (!append) {
    current.lock <- filelock::lock(perf.lock)
    write(header, file = perf.file)
    filelock::unlock(current.lock)
  }

  targets <- switch(class(target.subset),
    "numeric" = colnames(expr)[target.subset],
    "integer" = colnames(expr)[target.subset],
    "character" = target.subset,
    "NULL" = colnames(expr),
    NULL
  )

  message("\nTraining models")
  targets %>% furrr::future_map_chr(function(target, ...) {
    target.model <- build_model(
      views = clean.views, target = target,
      model.function = model.function,
      model.name = model.name,
      cv.folds = cv.folds,
      bypass.intra = bypass.intra,
      seed = seed, cached = cached, ...
    )

    combined.views <- target.model[["meta.model"]]

    model.lm <- methods::is(combined.views, "lm")

    coefs <- stats::coef(combined.views) %>% tidyr::replace_na(0)

    pvals <- if (model.lm) {
      # fix for missing pvals
      combined.views.summary <- summary(combined.views)
      pvals <- data.frame(c = stats::coef(combined.views)) %>%
        tibble::rownames_to_column("views") %>%
        dplyr::left_join(
          data.frame(p = stats::coef(combined.views.summary)[, 4]) %>%
            tibble::rownames_to_column("views"),
          by = "views"
        ) %>%
        dplyr::pull(p) %>%
        tidyr::replace_na(1)

      if (bypass.intra) append(pvals[-1], c(NA, 1), 0) else c(NA, pvals)
    } else {
      pvals <- ridge::pvals(combined.views)$pval[, combined.views$chosen.nPCs]
      if (bypass.intra) append(pvals, c(NA, 1), 0) else c(NA, pvals)
    }


    # coefficient values and p-values
    coeff <- c(
      if (bypass.intra) append(coefs, 0, 1) else coefs,
      pvals
    )

    current.lock <- filelock::lock(coef.lock)
    write(paste(target, paste(coeff, collapse = " ")),
      file = coef.file, append = TRUE
    )
    filelock::unlock(current.lock)

    # raw importances
    target.model[["model.importances"]] %>% purrr::walk2(
      view.abbrev,
      function(model.importance, abbrev) {
        targets <- names(model.importance)

        imps <- tibble::tibble(
          target = targets,
          imp = model.importance
        )

        readr::write_csv(
          imps,
          paste0(
            normalized.results.folder, .Platform$file.sep,
            "importances_", target, "_", abbrev, ".txt"
          )
        )
      }
    )


    # performance
    if (sum(target.model[["performance.estimate"]] < 0 |
      is.na(target.model[["performance.estimate"]])) > 0) {
      warning.message <-
        paste(
          "Negative performance detected and replaced with 0 for target",
          target
        )
      warning(warning.message)
    }

    performance.estimate <- target.model[["performance.estimate"]] %>%
      dplyr::mutate(dplyr::across(
        dplyr::ends_with("R2"),
        ~ pmax(., 0, na.rm = TRUE)
      )) %>%
      dplyr::mutate(dplyr::across(
        dplyr::ends_with("RMSE"),
        ~ pmin(., max(.), na.rm = TRUE)
      ))

    performance.summary <- c(
      performance.estimate %>% colMeans(),
      tryCatch(stats::t.test(
        performance.estimate %>%
          dplyr::pull(intra.RMSE),
        performance.estimate %>%
          dplyr::pull(multi.RMSE),
        alternative = "greater"
      )$p.value, error = function(e) {
        warning.message <- paste(
          "t-test of RMSE performance failed with error:",
          e$message
        )
        warning(warning.message)
        1
      }),
      tryCatch(stats::t.test(
        performance.estimate %>%
          dplyr::pull(intra.R2),
        performance.estimate %>%
          dplyr::pull(multi.R2),
        alternative = "less"
      )$p.value, error = function(e) {
        warning.message <- paste(
          "t-test of R2 performance failed with error:",
          e$message
        )
        warning(warning.message)
        1
      })
    )

    current.lock <- filelock::lock(perf.lock)
    # replace NaN p values with 1
    write(paste(target, paste(performance.summary %>% tidyr::replace_na(1), collapse = " ")),
      file = perf.file, append = TRUE
    )
    filelock::unlock(current.lock)

    return(target)
  }, ..., .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))

  return(normalized.results.folder)
}


#' Train sliding MISTy models
#'
#' Train local MISTy models by sliding a window across the sample as captured by the
#' view composition.
#'
#' @param views view composition.
#' @param positions a \code{data.frame}, \code{tibble} or a \code{matrix}
#'     with named coordinates in columns and rows for each spatial unit ordered
#'     as in the intraview.
#' @param window size of the
#' @param overlap overlap of consecutive windows (percentage).
#' @param results.folder path to the top level folder to store raw results.
#' @param minu minimum number of spatial units in the window.
#' @param ... all other parameters are passed to \code{\link{run_misty}()}.
#'
#' @return Path to the result folder(s) that can be passed to
#'     \code{\link{collect_results}()}.
#'
#' @examples
#' # TODO
#'
#' @export
run_sliding_misty <- function(views, positions, window, overlap = 50,
                              results.folder = "results", minu = 100,
                              ...) {
  x <- tibble::tibble(
    xl = seq(
      min(positions[, 1]),
      max(positions[, 1]),
      window - window * overlap / 100
    ),
    xu = xl + window
  ) %>%
    dplyr::filter(xl < max(positions[, 1])) %>%
    dplyr::mutate(xu = pmin(xu, max(positions[, 1]))) %>%
    round(2)

  y <- tibble::tibble(
    yl = seq(
      min(positions[, 2]),
      max(positions[, 2]),
      window - window * overlap / 100
    ),
    yu = yl + window
  ) %>%
    dplyr::filter(yl < max(positions[, 2])) %>%
    dplyr::mutate(yu = pmin(yu, max(positions[, 2]))) %>%
    round(2)

  tiles <- tidyr::expand_grid(x, y)

  # make nested plan here with sequential at the second level
  # retrieve current plan by simply running plan() without parameters
  old.plan <- future::plan()
  future::plan(list(old.plan, future::sequential))
  message("\nSliding")
  tiles %>% furrr::future_pwalk(\(xl, xu, yl, yu){
    selected.rows <- which(
      positions[, 1] >= xl & positions[, 1] <= xu &
        positions[, 2] >= yl & positions[, 2] <= yu
    )

    if (length(selected.rows) >= minu) {
      filtered.views <- views %>%
        filter_views(selected.rows)

      suppressMessages(run_misty(
        filtered.views,
        paste0(results.folder, "/", xl, "_", yl, "_", xu, "_", yu),
        ...
      ))
    }
  }, .progress = TRUE, .options = furrr::furrr_options(seed = TRUE))

  future::plan(old.plan)
  return(list.dirs(results.folder)[-1])
}
