boot_ddsc_ml_fixef <-
  function(model,
           nsim,
           seed,
           results,
           predictor,
           moderator,
           moderator_values,
           scaling_SDs,
           descriptives,
           level) {
    # run bootMer
    booted_fixef <-
      lme4::bootMer(
        x = model,
        FUN = lme4::fixef,
        nsim = 1000,
        use.u = FALSE,
        seed = seed,
        type = c("parametric"),
        verbose = FALSE
      )
    # obtain all the bootstrap estimates
    boot_est <- data.frame(booted_fixef$t)
    # rename the columns in the bootstrap datafile

    names(boot_est)[names(boot_est) %in% grep("Intercept", names(boot_est), value = TRUE)] <- "b0"
    names(boot_est)[names(boot_est) == predictor] <- "b1"
    names(boot_est)[names(boot_est) == moderator] <- "b2"
    names(boot_est)[names(boot_est) %in% grep(predictor, names(boot_est), value = TRUE) &
      names(boot_est) %in% grep(moderator, names(boot_est), value = TRUE)] <- "b3"
    # calculate the derivative estimates
    boot_est[, "w_11"] <- boot_est[, "b1"] + moderator_values[1] * boot_est[, "b3"]
    boot_est[, "w_21"] <- boot_est[, "b1"] + moderator_values[2] * boot_est[, "b3"]
    boot_est[, "r_xy1"] <- boot_est[, "w_11"] / (1 * scaling_SDs["SD_y1"])
    boot_est[, "r_xy2"] <- boot_est[, "w_21"] / (1 * scaling_SDs["SD_y2"])
    boot_est[, "b_11"] <- boot_est[, "w_11"] / (1 * scaling_SDs["SD_pooled"])
    boot_est[, "b_21"] <- boot_est[, "w_21"] / (1 * scaling_SDs["SD_pooled"])

    boot_est[, "r_xy1y2"] <-
      (boot_est[, "r_xy1"] * scaling_SDs["SD_y1"] - boot_est[, "r_xy2"] * scaling_SDs["SD_y2"]) /
        sqrt(scaling_SDs["SD_y1"]^2 + scaling_SDs["SD_y2"]^2 -
          2 * descriptives["means_y1", "means_y2"] * scaling_SDs["SD_y1"] * scaling_SDs["SD_y2"])

    boot_est[, "main_effect"] <- boot_est[, "b1"]
    boot_est[, "moderator_effect"] <- boot_est[, "b2"]
    boot_est[, "interaction"] <- boot_est[, "b3"]

    boot_est[, "q_b11_b21"] <- atanh(boot_est[, "b_11"]) - atanh(boot_est[, "b_21"])
    boot_est[, "q_rxy1_rxy2"] <- atanh(boot_est[, "r_xy1"]) - atanh(boot_est[, "r_xy2"])

    boot_est[, "cross_over_point"] <- (-1) * boot_est[, "b2"] / boot_est[, "b3"]

    boot_est[, "interaction_vs_main"] <- abs(boot_est[, "b3"]) - abs(boot_est[, "b1"])
    boot_est[, "interaction_vs_main_bscale"] <- abs(boot_est[, "b_11"] - boot_est[, "b_21"]) - abs((boot_est[, "b_11"] + boot_est[, "b_21"]) / 2)
    boot_est[, "interaction_vs_main_rscale"] <- abs(boot_est[, "r_xy1"] - boot_est[, "r_xy2"]) - abs((boot_est[, "r_xy1"] + boot_est[, "r_xy2"]) / 2)

    boot_est[, "dadas"] <- abs(boot_est[, "w_11"] - (boot_est[, "w_21"])) - abs(boot_est[, "w_11"] + (boot_est[, "w_21"]))
    boot_est[, "dadas_bscale"] <- abs(boot_est[, "b_11"] - (boot_est[, "b_21"])) - abs(boot_est[, "b_11"] + (boot_est[, "b_21"]))
    boot_est[, "dadas_rscale"] <- abs(boot_est[, "r_xy1"] - (boot_est[, "r_xy2"])) - abs(boot_est[, "r_xy1"] + (boot_est[, "r_xy2"]))

    boot_est[, "abs_diff"] <- abs(boot_est[, "w_11"] - (boot_est[, "w_21"]))
    boot_est[, "abs_sum"] <- abs(boot_est[, "w_11"] + (boot_est[, "w_21"]))

    boot_est[, "abs_diff_bscale"] <- abs(boot_est[, "b_11"] - (boot_est[, "b_21"]))
    boot_est[, "abs_sum_bscale"] <- abs(boot_est[, "b_11"] + (boot_est[, "b_21"]))

    boot_est[, "abs_diff_rscale"] <- abs(boot_est[, "r_xy1"] - (boot_est[, "r_xy2"]))
    boot_est[, "abs_sum_rscale"] <- abs(boot_est[, "r_xy1"] + (boot_est[, "r_xy2"]))

    # Calculate bootstrap summary statistics
    boot_results <- t(as.data.frame(sapply(
      boot_est,
      function(x) {
        c(
          Estimate = mean(x, na.rm = TRUE),
          SE = stats::sd(x, na.rm = TRUE),
          stats::quantile(x, c((1 - level) / 2, 1 - (1 - level) / 2), na.rm = TRUE)
        )
      }
    )))

    # Add column names
    colnames(boot_results) <- c("boot_est", "boot_se", "boot_LL", "boot_UL")
    # filter to same parameter estimates as in normal results
    boot_results <- boot_results[rownames(boot_results) %in% rownames(results), ]
    boot_results <- boot_results[match(rownames(results), rownames(boot_results)), ]

    results_with_boots <- cbind(results, boot_results)
    return(results_with_boots)
  }
