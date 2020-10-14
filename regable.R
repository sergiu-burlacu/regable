# library("magrittr")
# library("dplyr")
# library("purrr")
# library("tidyr")
# library("tibble")
# library("broom")
# library("kableExtra")
# library("huxtable")
# library("stringr")
# source("/home/sergiu/Documents/Google Drive/R Codes/regable.R")
# df <- mtcars
# m1 <- df %>% lm(mpg ~ cyl, data = .)
# m2 <- df %>% lm(mpg ~ cyl + wt, data = .)

label_variables <- function(..., labels = NULL, fancy_latex = FALSE) {
  
  regs <- list(...)
  # Get the list of covariates
  covariates <-regs %>% 
    map_dfr(tidy) %>%
    select(term) %>%
    distinct()
  
  if (is.null(labels)) {
    labels <- c("No Labels" = "")
  }
  if (typeof(labels) == "character") {
    labels <- labels %>% tibble(term = ., label = names(.))
  }
  
  find_factors <- function(x) {
    select_if(x, is.factor)
  }
  # Label factors
  factor_labels <- regs %>%
    map(model.frame) %>% 
    map_dfr(find_factors)
  if (length(factor_labels) > 0) {
    factor_labels %<>% gather(var, value) %>%
      filter(!is.na(value)) %>%
      mutate(
        term = paste0(var, value),
        label = ifelse(str_detect(
          var, "^factor"),
          paste0(var, value),
          value
        )
      ) %>%
      select(-var, -value) %>%
      distinct(term, .keep_all = TRUE)
    labels <- rbind(labels, factor_labels)
  }
  
  # If label was not supplied, replace with the name of the variables
  cov_labels <- left_join(
    covariates, labels,
    by = "term"
  ) %>%
    mutate(
      label = ifelse(is.na(label), term, label)
    )
  
  if (fancy_latex == T) {
    sep_interact <- " $\\times$ "
  } else {
    sep_interact <- " x "
  }
  labels_interact <- cov_labels %>% filter(str_detect(term, ":"))
  labels_not_interact <- setdiff(cov_labels, labels_interact)
  labels_interact %<>% mutate(
    word_1 = word(term, 1, sep = "\\:"),
    word_2 = word(term, 2, sep = "\\:"),
    word_3 = word(term, 3, sep = "\\:"),
    word_4 = word(term, 4, sep = "\\:")
  )
  labels_interact <- labels_interact %>%
    select(-label) %>%
    gather(var, value, -term) %>%
    left_join(labels_not_interact %>% rename(value = term), by = "value") %>%
    filter(!is.na(value)) %>%
    group_by(term) %>%
    mutate(label2 = paste0(label, collapse = sep_interact)) %>%
    ungroup() %>%
    filter(var == "word_1") %>%
    select(term, label2) %>%
    rename(label = label2)
  cov_labels <- bind_rows(labels_not_interact, labels_interact)
  
  dep_var_labels <- regs %>%
    map(model.frame) %>%
    map(names) %>%
    map_chr(1) %>%
    tibble(term = .) %>%
    left_join(labels, by = "term") %>%
    mutate(term = ifelse(!is.na(label), label, term)) %>%
    pull(term)
  
  label_list = list(cov_labels = cov_labels, dep_var_labels = dep_var_labels)
  return(label_list)
  
}


regable_coef <- function(
                         ..., labels = NULL, type = "kable", format = "latex", 
                         intercept = FALSE, omit = NULL, keep = NULL,
                         star_levels = c(0.10, 0.05, 0.01),
                         beta_format = "%.2f", se_format = "%.2f",
                         header = c("numbers", "dep_vars"),
                         caption = NULL,
                         fancy_latex = FALSE,
                         longtable = T, booktabs = T,
                         align = NULL) {


  # TIDY OUTPUT ------------------------------------------------------
  # Create list of models
  regs <- list(...)

  # Get nr of models
  nmodels <- length(regs)

  # Tidy the model output and calculate nr of significance stars
  tidy_output <-
    # Tidies regressions lists indexing each with reg_number
    map_dfr(regs, tidy, .id = "model") %>%
    # Create numeric id for model
    as_tibble() %>%
    mutate(model = model %>% as.numeric()) %>%
    # Add nr of stars
    rowwise() %>%
    mutate(sig_stars = sum(p.value < star_levels)) %>%
    ungroup()

  #  DROP OR KEEP VARIABLES --------------------------------------------
  #' TO DO:
  #' - this will drop factors and interaction terms, should replace exact matching with
  #' something else PRIORITY 3
  #' - is intercept called the same in all models? PRIORITY 1
  if (intercept == FALSE) {
    tidy_output %<>% filter(term != "(Intercept)")
  } else {
    keep <- c(keep, "(Intercept)")
  }
  if (!is.null(keep)) {
    omit <- NULL
    keep <- paste(keep, collapse = "|")
    tidy_output %<>% filter(str_detect(term, keep))
  }
  if (!is.null(omit)) {
    omit <- paste(omit, collapse = "|")
    tidy_output %<>% filter(!str_detect(term, omit))
  }



  # Get the list of covariates
  covariates <- tidy_output %>%
    select(term) %>%
    distinct() %>%
    pull()

  label_list <- label_variables(..., labels = labels, 
                                fancy_latex = fancy_latex)
  cov_labels <- label_list[[1]]
  dep_var_labels <- label_list[[2]]
  
  # Header
  if (typeof(header) == "character"){
    dep_vars <- c(" ", dep_var_labels)
    numbers <- c(" ", paste0("(", seq_along(regs), ")"))
    header_list <- list(dep_vars = dep_vars, numbers = numbers)
  } 
  if (typeof(header) == "list") {
    header_list <- header
  }
  
  # PREPARE FINAL TABLE -----------------------------------------------------
  #' Expand grid works well because if a variable is missing in a model, it will
  #' still include that variables with missing value. Maybe it works without it,
  #' i'm ok as long as it's working
  coef_table <- expand.grid(
    term = covariates, model = 1:nmodels, stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    left_join(tidy_output, by = c("model", "term")) %>%
    mutate(beta = estimate, se = std.error)


  # If kable table then we need to format digits prior to sending it to kable
  #' If format is numeric, pass it to sprintf
  #' If it's a character leave it as it is. More flexibility
  #' Huxtable does this by itself
  if (type != "hux") {
    # If beta_format is numeric, pass it to sprintf
    if (typeof(beta_format) == "double") {
      beta_format <- beta_format %>% sprintf("%s%df", "%.", .)
    }
    # If se_format is numeric, pass it to sprintf
    if (typeof(se_format) == "double") {
      se_format <- se_format %>% sprintf("%s%df", "%.", .)
    }
    coef_table %<>% mutate(
      beta = ifelse(
        !is.na(beta),
        beta %>% sprintf(beta_format, .),
        ""),
      se = ifelse(
        !is.na(se),
        se %>% sprintf(se_format, .),
        "")
    )
  }

  if (format == "latex" & fancy_latex == TRUE) {
    # Add significance stars
    if (!is.null(star_levels)) {
      coef_table %<>% mutate(
        beta = case_when(
          sig_stars == 0 ~ paste0("$", beta, "$"),
          beta == "" ~ beta,
          TRUE ~ paste0("$", beta, "^{", strrep("*", sig_stars), "}$")
          )
        )
    } else {
      coef_table %<>% mutate(
        beta = ifelse(beta == "", beta, paste0("$", beta, "$"))
        )
      }
  } else {
    if (!is.null(star_levels)) {
      coef_table %<>% 
        mutate(beta = ifelse(is.na(beta), "", as.character(beta))) %>%
        mutate(
        beta = case_when(
          sig_stars == 0 | beta == "" ~ beta,
          TRUE ~ paste0(beta, strrep("*", sig_stars))
        )
      )
    } 
  }

  coef_table %<>% mutate(
    se = if_else(
      is.na(std.error), "", paste0("(", se, ")")
    )
  ) %>%
    select(model, beta, se, term)

  # Gather and spread the final table -------------------------------
  # 'order according to how they appeared in the equation
  # Need to try to order interaction terms after their terms somehow
  coef_table %<>%
    # Transpose to the standar coefficient, se under format
    gather(key = "beta_se", value = "value", -model, -term, factor_key = TRUE) %>%
    spread(key = model, value = value, fill = "", sep = "_") %>%
    mutate(order2 = row_number()) %>%
    # Add labels
    left_join(cov_labels %>%
      mutate(order1 = row_number()), by = "term") %>%
    # if the intercept is around, make it last
    # should assign an option to label intercept PRIORITY 3
    mutate(order0 = if_else(term == "(Intercept)", 2, 1)) %>%
    arrange(order0, order1, order2) %>%
    select(-order0, -order1, -order2) %>%
    mutate(term = case_when(!is.na(label) ~ label, TRUE ~ term)) %>%
    mutate(term = if_else(beta_se == "beta", term, "")) %>%
    select(-label) %>%
    rename(row_names = term)


  if (type == "hux") {
    coef_table <- as_hux(coef_table)
    coef_table %<>%
      set_number_format(coef_table$beta_se == "beta", starts_with("model"), beta_format) %>%
      set_number_format(coef_table$beta_se == "se", starts_with("model"), se_format) %>%
      select(-beta_se)

    if (format == "latex" & fancy_latex == TRUE) {
      escape_contents(coef_table)[, 2:ncol(coef_table)] <- FALSE
    }

    if (length(header) > 0) {
        if ((typeof(header) == "character")) {
          row_1 <- header_list[[header[[1]]]]
        }
        if ((typeof(header) == "list")) {
          row_1 <- header_list[[1]]
        }
        coef_table <- rbind(row_1, coef_table, copy_cell_props = FALSE)
        bottom_border(coef_table)[1, ] <- 0.5
        
        if (length(header) == 2) {
          if ((typeof(header) == "character")) {
            row_2 <- header_list[[header[[2]]]]
          }
          if ((typeof(header) == "list")) {
            row_2 <- header_list[[2]]
          }
        coef_table <- rbind(row_2, coef_table, copy_cell_props = FALSE)
        }
    }

    align(coef_table)[, 2:ncol(coef_table)] <- "center"
    return(coef_table)
    
    
  } else if (type == "kable") {
    if (fancy_latex == TRUE) {
      escape <- FALSE
    } else {
      escape <- TRUE
    }
    coef_table %<>% select(-beta_se)
    
    if (length(header) > 0) {
      if ((typeof(header) == "character")) {
        row_1 <- header_list[[header[[1]]]]
      }
      if ((typeof(header) == "list")) {
        row_1 <- header_list[[1]]
      }
      coef_table %<>% rename_all(~row_1)
    }

    if (is.null(align)) {
      align = c("l", rep("c", length(regs)))
    }
    coef_table %<>%
      kable(coef_table,
        format = format, longtable = longtable, booktabs = booktabs,
        caption = caption, align = align,
        linesep = "", escape = escape
      )
    if (length(header) == 2) {
      if ((typeof(header) == "character")) {
        row_2 <- header_list[[header[[2]]]]
      }
      if ((typeof(header) == "list")) {
        row_2 <- header_list[[2]]
      }
      coef_table %<>% add_header_above(row_2, line = F)
    }
    return(coef_table)
  } else {
    return(coef_table %>% select(-beta_se))
  }
}

regable_stat <- function(..., type = "kable", format = "latex", stats = c("Observations" = "nobs", "Adjusted R2" = "adj.r.squared"),
                         stat_format = 2, add_stats = NULL,
                         fancy_latex = FALSE) {
  # Default Parameters
  # type = "kable"
  # format = "html"
  # stats = c("Observations" =  "nobs", "Adjusted R2" = "adj.r.squared")
  # stat_format = 2
  # regs <- list(model1, model2)
  # add_stats <- list(Controls = c("Yes", "Yes"), `P-value` = as_integer(c(1,2)))
  regs <- list(...)
  # Get nr of models
  nmodels <- length(regs)

  # Extra statistics
  if (typeof(add_stats) == "list") {
    add_stats %<>% as_tibble()
  }

  if (typeof(stat_format) == "double") {
    stat_format <- stat_format %>% sprintf("%s%df", "%.", .)
  }

  regs_glance <- regs %>% map_dfr(glance)
  stat_table <-
    regs_glance %>%
    mutate(nobs = map_int(regs, nobs, use.fallback = TRUE)) %>%
    select(stats)
  if (!is.null(add_stats)) {
    stat_table <- cbind(add_stats, stat_table)
  }
  stat_table_class <- map_chr(stat_table, class)

  if (type != "hux") {
    stat_table %<>% 
      mutate_at(
        vars(which(stat_table_class == "numeric")), ~sprintf(stat_format, .)
        ) %>%
      mutate_at(
        vars(which(stat_table_class != "numeric")), 
        ~as.character(.)  
        )

  }
  stat_table %<>%
    mutate(model = row_number()) %>%
    select(model, everything()) %>%
    gather(stat, value, -model, factor_key = TRUE) %>%
    spread(key = model, value = value, fill = "", sep = "_") %>%
    rename(row_names = stat)

  if (type == "hux") {
    stat_table <- as_hux(stat_table)
    number_format(stat_table)[ 1:nrow(stat_table), 2:ncol(stat_table)] <- stat_format
    number_format(stat_table)[stat_table_class == "integer"] <- 0
    align(stat_table)[, 2:ncol(stat_table)] <- "center"
    stat_table
  } else if (type == "kable") {
    col_numbers <- c(" ", paste0("(", seq_along(regs), ")"))
    stat_table %<>% rename_all(~col_numbers)
    stat_table %<>%
      kable(
        format = format, longtable = T, booktabs = T,
        align = c("l", rep("c", length(regs))),
        linesep = ""
      )
    return(stat_table)
  } else {
    return(stat_table)
  }
}



regable <- function(..., type = "kable", format = "latex", labels = NULL,
                    intercept = FALSE, omit = NULL, keep = NULL,
                    star_levels = c(0.10, 0.05, 0.01),
                    beta_format = 2, se_format = 2, caption = NULL,
                    stats = c("Observations" = "nobs", "Adjusted R2" = "adj.r.squared"),
                    stat_format = 2, add_stats = NULL,
                    header = c("numbers", "dep_vars"),
                    fancy_latex = FALSE,
                    longtable = T, booktabs = T,
                    align = NULL) {
  # type = "kable"
  # format = "latex"
  # intercept = FALSE
  # omit = NULL
  # keep = NULL
  # labels = NULL
  # star_levels = c(0.10, 0.05, 0.01)
  # beta_format = 2
  # stat_format = 2
  # se_format = 2
  # caption = NULL
  # stats = c("Observations" =  "nobs", "Adjusted R2" = "adj.r.squared")
  # add_stats = list(Controls = c("Yes", "Yes"), `P-value` = as_integer(c(1,2)))
  # fancy_latex = FALSE
  # labels = labels = c("Miles per galon" = "mpg", "Weight" = "wt", "Cylinder" = "cyl")
  if (fancy_latex == TRUE) {
    escape <- FALSE
  } else {
    escape <- TRUE
  }
  if (type == "kable") {
    type <- "skip_kable"
  }
  regs <- list(...)
  coef_table <- regable_coef(
    ...,
    type = type, format = format, labels = labels, intercept = intercept, omit = omit, keep = keep, star_levels = star_levels,
    beta_format = beta_format, se_format = se_format,
    header = header,
    fancy_latex = fancy_latex,
    longtable = longtable, booktabs = booktabs,
    align = align
  )
  stat_table <- regable_stat(
    ...,
    type = type, format = format, stats = stats,
    stat_format = stat_format, add_stats = add_stats
  )
  

  if (type == "hux") {
    bottom_border(coef_table)[nrow(coef_table), ] <- 0.5
    final_table <- rbind(coef_table, stat_table)
    top_border(final_table)[1, ] <- 1
    bottom_border(final_table)[nrow(final_table), ] <- 1
    return(final_table)
  } else if (type == "skip_kable") {
    nrow_coef <- nrow(coef_table)
    nrow_stat <- nrow(stat_table)
    final_table <- bind_rows(
      coef_table %>% mutate(row_names = as.character(row_names)),
      stat_table %>% mutate(row_names = as.character(row_names))
    )
    #
    label_list <- label_variables(..., labels = labels, fancy_latex = fancy_latex)
    dep_var_labels <- label_list[[2]]
    # Header
    if (typeof(header) == "character"){
      dep_vars <- c(" ", dep_var_labels)
      numbers <- c(" ", paste0("(", seq_along(regs), ")"))
      header_list <- list(dep_vars = dep_vars, numbers = numbers)
    } 
    if (typeof(header) == "list") {
      header_list <- header
    }
    
    if (length(header) > 0) {
      if ((typeof(header) == "character")) {
        row_1 <- header_list[[header[[1]]]]
      }
      if ((typeof(header) == "list")) {
        row_1 <- header_list[[1]]
      }
      final_table %<>% rename_all(~row_1)
    }
    
    if (is.null(align)) {
      align = c("l", rep("c", length(regs)))
    }
    final_table %<>%
      kable(final_table,
            format = format, longtable = longtable, booktabs = booktabs,
            caption = caption, align = align,
            linesep = "", escape = escape
      ) %>% row_spec(nrow_coef, hline_after = T)
    if (length(header) == 2) {
      if ((typeof(header) == "character")) {
        row_2 <- header_list[[header[[2]]]]
      }
      if ((typeof(header) == "list")) {
        row_2 <- header_list[[2]]
      }
      final_table %<>% add_header_above(row_2, line = F)
    }
    return(final_table)
  } else {
    final_table <- bind_rows(
      coef_table %>% mutate(row_names = as.character(row_names)),
      stat_table %>% mutate(row_names = as.character(row_names))
    )
    return(final_table)
  }
}
