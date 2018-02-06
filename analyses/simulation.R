### Setup R session ------------------------------------------------------------
options(
  repos = c("https://cran.rstudio.com/", "http://rweb.crmda.ku.edu/kran"),
  scipen = 999
)

needed_packages <- c("tidyverse", "portableParallelSeeds", "MplusAutomation",
  "tidyselect", "stringr", "forcats", "lubridate")
load_packages <- function(x) {
  if (!(x %in% rownames(installed.packages()))) install.packages(x)
  
  suppressPackageStartupMessages(require(x, character.only = TRUE))
}
sapply(needed_packages, load_packages)

rm(list = ls(all.names = TRUE)); gc()


### Utility functions-----------------------------------------------------------
safe_mplus <- safely(runModels)
quiet_read <- quietly(readModels)
only_if <- function(condition) {
  function(func) {
    if (condition) {
      func
    } else {
      function(., ...) .
    }
  }
}
logit <- function(x) {
  log(x / (1- x))
}
invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}
check_length <- function(x) {
  if (nchar(x) > 80) {
    breaks <- seq(80, nchar(x) + 80, 80)
    for (i in seq_along(breaks)) {
      cur_break <- breaks[i]
      if (nchar(x) > cur_break) {
        x <- strsplit(x, "") %>% unlist()
        replace <- which(x %in% c(" ", "+"))
        replace <- max(replace[which(replace < cur_break)])
        x[replace] <- paste(x[replace], "\n", sep = "")
        x <- paste(x, collapse = "")
      }
    }
    return(x)
  } else {
    return(x)
  }
}
getSection <- function(sectionHeader, outfiletext, headers = "standard",
  omit = NULL) {
  #encode the top-level major headers here, but allow for custom headers to be passed in
  #omit allows for one or more strings from headers not to be considered
  #just used for factor score statistics at the moment (these include a SAMPLE STATISTICS section)
  if (headers[1L] == "standard") headers <- c("INPUT INSTRUCTIONS", "SUMMARY OF ANALYSIS",
    "SUMMARY OF DATA FOR THE FIRST DATA SET", "SUMMARY OF DATA FOR THE FIRST REPLICATION",
    "SUMMARY OF MISSING DATA PATTERNS FOR THE FIRST REPLICATION",
    "SUMMARY OF MISSING DATA PATTERNS",
    "COVARIANCE COVERAGE OF DATA FOR THE FIRST REPLICATION",
    "SAMPLE STATISTICS", "SAMPLE STATISTICS FOR THE FIRST REPLICATION",
    "CROSSTABS FOR CATEGORICAL VARIABLES", "UNIVARIATE PROPORTIONS AND COUNTS FOR CATEGORICAL VARIABLES",
    "RANDOM STARTS RESULTS RANKED FROM THE BEST TO THE WORST LOGLIKELIHOOD VALUES",
    "TESTS OF MODEL FIT", "MODEL FIT INFORMATION", "CLASSIFICATION QUALITY",
    "FINAL CLASS COUNTS AND PROPORTIONS FOR THE LATENT CLASSES",
    "FINAL CLASS COUNTS AND PROPORTIONS FOR THE LATENT CLASS PATTERNS",
    "LATENT TRANSITION PROBABILITIES BASED ON THE ESTIMATED MODEL",
    "FINAL CLASS COUNTS AND PROPORTIONS FOR EACH LATENT CLASS VARIABLE",
    "CLASSIFICATION OF INDIVIDUALS BASED ON THEIR MOST LIKELY LATENT CLASS MEMBERSHIP",
    "Average Latent Class Probabilities for Most Likely Latent Class Membership \\(Row\\)",
    "Classification Probabilities for the Most Likely Latent Class Membership \\(Row\\)",
    "Classification Probabilities for the Most Likely Latent Class Membership \\(Column\\)",
    "Logits for the Classification Probabilities for the Most Likely Latent Class Membership \\(Row\\)",
    "Logits for the Classification Probabilities for the Most Likely Latent Class Membership \\(Column\\)",
    "MODEL RESULTS", "LOGISTIC REGRESSION ODDS RATIO RESULTS", "RESULTS IN PROBABILITY SCALE",
    "IRT PARAMETERIZATION IN TWO-PARAMETER LOGISTIC METRIC",
    "IRT PARAMETERIZATION IN TWO-PARAMETER PROBIT METRIC",
    "ALTERNATIVE PARAMETERIZATIONS FOR THE CATEGORICAL LATENT VARIABLE REGRESSION",
    "LATENT CLASS ODDS RATIO RESULTS", "LOGRANK OUTPUT", "STANDARDIZED MODEL RESULTS",
    "R-SQUARE", "QUALITY OF NUMERICAL RESULTS", "TECHNICAL OUTPUT", "TECHNICAL \\d+ OUTPUT",
    "TECHNICAL 5/6 OUTPUT",
    "TOTAL, TOTAL INDIRECT, SPECIFIC INDIRECT, AND DIRECT EFFECTS",
    "STANDARDIZED TOTAL, TOTAL INDIRECT, SPECIFIC INDIRECT, AND DIRECT EFFECTS", "CONFIDENCE INTERVALS OF MODEL RESULTS",
    "CONFIDENCE INTERVALS FOR THE LOGISTIC REGRESSION ODDS RATIO RESULTS",
    "CREDIBILITY INTERVALS OF MODEL RESULTS",
    "CONFIDENCE INTERVALS OF STANDARDIZED MODEL RESULTS",
    "CREDIBILITY INTERVALS OF STANDARDIZED MODEL RESULTS",
    "CONFIDENCE INTERVALS IN PROBABILITY SCALE",
    "CONFIDENCE INTERVALS OF TOTAL, TOTAL INDIRECT, SPECIFIC INDIRECT, AND DIRECT EFFECTS",
    "CONFIDENCE INTERVALS OF STANDARDIZED TOTAL, TOTAL INDIRECT, SPECIFIC INDIRECT,", #omitted "AND DIRECT EFFECTS"
    "EQUALITY TESTS OF MEANS ACROSS CLASSES USING POSTERIOR PROBABILITY-BASED",
    "EQUALITY TESTS OF MEANS ACROSS CLASSES USING THE BCH PROCEDURE",
    "EQUALITY TESTS OF MEANS ACROSS CLASSES USING THE 3-STEP PROCEDURE",
    "EQUALITY TESTS OF MEANS/PROBABILITIES ACROSS CLASSES",
    "THE FOLLOWING DATA SET\\(S\\) DID NOT RESULT IN A COMPLETED REPLICATION:",
    "RESIDUAL OUTPUT", "MODEL MODIFICATION INDICES", "MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES",
    "Available post-processing tools:",
    "FACTOR SCORE INFORMATION \\(COMPLETE DATA\\)", "SUMMARY OF FACTOR SCORES", "PLOT INFORMATION", "SAVEDATA INFORMATION",
    "RESULTS SAVING INFORMATION", "SAMPLE STATISTICS FOR ESTIMATED FACTOR SCORES", "DIAGRAM INFORMATION",
    "Beginning Time:\\s*\\d+:\\d+:\\d+", "MUTHEN & MUTHEN"
  )
  
  if (!is.null(omit)) headers <- headers[which(!headers %in% omit)] #drop omit
  
  #allow for syntax to include :: to specify a header that spans 2 rows. Example:
  #FINAL CLASS COUNTS AND PROPORTIONS FOR THE LATENT CLASSES
  #BASED ON THE ESTIMATED MODEL
  
  #note that this will identify a unique match for the target sectionHeader, but the search for
  #subsequent headers just uses the one-row list above. For now, this works in all cases I know of.
  if (grepl("::", sectionHeader, fixed=TRUE)) {
    firstLine <- sub("^(.*)::.*$", "\\1", sectionHeader, perl=TRUE)
    nextLine <- sub("^.*::(.*)$", "\\1", sectionHeader, perl=TRUE)
    candidates <- grep(firstLine, outfiletext, perl=TRUE) #will return NA if no match
    bothMatch <- grep(nextLine, outfiletext[candidates+1], perl=TRUE)[1] #return first match among candidates
    if (!is.na(bothMatch)) { beginSection <- candidates[bothMatch] + 1 #since it's a two-line header, skip the first to match typical case
    } else { beginSection <- NA } #could not find section with both rows
  } else {
    beginSection <- grep(sectionHeader, outfiletext, perl=TRUE)[1]  
  }
  
  #if section header cannot be found, then bail out
  if (is.na(beginSection)) return(NULL)
  
  #form alternation pattern for regular expression (currently adds leading and trailing spaces permission to each header)
  headerRegexpr <- paste("(", paste(gsub("(.*)", "^\\\\s*\\1\\\\s*$", headers, perl=TRUE), sep="", collapse="|"), ")", sep="")
  headerLines <- grep(headerRegexpr, outfiletext, perl=TRUE)
  subsequentHeaders <- which(headerLines > beginSection)
  
  if (length(subsequentHeaders) == 0) nextHeader <- length(outfiletext) #just return the whole enchilada
  else nextHeader <- headerLines[subsequentHeaders[1]] - 1
  
  section.found <- outfiletext[(beginSection+1):nextHeader]
  attr(section.found, "lines") <- beginSection:nextHeader
  
  return(section.found)
  
}
extractValue <- function(pattern, textToScan, filename, type = "int") {
  #regex pattern now allows for specification to search for value on some line before or after match
  #example: +2:the Observed and the Replicated Chi-Square Values
  
  offset <- 0
  if (grepl("^[+-]+\\d+:.*$", pattern, perl=TRUE)) {
    offset <- as.numeric(sub("^([+-]+\\d+):.*$", "\\1", pattern, perl=TRUE))
    pattern <- sub("^[+-]+\\d+:(.*)$", "\\1", pattern, perl=TRUE) #chop offset
  }
  
  #locate the matching line in the output file
  matchpos <- grep(pattern, textToScan, ignore.case=TRUE)
  matchlines <- textToScan[(matchpos+offset)]
  
  if (length(matchlines) > 1) {
    return(length(matchlines))
    #stop("More than one match found for parameter: ", pattern, "\n  ", filename)
    #return(matchlines) #not sure what I was thinking here... seems better to stop than warn and return lines
  }
  else if (length(matchlines) == 0) {
    #if the parameter of interest not found in this file, then return NA
    #warning(paste("Parameter not found: ", pattern, "\n  ", filename, sep=""))
    if (type == "int") return(NA_integer_)
    else if (type == "dec") return(NA_real_)
    else if (type == "str") return(NA_character_)
  }
  
  #different idea: concatenate pattern with var type and match on that
  #then sub just the pattern part from the larger line
  
  typePrefix <- substr(type, 1, 3)
  
  if (typePrefix == "int") {
    regexp <- "-*\\d+" #optional negative sign in front
  }
  else if (typePrefix == "dec") {
    #regexpr: -*\\d+\\.\\d+ : -* optional negative sign, \\d+ match at least one digit \\. match decimal sign \\d+ match decimal digits
    regexp <- "-*\\d+\\.\\d+"
  }
  else if (typePrefix == "str") {
    regexp <- paste(pattern, ".*", sep="")
  }
  
  #locate the match
  valueMatches <- gregexpr(regexp, matchlines[1], perl=TRUE)[[1]]
  
  if (type == "str") {
    #remove the tag portion of the string (e.g., "title:"), retaining rest of line
    returnVal <- as.character(sub(pattern, "", matchlines[1], ignore.case=TRUE))
  }
  else {
    #excessively tight syntax: replace dec[15] with 15, if number at end of type. Otherwise return just "dec".
    #then grep result for only numeric characters (\\d+). If grep is false (i.e., no numerals in substitution,
    #then no index was specified in type, so type must be simply "dec", "int", or "str" (as opposed to "int[15]"), so set as 1
    if (!grepl("^\\d+$", whichMatch <- sub("^.*\\[(\\d+)\\]$", "\\1", type, perl=TRUE), perl=TRUE)) whichMatch <- 1
    else whichMatch <- as.numeric(whichMatch)
    
    #pull from the start of the match through match.length, which is the length of characters that matched
    #need to subtract one from the start + length offset to grab the correct number of characters
    #(e.g., if the match runs from 40-44, the start will be 40, with length 5, but 40 + 5 would be 6 characters, hence -1
    returnVal <- as.numeric(substr(matchlines[1], valueMatches[whichMatch], valueMatches[whichMatch] + attr(valueMatches, "match.length")[whichMatch] - 1))
    
  }
  
  return(returnVal)
}
master_round <- function(x, cut = 0.5) {
  ifelse(x >= cut, 1, 0)
}
calc_kappa <- function(x, y, min_score, max_score) {
  score1 <- x
  score2 <- y
  
  if(length(unique(score1)) == 1 && length(unique(score2)) == 1 &&
      unique(score1) == unique(score2)) {
    return(1)
  }
  
  if (missing(min_score)) {
    min_score <- min(min(score1),min(score2))
  }
  if (missing(max_score)) {
    max_score <- max(max(score1),max(score2))
  }
  
  score1 <- factor(score1, levels = min_score:max_score)
  score2 <- factor(score2, levels = min_score:max_score)
  
  # pairwise frequencies
  confusion.mat <- table(data.frame(score1, score2))
  confusion.mat <- confusion.mat / sum(confusion.mat)
  
  # get expected pairwise frequencies under independence
  histogram.a <- table(score1) / length(table(score1))
  histogram.b <- table(score2) / length(table(score2))
  expected.mat <- histogram.a %*% t(histogram.b)
  expected.mat <- expected.mat / sum(expected.mat)
  
  # get weights
  labels <- as.numeric(as.vector(names(table(score1))))
  if (min_score == 1 && max_score == 8) {
    profs <- c("000", "100", "010", "001", "110", "101", "011", "111")
    weights <- matrix(data = NA, nrow = 8, ncol = 8)
    for (r in seq_along(profs)) {
      for (c in seq_along(profs)) {
        row <- profs[r] %>% str_split("") %>% flatten_chr()
        col <- profs[c] %>% str_split("") %>% flatten_chr()
        
        weights[r, c] <- length(which(!(row == col)))^2
      }
    }
  } else if (min_score == 1 && max_score == 16) {
    profs <- c("0000", "1000", "0100", "0010", "0001", "1100", "1010", "1001",
      "0110", "0101", "0011", "1110", "1101", "1011", "0111", "1111")
    weights <- matrix(data = NA, nrow = 16, ncol = 16)
    for (r in seq_along(profs)) {
      for (c in seq_along(profs)) {
        row <- profs[r] %>% str_split("") %>% flatten_chr()
        col <- profs[c] %>% str_split("") %>% flatten_chr()
        
        weights[r, c] <- length(which(!(row == col)))^2
      }
    }
  } else {
    weights <- outer(labels, labels, FUN = function(x,y) (x-y)^2)
  }
  
  # calculate kappa
  kappa <- 1 - sum(weights * confusion.mat) / sum(weights * expected.mat)
  kappa
}
calc_mean <- function (values, na.rm = TRUE) {
  max999 <- function(x) sign(x) * min(0.999, abs(x))
  min001 <- function(x) sign(x) * max(0.001, abs(x))
  values <- map_dbl(.x = values, .f = max999)
  values <- map_dbl(.x = values, .f = min001)
  
  if (na.rm) {
    values <- values[which(!is.na(values))]
  }
  
  r2z <- function(x) 0.5 * log((1 + x) / (1 - x))
  z2r <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)
  values <- map_dbl(.x = values, .f = r2z)
  values <- mean(values)
  values <- z2r(values)
  values
}
measure_fit <- function(x) {
  x %>%
    gather(fit_measure, value, -(model:dataset)) %>%
    filter(!is.na(value)) %>%
    group_by(dataset, fit_measure) %>%
    top_n(-1, wt = value) %>%
    group_by(fit_measure, model, dataset) %>%
    mutate(measure_picks = n()) %>%
    group_by(dataset, fit_measure) %>%
    mutate(dataset_picks = n()) %>%
    ungroup() %>%
    mutate(
      weight = measure_picks / dataset_picks,
      total_reps = max(dataset),
      pct = weight / total_reps
    ) %>%
    select(model, fit_measure, pct)
}


### Simulation functions -------------------------------------------------------
gen_strc <- function(num_resp, num_att) {
  num_param <- (2^num_att) - 1
  params <- rnorm(num_param, mean = 0, sd = 2)
  
  param_names <- seq_len(num_att) %>%
    map(., .f = function(num, num_att) {
      combs <- combn(seq_len(num_att), num)
      comb_names <- list_along(seq_len(ncol(combs)))
      for (i in seq_along(comb_names)) {
        comb_names[i] <- paste0("G_", num, paste(combs[,i], collapse = ""))
      }
      flatten_chr(comb_names)
    }, num_att = num_att) %>%
    flatten_chr()
  
  prof_list <- data_frame(att1 = c(0, 1))
  for (i in seq_len(num_att)) {
    prof_list[,paste0("att", i)] <- c(0, 1)
  }
  att_profiles <- expand.grid(prof_list) %>%
    mutate(row_sum = rowSums(.)) %>%
    group_by(row_sum) %>%
    nest(.key = profiles) %>%
    mutate(new_profiles = map(profiles, .f = function(x) {
      arrange_all(x, desc)
    })) %>%
    select(-profiles) %>%
    unnest() %>%
    select(-row_sum)
  
  if (num_att == 3) {
    att_profiles <- att_profiles %>%
      mutate(
        sum = att1 * params[1] + att2 * params[2] + att3 * params[3] +
          att1 * att2 * params[4] + att1 * att3 * params[5] +
          att2 * att3 * params[6] + att1 * att2 * att3 * params[7],
        exp_sum = exp(sum),
        prob = exp_sum / sum(exp_sum)
      )
    resp_profiles <- att_profiles %>%  
      sample_n(size = num_resp, replace = TRUE, weight = prob) %>%
      select(-(sum:prob)) %>%
      rowid_to_column(var = "respondent")
  } else if (num_att == 4) {
    att_profiles <- att_profiles %>%
      mutate(
        sum = att1 * params[1] + att2 * params[2] + att3 * params[3] +
          att4 * params[4] + att1 * att2 * params[5] + att1 * att3 * params[6] +
          att1 * att4 * params[7] + att2 * att3 * params[8] +
          att2 * att4 * params[9] + att3 * att4 * params[10] +
          att1 * att2 * att3 * params[11] + att1 * att2 * att4 * params[12] +
          att1 * att3 * att4 * params[13] + att2 * att3 * att4 * params[14] +
          att1 * att2 * att3 * att4 * params[15],
        exp_sum = exp(sum),
        prob = exp_sum / sum(exp_sum)
      )
    resp_profiles <- att_profiles %>%  
      sample_n(size = num_resp, replace = TRUE, weight = prob) %>%
      select(-(sum:prob)) %>%
      rowid_to_column(var = "respondent")
  }
  
  strc_params <- data_frame(
    parameter = param_names,
    value = params
  )
  
  list(
    strc_params = strc_params,
    strc_probs = att_profiles,
    resp_profiles = resp_profiles
  )
}
gen_meas <- function(num_att) {
  # generate q-matrix
  prof_list <- data_frame(att1 = c(0, 1))
  for (i in seq_len(num_att)) {
    prof_list[,paste0("att", i)] <- c(0, 1)
  }
  all_q <- expand.grid(prof_list) %>%
    mutate(row_sum = rowSums(.)) %>%
    filter(row_sum %in% c(1, 2)) %>%
    group_by(row_sum) %>%
    nest(.key = profiles) %>%
    mutate(new_profiles = map(profiles, .f = function(x) {
      arrange_all(x, desc)
    })) %>%
    select(-profiles) %>%
    unnest()
  reject <- TRUE
  while (reject) {
    q_matrix <- bind_rows(
      all_q %>%
        filter(row_sum == 1) %>%
        sample_n(size = 20, replace = TRUE) %>%
        select(-row_sum),
      all_q %>%
        sample_n(size = 10, replace = TRUE) %>%
        select(-row_sum)
    ) %>%
      rowid_to_column(var = "item")
    
    if (min(colSums(q_matrix)) >= 3) reject <- FALSE
  }
  
  # generate measurement model parameters
  intercepts <- data_frame(
    item = seq_len(30),
    parameter = paste0("L", seq_len(30), "_0"),
    type = "intercept"
  ) %>%
    mutate(value = runif(n = nrow(.), min = -3, max = 3))
  main_effects <- q_matrix %>%
    gather(key = "attribute", value = "measure", att1:last_col()) %>%
    filter(measure == 1) %>%
    mutate(
      attribute = str_replace(attribute, "att", ""),
      parameter = paste0("L", item, "_1", attribute),
      type = "main_effect"
    ) %>%
    mutate(value = runif(n = nrow(.), min = 0, max = 5)) %>%
    select(item, parameter, type, value)
  interactions <- q_matrix %>%
    mutate(total_att = q_matrix %>% select(contains("att")) %>% rowSums()) %>%
    filter(total_att == 2) %>%
    left_join(main_effects, by = "item") %>%
    select(item, parameter, value) %>%
    group_by(item) %>%
    nest() %>%
    mutate(inter_param = map(data, function(x) {
      me <- pull(x, parameter) %>%
        strsplit(., "_") %>% unlist()
      item <- str_subset(me, "L") %>% unique()
      attributes <- me[!(str_detect(me, "L"))] %>%
        str_replace("1", "") %>%
        paste(collapse = "")
      
      data_frame(
        parameter = paste0(item, "_2", attributes),
        type = "interaction",
        value = runif(n = 1, min = (-1 * min(x$value)), max = 2)
      )
    })) %>%
    select(-data) %>%
    unnest()

  list(
    q_matrix = q_matrix,
    meas_params = bind_rows(intercepts, main_effects, interactions)
  )
}
gen_data <- function(num_resp, num_att) {
  strc <- gen_strc(num_resp = num_resp, num_att = num_att)
  meas <- gen_meas(num_att = num_att)
  
  item_parameters <- meas$meas_params %>%
    group_by(item) %>%
    nest()
  
  good_data <- FALSE
  while (!good_data) {
    response_data <- strc$resp_profiles %>%
      group_by(respondent) %>%
      nest(.key = "prof") %>%
      mutate(resp_data = map(prof, function(x, item_parameters) {
        item_parameters %>%
          mutate(
            log_odds = map_dbl(data, function(params, prof) {
              master <- prof %>%
                gather(key = "attribute", value = "master") %>%
                mutate(
                  attribute = as.numeric(str_replace(attribute, "att", ""))
                )
              params %>%
                mutate(
                  include = map_lgl(parameter, function(p, master) {
                    att <- str_split(p, "_") %>% flatten_chr() %>%
                      .[2]
                    if (att == "0") return(TRUE)
                    att <- str_sub(att, start = 2, end = -1) %>%
                      str_split(pattern = "") %>% flatten_chr() %>% as.numeric()
                    if (all(master$master[which(master$attribute %in% att)] == 1)) {
                      TRUE
                    } else {
                      FALSE
                    }
                  },
                    master = master)
                ) %>%
                filter(include) %>%
                pull(value) %>%
                sum()
            }, prof = x),
            prob = invlogit(log_odds),
            random = runif(n = nrow(.), min = 0, max = 1),
            resp = case_when(
              random <= prob ~ 1,
              TRUE ~ 0
            ),
            item = paste0("X", item)
          ) %>%
          select(item, resp)
      },
        item_parameters = item_parameters)) %>%
      select(-prof) %>%
      unnest() %>%
      mutate(item = forcats::fct_inorder(item)) %>%
      spread(key = "item", value = "resp") %>%
      remove_rownames()
    
    data_check <- gather(response_data, item, response, -respondent) %>%
      select(item, response) %>%
      distinct() %>%
      group_by(item) %>%
      summarize(count = n()) %>%
      filter(count < 2)
    
    if (nrow(data_check) == 0) {
      good_data <- TRUE
    }
  }
  
  list(
    strc_params = strc$strc_params,
    strc_probs = strc$strc_probs,
    resp_profiles = strc$resp_profiles,
    q_matrix = meas$q_matrix,
    meas_params = meas$meas_params,
    response_data = response_data
  )
}
gen_mplus <- function(data, qmat, cores = parallel::detectCores() - 1) {
  mplus_syntax <- list()
  tracker <- 1
  
  mplus_syntax[[tracker]] <- paste("TITLE:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("Dissertation - Simulation")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("DATA:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("FILE IS simdata.csv;")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("VARIABLE:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("NAMES ARE")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("ID ", paste("X",
    seq_len(nrow(qmat)), sep = "", collapse = " "), ";", sep = "")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("USEVARIABLE ARE")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste(paste("X", seq_len(nrow(qmat)),
    sep = "", collapse = " "), ";", sep = "")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("CATEGORICAL ARE")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste(paste("X", seq_len(nrow(qmat)),
    sep = "", collapse = " "), ";", sep = "")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("CLASSES = c(", 2^ncol(qmat), ");", sep = "")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("IDvariable=ID;")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("MISSING=.;")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("ANALYSIS:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("TYPE=MIXTURE;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("STARTS=0;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("PROCESSORS=", cores, ";", sep = "")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("MITERATIONS=10000;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("LOGCRITERION=0.001;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("MCONVERGENCE=0.001;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("RLOGCRITERION=0.001;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("MCCONVERGENCE=0.001;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("MUCONVERGENCE=0.001;")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("MODEL:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("%OVERALL%")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  # Create item to profile table
  prof_list <- data_frame(att1 = c(0, 1))
  for (i in seq_along(qmat)) {
    prof_list[,paste0("att", i)] <- c(0, 1)
  }
  att_profiles <- expand.grid(prof_list) %>%
    mutate(row_sum = rowSums(.)) %>%
    group_by(row_sum) %>%
    nest(.key = profiles) %>%
    mutate(new_profiles = map(profiles,
      .f = function(x) arrange_all(x, desc))) %>%
    select(-profiles) %>%
    unnest() %>%
    select(-row_sum)
  itemprofile <- as_data_frame(
    matrix(data = NA, nrow = nrow(qmat), ncol = nrow(att_profiles))
  )
  for (i in seq_len(nrow(itemprofile))) {
    cur_q <- qmat[i,] %>% as.numeric()
    intercept <- paste0("L", i, "_0")
    for (a in seq_along(itemprofile)) {
      cur_profile <- att_profiles[a,] %>% flatten_dbl()
      
      prof_effects <- list()
      
      # intercept
      prof_effects[[1]] <- intercept
      
      # main effects
      for (attrib in seq_along(cur_profile)) {
        if (cur_profile[attrib] == 1 & cur_q[attrib] == 1) {
          prof_effects[[length(prof_effects) + 1]] <- paste0("L", i, "_1",
            attrib)
        }
      }
      
      # interactions
      interact <- which(cur_profile == 1 & cur_q == 1)
      if (length(interact) > 1) {
        for (inter in 2:length(interact)) {
          combos <- combn(x = interact, m = inter)
          for (c in 1:ncol(combos)) {
            cur_atts <- paste(combos[,c], collapse = "")
            prof_effects[[length(prof_effects) + 1]] <- paste0("L", i, "_",
              inter, cur_atts)
          }
        }
      }
      
      prof_effects <- flatten_chr(prof_effects) %>% paste(collapse = "+")
      itemprofile[i,a] <- prof_effects
    }
  }
  
  # Create threshold labels
  label_lookup <- list()
  for (i in seq_len(nrow(itemprofile))) {
    effects <- itemprofile[i,] %>% flatten_chr() %>% unique()
    
    item_effects <- data_frame(
      item = i,
      effect = effects,
      label = paste0("t", i, "_", seq_along(effects))
    )
    label_lookup[[i]] <- item_effects
  }
  label_lookup <- bind_rows(label_lookup)
  
  label_table <- as_data_frame(
    matrix(data = NA, nrow = nrow(itemprofile), ncol = ncol(itemprofile))
  )
  for (i in seq_len(nrow(label_table))) {
    for (j in seq_along(label_table)) {
      label_table[i,j] <- label_lookup$label[which(label_lookup$effect ==
          flatten_chr(itemprofile[i,j]))]
    }
  }
  
  # Define latent variable means
  for (i in seq_len(nrow(att_profiles) - 1)) {
    mplus_syntax[[tracker]] <- paste("[c#", i, "] (m", i, ");",
      " ! Latent variable mean for class ", i, sep = "")
    tracker <- tracker + 1
  }
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  # Define item labels for each class
  for (a in seq_along(label_table)) {
    cur_profile <- att_profiles[a,] %>% flatten_dbl() %>% paste(collapse = "")
    
    mplus_syntax[[tracker]] <- paste0("%c#", a, "%  ! for latent class ", a,
      "=", paste(cur_profile, collapse = ""))
    tracker <- tracker + 1
    
    for (i in seq_len(nrow(label_table))) {
      mplus_syntax[[tracker]] <- paste0("[X", i, "$1]  (",
        flatten_chr(label_table[i,a]), ");")
      tracker <- tracker + 1
    }
  }
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("MODEL CONSTRAINT:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("! STRUCTURAL MODEL PORTION")
  tracker <- tracker + 1
  
  strc_effects <- "g_0"
  for (a in seq_along(qmat)) {
    cur_effects <- combn(x = seq_along(qmat), m = a) %>%
      t() %>%
      as_data_frame() %>%
      apply(1, paste, collapse = "")
    
    strc_effects <- c(strc_effects, paste0("g_", a, cur_effects))
  }
  
  mplus_syntax[[tracker]] <- paste0("NEW(", paste(strc_effects, collapse = " "),
    ");")
  tracker <- tracker + 1
  
  intercept <- strc_effects[-1]
  
  for (a in seq_len(nrow(att_profiles) - 1)) {
    cur_profile <- att_profiles[a,] %>% flatten_dbl()
    cur_effects <- paste0("g_", 1, which(cur_profile == 1))
    if (sum(cur_profile) >= 2) {
      inter2 <- paste0("g_", 2, combn(x = which(cur_profile == 1), m = 2) %>%
          t() %>% as.data.frame() %>% mutate(name = paste0(V1, V2)) %>%
          pull(name))
      cur_effects <- c(cur_effects, inter2)
    }
    if (sum(cur_profile) >= 3) {
      inter3 <- paste0("g_", 3, combn(x = which(cur_profile == 1), m = 3) %>%
          t() %>% as.data.frame() %>% mutate(name = paste0(V1, V2, V3)) %>%
          pull(name))
      cur_effects <- c(cur_effects, inter3)
    }
    
    cur_effects <- cur_effects[which(nchar(cur_effects) > 3)]
    mplus_syntax[[tracker]] <- paste0("m", a, "=", paste(cur_effects,
      collapse = "+"), "-(", paste(intercept, collapse = "+"), ");")
    tracker <- tracker + 1
  }
  
  mplus_syntax[[tracker]] <- paste0("g_0=-(", paste(intercept, collapse = "+"),
    ");")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  # Model constraints syntax
  for (i in seq_len(nrow(itemprofile))) {
    cur_effects <- itemprofile[i,] %>%
      flatten_chr() %>%
      as.list() %>%
      map(.x = ., .f = function(x) unlist(strsplit(x, "+", fixed = TRUE))) %>%
      flatten_chr() %>%
      unique() %>%
      sort()
    
    cur_q <- qmat[i,] %>% as.numeric() %>% paste(collapse = "")
    
    mplus_syntax[[tracker]] <- paste0("! ITEM ", i, ":")
    tracker <- tracker + 1
    mplus_syntax[[tracker]] <- paste0("! Q-matrix Entry ", cur_q)
    tracker <- tracker + 1
    mplus_syntax[[tracker]] <- paste0("NEW (", paste(cur_effects,
      collapse = " "), ");")
    tracker <- tracker + 1
    
    label <- label_lookup %>%
      filter(item == i)
    
    for (l in 1:nrow(label)) {
      mplus_syntax[[tracker]] <- paste0(label$label[l], "=-(", label$effect[l],
        ");")
      tracker <- tracker + 1
    }
    
    mplus_syntax[[tracker]] <- paste0("! Order constraints")
    tracker <- tracker + 1
    
    maineffects <- grep("_1", cur_effects, value = TRUE)
    for (m in seq_along(maineffects)) {
      mplus_syntax[[tracker]] <- paste0(maineffects[m], ">0;")
      tracker <- tracker + 1
    }
    
    interaction_levels <- str_split(cur_effects, "_") %>%
      map(., function(x) {
        x %>% .[2] %>% str_sub(start = 1, end = 1) %>% as.numeric()
      }) %>%
      flatten_dbl() %>%
      .[which(. >= 2)] %>%
      unique()
    
    if (length(interaction_levels) == 0) {
      mplus_syntax[[tracker]] <- ""
      tracker <- tracker + 1
      next
    }
    
    effect_info <- data_frame(
      param = cur_effects
    ) %>%
      mutate(
        level = map_dbl(param, function(x) {
          str_split(x, "_") %>%
            flatten_chr() %>%
            .[2] %>%
            str_sub(start = 1, end = 1) %>%
            as.numeric()
        }),
        attributes = map_chr(param, function(x) {
          str_split(x, "_") %>%
            flatten_chr() %>%
            .[2] %>%
            str_sub(start = 2, end = -1)
        })
      )
    for (inter_level in min(interaction_levels):max(interaction_levels)) {
      interactions <- str_subset(cur_effects, paste0("_", inter_level))
      for (inter in seq_along(interactions)) {
        cur_inter <- interactions[inter]
        involved_atts <- str_split(cur_inter, "_") %>%
          flatten_chr() %>%
          .[2] %>%
          str_sub(start = 2, end = -1) %>%
          str_split("") %>%
          flatten_chr()
        
        for (att in seq_along(involved_atts)) {
          params <- effect_info %>%
            filter(
              level < inter_level,
              str_detect(attributes, involved_atts[att])
            ) %>%
            pull(param)
          mplus_syntax[[tracker]] <- paste0(cur_inter, ">-(",
            paste(params, collapse = "+"), ");")
          tracker <- tracker + 1
        }
      }
    }
    
    mplus_syntax[[tracker]] <- ""
    tracker <- tracker + 1
  }
  
  mplus_syntax[[tracker]] <- paste("OUTPUT:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("TECH1 TECH5 TECH8 TECH10;")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("SAVEDATA:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("FORMAT is f10.5;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("FILE is sim.dat;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("SAVE = CPROBABILITIES;")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("RESULTS ARE sim.res;")
  tracker <- tracker + 1
  
  return(mplus_syntax)
}
eval_mplus <- function(safe_results, path, true_data) {
  model_summary <- quiet_read(target = path) %>%
    .[["result"]] %>%
    .[["summaries"]] %>%
    as_data_frame()
  
  outdir <- str_split(path, "/") %>%
    flatten_chr() %>%
    .[-length(.)] %>%
    paste(collapse = "/")
  model_name <- str_split(path, "/") %>%
    flatten_chr() %>%
    .[length(.)]
  if (!is.null(safe_results[[2]]) | ncol(model_summary) == 5 | 
      !file.exists(paste0(outdir, "/sim.dat"))) {
    return(NULL)
  } else {
    outfiletext <- scan(path, what = "character",
      sep="\n", strip.white = FALSE, blank.lines.skip = FALSE, quiet = TRUE)
    Tech10 <- getSection("^\\s*TECHNICAL 10 OUTPUT\\s*$", outfiletext)
    uni_chi <- extractValue(
      pattern = "^\\s*Overall Univariate Pearson Chi-Square\\s*",
      Tech10, model_name, type = "dec")
    uni_log_chi <- extractValue(
      pattern="^\\s*Overall Univariate Log-Likelihood Chi-Square\\s*",
      Tech10, model_name, type = "dec")
    bi_chi <- extractValue(
      pattern = "^\\s*Overall Bivariate Pearson Chi-Square\\s*",
      Tech10, model_name, type = "dec")
    bi_log_chi <- extractValue(
      pattern="^\\s*Overall Bivariate Log-Likelihood Chi-Square\\s*",
      Tech10, model_name, type = "dec")
    uni_chi_df <- extractValue(
      pattern = "^\\s*Univariate Pearson Chi-Square\\s*",
      Tech10, model_name, type = "dec")
    bi_chi_df <- extractValue(
      pattern = "^\\s*Bivariate Pearson Chi-Square\\s*",
      Tech10, model_name, type = "dec")
    
    modfit <- getSection("^\\s*MODEL FIT INFORMATION\\s*$", outfiletext)
    aic <- extractValue(pattern = "^\\s*Akaike\\s*", modfit, model_name,
      type = "dec")
    bic <- extractValue(pattern = "^\\s*Bayesian\\s*", modfit, model_name,
      type = "dec")
    abic <- extractValue(pattern = "^\\s*Sample-Size Adjusted BIC\\s*", modfit,
      model_name, type = "dec")
    log_lik <- extractValue(pattern = "^\\s*H0 Value\\s*", modfit,
      model_name, type = "dec")
    
    quality <- getSection("^\\s*QUALITY OF NUMERICAL RESULTS\\s*$", outfiletext)
    cond_number <- extractValue(
      pattern = "^\\s*Condition Number for the Information Matrix\\s*", quality,
      model_name, type = "str") %>%
      as.numeric()
    
    model_params <- quiet_read(target = path) %>%
      .[["result"]] %>%
      .[["parameters"]] %>%
      .[["unstandardized"]] %>%
      as_data_frame() %>%
      filter(paramHeader == "New.Additional.Parameters") %>%
      select(param, est, se, est_se, pval) %>%
      mutate(
        item = str_sub(param, 2, regexpr("_", param) - 1) %>% as.numeric(),
        parameter = recode(
          str_sub(param, regexpr("_", param) + 1, regexpr("_", param) + 1),
          `0` = "intercept", `1` = "main_effect", `2` = "twoway_interaction",
          `3` = "threeway_interaction", `4` = "fourway_interaction"
        ),
        param = case_when(
          parameter == "fourway_interaction" & !str_detect(param, "4$") ~ paste0(param, "4"),
          TRUE ~ param
        )
      ) %>%
      select(param, item, parameter, everything()) %>%
      mutate_at(vars(est:pval), as.numeric)
    
    model_profiles <- read_table(file = file.path(outdir, "sim.dat"), na = "*",
      col_names = c(
        paste0("item", seq_len(30)),
        "id",
        paste0("profprob_", seq_len(2^(ncol(true_data$resp_profiles) - 1))),
        "profile"
      ), col_types = cols(.default = col_double())) %>%
      select(-(seq_len(30)), -profile, -id) %>%
      mutate_if(is.character, as.numeric) %>%
      remove_rownames()
    for (i in seq_len(ncol(true_data$resp_profiles) - 1)) {
      att_num <- rlang::sym(paste0("att", i))
      cur_profiles <- true_data$strc_probs %>%
        rowid_to_column(var = "class") %>%
        select(-(sum:prob)) %>%
        filter((!!att_num) == 1) %>%
        pull(class)
      class_names <- rlang::syms(as.list(paste0("profprob_", cur_profiles)))
      class_probs <- model_profiles %>%
        select(!!!class_names) %>%
        rowSums()
      model_profiles[, paste0("att", i)] <- class_probs
    }
    
    est_profile <- model_profiles %>%
      select(grep("att", colnames(.))) %>%
      t() %>%
      as_data_frame() %>%
      as.list() %>%
      map(.x = ., .f = master_round, cut = 0.5)
    
    # model performance
    parameter_recovery <- bind_rows(true_data$strc_params,
      true_data$meas_params) %>%
      full_join(select(model_params, param, est, se, est_se, pval),
        by = c("parameter" = "param")) %>%
      mutate(
        mismatch = case_when(
          parameter == "G_0" ~ FALSE,
          is.na(value) ~ TRUE,
          is.na(est) ~ TRUE,
          TRUE ~ FALSE
        ),
        type = map_chr(parameter, function(x) {
          initial <- str_split(x, "_") %>% flatten_chr()
          if (initial[1] == "G") return("structural")
          level <- initial[2] %>%
            str_sub(start = 1, end = 1)
          if (level == 0) {
            return("intercept")
          } else if (level == 1) {
            return("main_effect")
          } else {
            return(paste0("interaction", level))
          }
        }),
        item = map_dbl(parameter, function(x) {
          initial <- str_split(x, "_") %>% flatten_chr()
          if (initial[1] == "G") return(NA_real_)
          initial[1] %>% str_replace("L", "") %>% as.numeric()
        })
      ) %>%
      replace_na(list(value = 0, est = 0)) %>%
      mutate(value = case_when(parameter == "G_0" ~ NA_real_, TRUE ~ value)) %>%
      arrange(item, parameter) %>%
      mutate(
        bias = est - value,
        sq_error = bias^2
      )
    
    real_profile <- true_data$resp_profiles %>%
      select(-respondent) %>%
      t() %>%
      as_data_frame() %>%
      as.list()
    all_profs <- true_data$strc_probs %>%
      select(-(sum:prob)) %>%
      t() %>%
      as_data_frame() %>%
      as.list()
    est_master <- map_dbl(.x = est_profile, .f = function(x, all_profs) {
      which(map_lgl(all_profs, identical, y = x)) %>% as.numeric()
    },
      all_profs = all_profs)
    real_master <- map_dbl(.x = real_profile, .f = function(x, all_profs) {
      which(map_lgl(all_profs, identical, y = x)) %>% as.numeric()
    },
      all_profs = all_profs)
    kappa <- calc_kappa(x = real_master, y = est_master, min_score = 1,
      max_score = length(all_profs))
    ccr <- length(which(real_master == est_master)) / length(real_master)
    
    kappa_values <- map2_dbl(.x = real_profile, .y = est_profile,
      .f = calc_kappa, min_score = 0, max_score = 1)
    ccr_values <- map2_dbl(.x = real_profile, .y = est_profile,
      .f = function(x, y) length(which(x == y)) / length(x))
    
    # create results
    results <- data_frame(
      summary = list(model_summary),
      sim_data = list(true_data),
      profiles_probs = list(model_probs = model_profiles),
      est_profiles = list(model_profs = est_profile),
      param_recovery = list(param_recovery = parameter_recovery),
      correct_strc = ifelse(any(parameter_recovery %>%
          filter(type == "structural") %>% pull(mismatch)), FALSE, TRUE),
      correct_meas = ifelse(any(parameter_recovery %>%
          filter(type != "structural") %>% pull(mismatch)), FALSE, TRUE),
      attribute_kappa = calc_mean(kappa_values),
      attribute_ccr = calc_mean(ccr_values),
      profile_kappa = kappa,
      profile_ccr = ccr,
      num_est_profile = length(unique(est_profile)),
      mplus_condition_number = cond_number,
      log_lik,
      aic,
      bic,
      abic,
      uni_chi,
      uni_log_chi,
      uni_chi_df,
      bi_chi,
      bi_log_chi,
      bi_chi_df
    )
    
    recover_sum <- parameter_recovery %>%
      group_by(type) %>%
      summarize(
        bias = mean(bias, na.rm = TRUE),
        mse = mean(sq_error, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      gather(error, value, bias:mse) %>%
      arrange(type) %>%
      mutate(full_type = paste0(type, "_", error)) %>%
      select(full_type, value) %>%
      spread(key = full_type, value = value)
    
    results <- bind_cols(results, recover_sum)
    return(results)
  }
}
est_mplus <- function(data, qmat, cores = parallel::detectCores() - 1, cond,
  rep, ospec) {
  cur_path <- file.path("analyses", "replications", paste0("cond", cond),
    paste0("rep_", sprintf("%03d", rep)), paste0("over_",
      sprintf("%02d", ospec)))
  
  # saturated model ------------------------------------------------------------
  satr_mplus <- gen_mplus(data = data$response_data, qmat = select(qmat, -item),
    cores = cores)
  print_mplus <- map(satr_mplus, check_length)
  
  inp_path <- file.path(cur_path, "satr", "satr_model.inp")
  inp_file <- file(inp_path)
  writeLines(paste(print_mplus), con = inp_file, sep = "\n")
  close(inp_file)
  
  mplus_output <- safe_mplus(file.path(getwd(), cur_path, "satr",
    "satr_model.inp"))

  satr_model <- eval_mplus(safe_results = mplus_output,
    path = file.path(getwd(), cur_path, "satr", "satr_model.out"),
    true_data = data)
  
  if (is.null(satr_model)) {
    satr_model <- data_frame(
      rep = rep,
      over_spec = ospec / 100,
      model = "satr",
      rm_param = NA,
      reduction = NA,
      converge = FALSE,
      satr_converge = FALSE
    )
    satr_nonsig <- quiet_read(target = file.path(getwd(), cur_path, "satr",
      "satr_model.out")) %>%
      .[["result"]] %>%
      .[["parameters"]] %>%
      .[["unstandardized"]] %>%
      as_data_frame() %>%
      filter(paramHeader == "New.Additional.Parameters") %>%
      select(param, est) %>%
      mutate(
        item = str_sub(param, 2, regexpr("_", param) - 1) %>% as.numeric(),
        parameter = recode(
          str_sub(param, regexpr("_", param) + 1, regexpr("_", param) + 1),
          `0` = "intercept", `1` = "main_effect", `2` = "twoway_interaction",
          `3` = "threeway_interaction", `4` = "fourway_interaction"
        ),
        param = case_when(
          parameter == "fourway_interaction" & !str_detect(param, "4$") ~ paste0(param, "4"),
          TRUE ~ param
        )
      ) %>%
      select(param, type = parameter) %>%
      filter(str_detect(type, "threeway|fourway")) %>%
      select(parameter = param)
  } else {
    satr_model = bind_cols(
      data_frame(rep = rep, over_spec = ospec / 100, model = "satr",
        rm_param = NA, reduction = NA, converge = TRUE, satr_converge = TRUE),
      satr_model
    )
    satr_nonsig <- satr_model$param_recovery$param_recovery %>%
      filter(type != "intercept", (pval > 0.05 | se == 0))
  }
  
  # simultaneous reduction------------------------------------------------------
  strc_remove <- satr_nonsig %>%
    filter(str_detect(parameter, "G")) %>%
    pull(parameter) %>%
    tolower()
  meas_remove <- satr_nonsig %>%
    filter(str_detect(parameter, "L")) %>%
    pull(parameter)
  remove <- c(strc_remove, meas_remove)
  
  siml_mplus <- satr_mplus
  for (i in seq_along(remove)) {
    siml_mplus <- gsub(paste0("+", remove[i]), "", siml_mplus, fixed = TRUE)
    siml_mplus <- gsub(paste0(remove[i], "+"), "", siml_mplus, fixed = TRUE)
    siml_mplus <- gsub(paste0(" ", remove[i]), "", siml_mplus, fixed = TRUE)
    siml_mplus <- gsub(paste0(remove[i], " "), "", siml_mplus, fixed = TRUE)
    
    del_rows <- grep(remove[i], siml_mplus)
    if (length(del_rows) > 0) siml_mplus <- siml_mplus[-del_rows]
  }
  siml_mplus <- as.list(siml_mplus)
  print_mplus <- map(siml_mplus, check_length)
  
  inp_path <- file.path(cur_path, "siml", "siml_model.inp")
  inp_file <- file(inp_path)
  writeLines(paste(print_mplus), con = inp_file, sep = "\n")
  close(inp_file)
  
  mplus_output <- safe_mplus(file.path(getwd(), cur_path, "siml",
    "siml_model.inp"))
  
  siml_model <- eval_mplus(safe_results = mplus_output,
    path = file.path(getwd(), cur_path, "siml", "siml_model.out"),
    true_data = data)
  
  if (is.null(siml_model)) {
    siml_model <- data_frame(
      rep = rep,
      over_spec = ospec / 100,
      model = "siml",
      rm_param = list(filter(satr_nonsig, parameter %in% toupper(remove))),
      reduction = ifelse(satr_model$converge, "pval", "rule"),
      converge = FALSE,
      satr_converge = satr_model$satr_converge
    )
  } else {
    siml_model = bind_cols(
      data_frame(rep = rep, over_spec = ospec / 100, model = "siml",
        rm_param = list(filter(satr_nonsig, parameter %in% toupper(remove))),
        reduction = ifelse(satr_model$converge, "pval", "rule"),
        converge = TRUE, satr_converge = satr_model$satr_converge),
      siml_model
    )
  }
  
  # measurement reduction ------------------------------------------------------
  remove <- meas_remove
  
  meas_mplus <- satr_mplus
  for (i in seq_along(remove)) {
    meas_mplus <- gsub(paste0("+", remove[i]), "", meas_mplus, fixed = TRUE)
    meas_mplus <- gsub(paste0(remove[i], "+"), "", meas_mplus, fixed = TRUE)
    meas_mplus <- gsub(paste0(" ", remove[i]), "", meas_mplus, fixed = TRUE)
    meas_mplus <- gsub(paste0(remove[i], " "), "", meas_mplus, fixed = TRUE)
    
    del_rows <- grep(remove[i], meas_mplus)
    if (length(del_rows) > 0) meas_mplus <- meas_mplus[-del_rows]
  }
  meas_mplus <- as.list(meas_mplus)
  print_mplus <- map(meas_mplus, check_length)
  
  inp_path <- file.path(cur_path, "meas", "meas_model.inp")
  inp_file <- file(inp_path)
  writeLines(paste(print_mplus), con = inp_file, sep = "\n")
  close(inp_file)
  
  mplus_output <- safe_mplus(file.path(getwd(), cur_path, "meas",
    "meas_model.inp"))
  
  meas_model <- eval_mplus(safe_results = mplus_output,
    path = file.path(getwd(), cur_path, "meas", "meas_model.out"),
    true_data = data)
  
  if (is.null(meas_model)) {
    meas_model <- data_frame(
      rep = rep,
      over_spec = ospec / 100,
      model = "meas",
      rm_param = list(filter(satr_nonsig, parameter %in% toupper(remove))),
      reduction = ifelse(satr_model$converge, "pval", "rule"),
      converge = FALSE,
      satr_converge = satr_model$satr_converge
    )
  } else {
    meas_model = bind_cols(
      data_frame(rep = rep, over_spec = ospec / 100, model = "meas",
        rm_param = list(filter(satr_nonsig, parameter %in% toupper(remove))),
        reduction = ifelse(satr_model$converge, "pval", "rule"),
        converge = TRUE, satr_converge = satr_model$satr_converge),
      meas_model
    )
    meas_nonsig <- meas_model$param_recovery$param_recovery %>%
      filter(type != "intercept", (pval > 0.05 | se == 0))
  }
  
  # structural reduction -------------------------------------------------------
  remove <- strc_remove
  
  strc_mplus <- satr_mplus
  for (i in seq_along(remove)) {
    strc_mplus <- gsub(paste0("+", remove[i]), "", strc_mplus, fixed = TRUE)
    strc_mplus <- gsub(paste0(remove[i], "+"), "", strc_mplus, fixed = TRUE)
    strc_mplus <- gsub(paste0(" ", remove[i]), "", strc_mplus, fixed = TRUE)
    strc_mplus <- gsub(paste0(remove[i], " "), "", strc_mplus, fixed = TRUE)
    
    del_rows <- grep(remove[i], strc_mplus)
    if (length(del_rows) > 0) strc_mplus <- strc_mplus[-del_rows]
  }
  strc_mplus <- as.list(strc_mplus)
  print_mplus <- map(strc_mplus, check_length)
  
  inp_path <- file.path(cur_path, "strc", "strc_model.inp")
  inp_file <- file(inp_path)
  writeLines(paste(print_mplus), con = inp_file, sep = "\n")
  close(inp_file)
  
  mplus_output <- safe_mplus(file.path(getwd(), cur_path, "strc",
    "strc_model.inp"))
  
  strc_model <- eval_mplus(safe_results = mplus_output,
    path = file.path(getwd(), cur_path, "strc", "strc_model.out"),
    true_data = data)
  
  if (is.null(strc_model)) {
    strc_model <- data_frame(
      rep = rep,
      over_spec = ospec / 100,
      model = "strc",
      rm_param = list(filter(satr_nonsig, parameter %in% toupper(remove))),
      reduction = ifelse(satr_model$converge, "pval", "rule"),
      converge = FALSE,
      satr_converge = satr_model$satr_converge
    )
  } else {
    strc_model = bind_cols(
      data_frame(rep = rep, over_spec = ospec / 100, model = "strc",
        rm_param = list(filter(satr_nonsig, parameter %in% toupper(remove))),
        reduction = ifelse(satr_model$converge, "pval", "rule"),
        converge = TRUE, satr_converge = satr_model$satr_converge),
      strc_model
    )
    strc_nonsig <- strc_model$param_recovery$param_recovery %>%
      filter(type != "intercept", (pval > 0.05 | se == 0))
  }
  
  # measurement-structural reduction -------------------------------------------
  if (meas_model$converge) {
    remove <- meas_nonsig %>%
      filter(str_detect(parameter, "G")) %>%
      pull(parameter) %>%
      tolower()
    
    measstrc_mplus <- meas_mplus
    for (i in seq_along(remove)) {
      measstrc_mplus <- gsub(paste0("+", remove[i]), "", measstrc_mplus,
        fixed = TRUE)
      measstrc_mplus <- gsub(paste0(remove[i], "+"), "", measstrc_mplus,
        fixed = TRUE)
      measstrc_mplus <- gsub(paste0(" ", remove[i]), "", measstrc_mplus,
        fixed = TRUE)
      measstrc_mplus <- gsub(paste0(remove[i], " "), "", measstrc_mplus,
        fixed = TRUE)
      
      del_rows <- grep(remove[i], measstrc_mplus)
      if (length(del_rows) > 0) measstrc_mplus <- measstrc_mplus[-del_rows]
    }
    measstrc_mplus <- as.list(measstrc_mplus)
    print_mplus <- map(measstrc_mplus, check_length)
    
    inp_path <- file.path(cur_path, "meas-strc", "meas-strc_model.inp")
    inp_file <- file(inp_path)
    writeLines(paste(print_mplus), con = inp_file, sep = "\n")
    close(inp_file)
    
    mplus_output <- safe_mplus(file.path(getwd(), cur_path, "meas-strc",
      "meas-strc_model.inp"))
    
    measstrc_model <- eval_mplus(safe_results = mplus_output,
      path = file.path(getwd(), cur_path, "meas-strc", "meas-strc_model.out"),
      true_data = data)
    
    if (is.null(measstrc_model)) {
      measstrc_model <- data_frame(
        rep = rep,
        over_spec = ospec / 100,
        model = "meas-strc",
        rm_param = list(bind_rows(
          meas_model$rm_param[[1]],
          filter(meas_nonsig, parameter %in% toupper(remove))
        )),
        reduction = paste0(meas_model$reduction, ",",
          ifelse(meas_model$converge, "pval", "rule")),
        converge = FALSE,
        satr_converge = satr_model$satr_converge
      )
    } else {
      measstrc_model = bind_cols(
        data_frame(rep = rep, over_spec = ospec / 100, model = "meas-strc",
          rm_param = list(bind_rows(
            meas_model$rm_param[[1]],
            filter(meas_nonsig, parameter %in% toupper(remove))
          )),
          reduction = paste0(meas_model$reduction, ",",
            ifelse(meas_model$converge, "pval", "rule")),
          converge = TRUE, satr_converge = satr_model$satr_converge
        ),
        measstrc_model
      )
    }
  } else {
    measstrc_model <- data_frame(
      rep = rep,
      over_spec = ospec / 100,
      model = "meas-strc",
      rm_param = NA,
      reduction = NA,
      converge = NA,
      satr_converge = satr_model$satr_converge
    )
  }
  
  # structural-measurement reduction -------------------------------------------
  if (strc_model$converge) {
    remove <- strc_nonsig %>%
      filter(str_detect(parameter, "L")) %>%
      pull(parameter)
    
    strcmeas_mplus <- strc_mplus
    for (i in seq_along(remove)) {
      strcmeas_mplus <- gsub(paste0("+", remove[i]), "", strcmeas_mplus,
        fixed = TRUE)
      strcmeas_mplus <- gsub(paste0(remove[i], "+"), "", strcmeas_mplus,
        fixed = TRUE)
      strcmeas_mplus <- gsub(paste0(" ", remove[i]), "", strcmeas_mplus,
        fixed = TRUE)
      strcmeas_mplus <- gsub(paste0(remove[i], " "), "", strcmeas_mplus,
        fixed = TRUE)
      
      del_rows <- grep(remove[i], strcmeas_mplus)
      if (length(del_rows) > 0) strcmeas_mplus <- strcmeas_mplus[-del_rows]
    }
    strcmeas_mplus <- as.list(strcmeas_mplus)
    print_mplus <- map(strcmeas_mplus, check_length)
    
    inp_path <- file.path(cur_path, "strc-meas", "strc-meas_model.inp")
    inp_file <- file(inp_path)
    writeLines(paste(print_mplus), con = inp_file, sep = "\n")
    close(inp_file)
    
    mplus_output <- safe_mplus(file.path(getwd(), cur_path, "strc-meas",
      "strc-meas_model.inp"))
    
    strcmeas_model <- eval_mplus(safe_results = mplus_output,
      path = file.path(getwd(), cur_path, "strc-meas", "strc-meas_model.out"),
      true_data = data)
    
    if (is.null(strcmeas_model)) {
      strcmeas_model <- data_frame(
        rep = rep,
        over_spec = ospec / 100,
        model = "strc-meas",
        rm_param = list(bind_rows(
          strc_model$rm_param[[1]],
          filter(strc_nonsig, parameter %in% toupper(remove))
        )),
        reduction = paste0(strc_model$reduction, ",",
          ifelse(strc_model$converge, "pval", "rule")),
        converge = FALSE,
        satr_converge = satr_model$satr_converge
      )
    } else {
      strcmeas_model = bind_cols(
        data_frame(rep = rep, over_spec = ospec / 100, model = "strc-meas",
          rm_param = list(bind_rows(
            strc_model$rm_param[[1]],
            filter(strc_nonsig, parameter %in% toupper(remove))
          )),
          reduction = paste0(strc_model$reduction, ",",
            ifelse(strc_model$converge, "pval", "rule")),
          converge = TRUE, satr_converge = satr_model$satr_converge
        ),
        strcmeas_model
      )
    }
  } else {
    strcmeas_model <- data_frame(
      rep = rep,
      over_spec = ospec / 100,
      model = "strc-meas",
      rm_param = NA,
      reduction = NA,
      converge = NA,
      satr_converge = satr_model$satr_converge
    )
  }
  
  # combine reduction models ---------------------------------------------------
  all_models <- bind_rows(satr_model, siml_model, meas_model, strc_model,
    measstrc_model, strcmeas_model) %>%
    rename(dataset = rep)
  
  write_rds(all_models, path = paste0(cur_path, "/all_models.rds"),
    compress = "gz")
  return(all_models)
}
single_sim <- function(cond, rep, attributes, sample_size, seeds) {
  # create output directory
  outdir <- file.path("analyses", "replications", paste0("cond", cond),
    paste0("rep_", sprintf("%03d", rep)))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  outdir <- normalizePath(outdir)
  outdir <- paste0(outdir, ifelse(Sys.info()["sysname"] == "Windows", "\\",
    "/"))
  if (Sys.info()["sysname"] == "Windows") {
    outdir <- gsub("\\", "/", outdir, fixed = TRUE)
  }
  
  # create estimation folders
  new_path <- crossing(
    over_spec = paste0("over_", c("00", "10", "20")),
    reduction = c("meas", "meas-strc", "satr", "siml", "strc", "strc-meas")
  ) %>%
    mutate(path = paste0(outdir, file.path(over_spec, reduction), "/")) %>%
    pull(path)
  walk(new_path, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))
  
  # set seed
  setSeeds(seeds, run = rep)
  useStream(cond)
  
  # generate data
  if (file.exists(paste0(outdir, "all_data.rds"))) {
    all_data <- read_rds(paste0(outdir, "all_data.rds"))
  } else {
    all_data <- gen_data(num_resp = sample_size, num_att = attributes)
  }
  write_rds(all_data, path = paste0(outdir, "all_data.rds"), compress = "gz")
  walk(new_path, function(x, data) {
    write_csv(data, paste0(x, "simdata.csv"), na = ".", col_names = FALSE)
  }, data = all_data$response_data)

  # estimate mplus with true q-matrix
  qmat_true <- all_data$q_matrix
  if (file.exists(paste0(outdir, "over_00/all_models.rds"))) {
    results_true <- read_rds(paste0(outdir, "over_00/all_models.rds"))
  } else {
    results_true <- est_mplus(data = all_data, qmat = qmat_true,
      cores = parallel::detectCores() - 1, cond = cond, rep = rep, ospec = 0)
  }
  
  # estimate mplus with 0.1 over specified q-matrix
  qmat_1 <- all_data$q_matrix
  for (i in 1:nrow(qmat_1)) {
    for (j in 1:ncol(qmat_1)) {
      if (qmat_1[i,j] != 0) {
        next
      } else {
        rand <- runif(1, min = 0, max = 1)
        if (rand <= 0.1) qmat_1[i,j] <- 1
      }
    }
  }
  if (file.exists(paste0(outdir, "over_10/all_models.rds"))) {
    results_1 <- read_rds(paste0(outdir, "over_10/all_models.rds"))
  } else {
    results_1 <- est_mplus(data = all_data, qmat = qmat_1,
      cores = parallel::detectCores() - 1, cond = cond, rep = rep, ospec = 10)
  }
  
  # estimate mplus with 0.2 over specified q-matrix
  qmat_2 <- all_data$q_matrix
  for (i in 1:nrow(qmat_2)) {
    for (j in 1:ncol(qmat_2)) {
      if (qmat_2[i,j] != 0) {
        next
      } else {
        rand <- runif(1, min = 0, max = 1)
        if (rand <= 0.2) qmat_2[i,j] <- 1
      }
    }
  }
  if (file.exists(paste0(outdir, "over_20/all_models.rds"))) {
    results_2 <- read_rds(paste0(outdir, "over_20/all_models.rds"))
  } else {
    results_2 <- est_mplus(data = all_data, qmat = qmat_2,
      cores = parallel::detectCores() - 1, cond = cond, rep = rep, ospec = 20)
  }
  
  # combine models
  all_results <- bind_rows(results_true, results_1, results_2)
  write_rds(all_results, path = paste0(outdir, "all_models.rds"),
    compress = "gz")
  return(all_results)
}
single_cond <- function(cond, attributes, sample_size, reps, seeds) {
  if (length(reps) == 1) {
    data_list <- data_frame(
      cond = cond,
      rep = seq_len(reps),
      attributes = attributes,
      sample_size = sample_size
    ) %>%
      as.list()
    
    outlist <- file.path("analyses", "replications", paste0("cond", cond),
      paste0("rep_", sprintf("%03d", seq_len(reps))))
    walk(.x = outlist, .f = function(x) {
      if (!dir.exists(x)) dir.create(x, recursive = TRUE)
    })
  } else {
    data_list <- data_frame(
      cond = cond,
      rep = reps,
      attributes = attributes,
      sample_size = sample_size
    ) %>%
      as.list()
    
    outlist <- file.path("analyses", "replications", paste0("cond", cond),
      paste0("rep_", sprintf("%03d", reps)))
    walk(.x = outlist, .f = function(x) {
      if (!dir.exists(x)) dir.create(x, recursive = TRUE)
    })
  }
  
  cond_results <- pmap_df(.l = data_list, .f = single_sim, seeds = seeds)

  cond_results <- bind_cols(
    data_frame(cond, attributes, sample_size) %>%
      sample_n(size = nrow(cond_results), replace = TRUE),
    cond_results
  )
  
  cond_recovery <- cond_results %>%
    select(cond, attributes, sample_size, over_spec, model, satr_converge,
      converge, correct_strc, correct_meas, attribute_kappa, attribute_ccr,
      profile_kappa, profile_ccr, num_est_profile,
      which(str_detect(colnames(.), "bias")),
      which(str_detect(colnames(.), "mse"))) %>%
    group_by(cond, attributes, sample_size, over_spec, model, satr_converge) %>%
    mutate(count = n(), converge = sum(converge, na.rm = TRUE) / count) %>%
    mutate_at(vars(contains("correct")), mean, na.rm = TRUE) %>%
    mutate_at(vars(contains("kappa")), calc_mean) %>%
    mutate_at(vars(contains("ccr")), calc_mean) %>%
    mutate_at(vars(contains("num_est")), mean, na.rm = TRUE) %>%
    mutate_at(vars(contains("bias")), median, na.rm = TRUE) %>%
    mutate_at(vars(contains("mse")), median, na.rm = TRUE) %>%
    distinct() %>%
    select(cond, attributes, sample_size, over_spec, model, converge,
      everything()) %>%
    arrange(cond, attributes, sample_size, over_spec, model)
    
  cond_fit <- cond_results %>%
    select(cond, attributes, sample_size, over_spec, model, dataset,
      satr_converge, aic, bic, abic) %>%
    group_by(cond, attributes, sample_size, over_spec, satr_converge) %>%
    nest(.key = "raw_data") %>%
    mutate(fit_data = map(raw_data, measure_fit)) %>%
    select(-raw_data) %>%
    unnest()
  
  write_path <- file.path("analyses", "replications", paste0("cond", cond))
  write_rds(cond_results, path = paste0(write_path, "/all_results.rds"),
    compress = "gz")
  write_rds(cond_recovery, path = paste0(write_path, "/recovery.rds"),
    compress = "gz")
  write_rds(cond_fit, path = paste0(write_path, "/fit.rds"), compress = "gz")
  
  list(
    all_results = cond_results,
    recovery = cond_recovery,
    fit = cond_fit
  )
}
full_sim <- function(conditions, reps, seeds) {
  sim_results <- pmap(.l = as.list(conditions), .f = single_cond,
    reps = reps, seeds = seeds)
  
  write_path <- file.path("analyses", "replications")
  write_rds(sim_results, path = paste0(write_path, "/full_results.rds"),
    compress = "gz")
  
  list(
    all_results = map(.x = sim_results, .f = function(x) x$all_results),
    recovery = map_df(.x = sim_results, .f = function(x) x$recovery),
    fit = map_df(.x = sim_results, .f = function(x) x$fit)
  )
}


### Setup simulation -----------------------------------------------------------
# define conditions
conditions <- crossing(
    attributes = c(3, 4),
    sample_size = c(500, 1000, 5000)
  ) %>%
  rowid_to_column(var = "cond")
reps <- 100

# create seeds
# projSeeds <- seedCreator(nReps = reps, streamsPerRep = nrow(conditions),
#   seed = 9416, file = "analyses/random-seeds.rds")
projSeeds <- read_rds("analyses/random-seeds.rds")

# filter conditions
conditions <- conditions %>%
  filter(between(cond, 1, 6))

# run simulation
sim_results <- full_sim(conditions = conditions, reps = reps, seeds = projSeeds)
write_rds(sim_results, path = "analyses/sim-results.rds", compress = "gz")
write_rds(sim_results$recovery, path = "analyses/sim-recovery.rds",
  compress = "gz")
write_rds(sim_results$fit, path = "analyses/sim-fit.rds", compress = "gz")

allreps <- sim_results %>%
  pluck("all_results") %>%
  bind_rows() %>%
  select(cond, attributes, sample_size, dataset, over_spec, model, rm_param,
    reduction, converge, satr_converge, correct_strc, correct_meas,
    attribute_kappa, attribute_ccr, profile_kappa, profile_ccr, aic, bic, abic,
    uni_chi, uni_chi_df, bi_chi, bi_chi_df, contains("bias"), contains("mse"))
write_rds(allreps, path = "analyses/sim-allreps.rds", compress = "gz")
