### Setup R session ------------------------------------------------------------
options(
  repos = c("https://cran.rstudio.com/", "http://rweb.crmda.ku.edu/kran"),
  scipen = 999
)

needed_packages <- c("tidyverse", "MplusAutomation")
load_packages <- function(x) {
  if (!(x %in% rownames(installed.packages()))) install.packages(x)
  
  suppressPackageStartupMessages(require(x, character.only = TRUE))
}
sapply(needed_packages, load_packages)


### Read in data ---------------------------------------------------------------
dtmr <- read_csv("data/dtmr-data.csv",
  col_names = c("id", paste0("X", seq_len(28)), "MCred", "ExpGd", "GdCert79",
    "GdCert78"), col_types = cols(.default = col_integer()), na = ".")

qmatrix <- read_csv("data/dtmr-qmatrix.csv",
  col_types = cols(.default = col_number(), Item = col_character()))


### Count missing data ---------------------------------------------------------
missing_data <- dtmr %>%
  select(id:X28) %>%
  gather(item, response, X1:X28)

pct_missing <- length(which(is.na(missing_data$response))) / nrow(missing_data)


### Define functions -----------------------------------------------------------
gen_nonfung_mplus <- function(data, qmat,
  cores = parallel::detectCores() - 1) {
  mplus_syntax <- list()
  tracker <- 1
  
  mplus_syntax[[tracker]] <- paste("TITLE:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("Dissertation - DTMR")
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- ""
  tracker <- tracker + 1
  
  mplus_syntax[[tracker]] <- paste("DATA:")
  tracker <- tracker + 1
  mplus_syntax[[tracker]] <- paste("FILE IS dtmr_mplus.csv;")
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
    mutate(new_profiles = map(profiles, .f = function(x) arrange_all(x, desc))) %>%
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
      prof_effects[[1]] <- intercept
      
      for (attrib in seq_along(cur_profile)) {
        if (cur_profile[attrib] == 1 & cur_q[attrib] == 1) {
          prof_effects[[length(prof_effects) + 1]] <- paste0("L", i, "_1",
            attrib)
        }
      }
      
      interact <- which(cur_profile == 1 & cur_q == 1)
      if (length(interact) > 1) {
        prof_effects[[length(prof_effects) + 1]] <- paste0("L", i, "_2",
          paste(interact, collapse = ""))
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
    
    interactions <- grep("_2", cur_effects, value = TRUE)
    for (inter in seq_along(interactions)) {
      for (m in seq_along(maineffects)) {
        mplus_syntax[[tracker]] <- paste0(interactions[inter], ">-",
          maineffects[m], ";")
        tracker <- tracker + 1
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


### Write data files -----------------------------------------------------------
model_paths <- c("satr/", "siml/", "meas/", "strc/", "meas-strc/", "strc-meas/")
model_paths <- paste0("analyses/dtmr-results/", model_paths)
walk(model_paths, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
  write_csv(select(dtmr, id:X28), paste0(x, "dtmr_mplus.csv"), na = ".",
    col_names = FALSE)
})


### Saturated model ------------------------------------------------------------
satr_mplus <- gen_nonfung_mplus(data = select(dtmr, id:X28),
  qmat = select(qmatrix, -Item), cores = 3)
print_mplus <- map(satr_mplus, check_length)

inp_file <- file("analyses/dtmr-results/satr/satr_model.inp")
writeLines(paste(print_mplus), con = inp_file, sep = "\n")
close(inp_file)

runModels(paste0(getwd(), "/analyses/dtmr-results/satr"))

satr_params <- extractModelParameters(target = paste0(getwd(),
  "/analyses/dtmr-results/satr"))[[1]] %>%
  as_data_frame() %>%
  filter(paramHeader == "New.Additional.Parameters")

write_csv(satr_params, "analyses/dtmr-results/satr/satr_params.csv")

rm(list = setdiff(ls(), c("satr_mplus", "satr_params", "check_length"))); gc()


### Simultaneous reduction -----------------------------------------------------
strc_remove <- satr_params %>%
  filter(grepl("G", param)) %>%
  filter(pval > 0.05) %>%
  pull(param) %>%
  tolower()
meas_remove <- satr_params %>%
  filter(grepl("L", param), !grepl("_0", param)) %>%
  filter((se == 0 | pval > 0.05) & !grepl("_0", param)) %>%
  pull(param)
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

inp_file <- file("analyses/dtmr-results/siml/siml_model.inp")
writeLines(paste(print_mplus), con = inp_file, sep = "\n")
close(inp_file)

runModels(paste0(getwd(), "/analyses/dtmr-results/siml"))

siml_params <- extractModelParameters(target = paste0(getwd(),
  "/analyses/dtmr-results/siml"))[[1]] %>%
  as_data_frame() %>%
  filter(paramHeader == "New.Additional.Parameters")

write_csv(siml_params, "analyses/dtmr-results/siml/siml_params.csv")

rm(list = setdiff(ls(), c("satr_mplus", "check_length",
  "meas_remove", "strc_remove"))); gc()


### Measurement reduction ------------------------------------------------------
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

inp_file <- file("analyses/dtmr-results/meas/meas_model.inp")
writeLines(paste(print_mplus), con = inp_file, sep = "\n")
close(inp_file)

runModels(paste0(getwd(), "/analyses/dtmr-results/meas"))

meas_params <- extractModelParameters(target = paste0(getwd(),
  "/analyses/dtmr-results/meas"))[[1]] %>%
  as_data_frame() %>%
  filter(paramHeader == "New.Additional.Parameters")

write_csv(meas_params, "analyses/dtmr-results/meas/meas_params.csv")

rm(list = setdiff(ls(), c("satr_mplus", "meas_mplus", "check_length",
  "strc_remove", "meas_params"))); gc()


### Structural reduction -------------------------------------------------------
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

inp_file <- file("analyses/dtmr-results/strc/strc_model.inp")
writeLines(paste(print_mplus), con = inp_file, sep = "\n")
close(inp_file)

runModels(paste0(getwd(), "/analyses/dtmr-results/strc"))

strc_params <- extractModelParameters(target = paste0(getwd(),
  "/analyses/dtmr-results/strc"))[[1]] %>%
  as_data_frame() %>%
  filter(paramHeader == "New.Additional.Parameters")

write_csv(strc_params, "analyses/dtmr-results/strc/strc_params.csv")

rm(list = setdiff(ls(), c("meas_mplus", "strc_mplus", "check_length",
  "meas_params", "strc_params"))); gc()


### Measurement-structural reduction -------------------------------------------
remove <- meas_params %>%
  filter(grepl("G", param)) %>%
  filter(pval > 0.05) %>%
  pull(param) %>%
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

inp_file <- file("analyses/dtmr-results/meas-strc/meas-strc_model.inp")
writeLines(paste(print_mplus), con = inp_file, sep = "\n")
close(inp_file)

runModels(paste0(getwd(), "/analyses/dtmr-results/meas-strc"))

measstrc_params <- extractModelParameters(target = paste0(getwd(),
  "/analyses/dtmr-results/meas-strc"))[[1]] %>%
  as_data_frame() %>%
  filter(paramHeader == "New.Additional.Parameters")

write_csv(measstrc_params,
  "analyses/dtmr-results/meas-strc/meas-strc_params.csv")

rm(list = setdiff(ls(), c("strc_mplus", "check_length", "strc_params"))); gc()


### Structural-measurement reduction -------------------------------------------
remove <- strc_params %>%
  filter(grepl("L", param), !grepl("_0", param)) %>%
  filter((se == 0 | pval > 0.05) & !grepl("_0", param)) %>%
  pull(param)

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

inp_file <- file("analyses/dtmr-results/strc-meas/strc-meas_model.inp")
writeLines(paste(print_mplus), con = inp_file, sep = "\n")
close(inp_file)

runModels(paste0(getwd(), "/analyses/dtmr-results/strc-meas"))

strcmeas_params <- extractModelParameters(target = paste0(getwd(),
  "/analyses/dtmr-results/strc-meas"))[[1]] %>%
  as_data_frame() %>%
  filter(paramHeader == "New.Additional.Parameters")

write_csv(strcmeas_params,
  "analyses/dtmr-results/strc-meas/strc-meas_params.csv")

rm(list = ls()); gc()


### Combine results ------------------------------------------------------------
param_paths <- dir("analyses/dtmr-results", recursive = TRUE)[grep("params.csv",
  dir("analyses/dtmr-results", recursive = TRUE), fixed = TRUE)]

all_params <- map_df(param_paths, function(x) {
  reduction <- strsplit(x, "/") %>% unlist() %>% .[1]
  
  read_csv(paste0("analyses/dtmr-results/", x), col_types = cols(
    .default = col_number(), paramHeader = col_character(),
    param = col_character(), LatentClass = col_character()
  )) %>%
    rename(param_name = param) %>%
    mutate(
      reduction = reduction,
      param_type = case_when(grepl("G", param_name) ~ "strc", TRUE ~ "meas"),
      item = map_dbl(param_name, function(x) {
        if (grepl("G", x)) return(NA_real_)
        strsplit(x, "_") %>% unlist() %>% .[1] %>% gsub("L", "", x = .) %>%
          as.numeric()
      }),
      param = map_chr(param_name, function(x) {
        strsplit(x, "_") %>% unlist() %>% .[2]
      })
    ) %>%
    select(reduction, param_name, param_type, item, param, est, se, est_se,
      pval)
})

saveRDS(all_params, file = "data/dtmr-results.rds")
