# Rescorla-Wagner Shiny App (single-file)
# Features:
# - dynamic CS editor (add/remove CS, per-phase alpha and initial V)
# - trial-type editor (checkbox UI for CS presence + lambda, add/remove types)
# - phase editor (add/remove phases, select multiple trial-types per phase, set n trials & randomization)
# - trial generation from phases and trial-types
# - simulation with phase-dependent alphas (per CS)
# - plot with markers for cue presence, export & save/load

# (C) Annika Ziereis 2025, annika.ziereis@uni-goettingen.de, annika.ziereis@web.de

if (!require("pacman")) install.packages("pacman")
pacman::p_load("shiny", "DT", "dplyr", "tidyr", "ggplot2", "purrr", "openxlsx", "jsonlite")


library(shiny)
library(DT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(openxlsx)
library(jsonlite)

# -------- helper utilities --------
parse_alpha <- function(a_raw){
  if(is.null(a_raw)) return(list(general = NA_real_))
  a_raw <- trimws(as.character(a_raw))
  if (a_raw == "") return(list(general = NA_real_))
  # numeric simple
  if (grepl('^[0-9.]+$', a_raw)) return(list(general = as.numeric(a_raw)))
  # JSON
  if (startsWith(a_raw, "{") && endsWith(a_raw, "}")){
    try({ x <- jsonlite::fromJSON(a_raw); return(as.list(x)) }, silent = TRUE)
  }
  # try eval R code
  try({ parsed <- eval(parse(text = a_raw));
  if (is.numeric(parsed) && !is.null(names(parsed))) return(as.list(parsed))
  if (is.list(parsed)) return(parsed)
  }, silent = TRUE)
  # fallback try single number
  num <- suppressWarnings(as.numeric(a_raw))
  if (!is.na(num)) return(list(general = num))
  return(list(general = NA_real_))
}

normalize_cs_list <- function(cs_list){
  lapply(cs_list, function(x){
    if(is.null(x$V0) && !is.null(x$v0)) x$V0 <- x$v0
    if(is.null(x$V0)) x$V0 <- 0
    if(is.null(x$alpha)) x$alpha <- list(general = 0.5)
    # if character -> parse
    if (is.character(x$alpha)) x$alpha <- parse_alpha(x$alpha)
    # scalar numeric -> wrap
    if (is.numeric(x$alpha) && is.null(names(x$alpha))) x$alpha <- list(general = as.numeric(x$alpha))
    # named numeric -> to list
    if (is.numeric(x$alpha) && !is.null(names(x$alpha))) x$alpha <- as.list(x$alpha)
    # ensure values numeric where possible
    if (is.list(x$alpha)) x$alpha <- lapply(x$alpha, function(v){ if(is.character(v)) { nv <- suppressWarnings(as.numeric(v)); if(!is.na(nv)) return(nv) else return(v)} else return(v) })
    x
  })
}

# generate trials from phases structure and trialtypes list
build_types_df_from_list <- function(tt_list, cs_names){
  if(is.null(tt_list) || length(tt_list)==0) return(NULL)
  rows <- lapply(names(tt_list), function(nm){
    tt <- tt_list[[nm]]
    cs_present <- setNames(rep(0, length(cs_names)), cs_names)
    if(!is.null(tt$cs_present)){
      for(cs in intersect(names(tt$cs_present), cs_names)) cs_present[cs] <- as.numeric(isTRUE(tt$cs_present[[cs]]))
    }
    df <- data.frame(name = nm, lambda = as.numeric(ifelse(is.null(tt$lambda), 0, tt$lambda)), stringsAsFactors = FALSE)
    for(cs in cs_names) df[[cs]] <- cs_present[cs]
    df
  })
  bind_rows(rows)
}

generate_trials <- function(phases_list, types_df, cs_names){
  # phases_list: named list with each element list(n_trials, trialtypes (character vector), randomize logical)
  out <- list(); counter <- 1
  # convert types_df to have rownames by name
  if(is.null(types_df) || nrow(types_df)==0) return(NULL)
  for (pname in names(phases_list)){
    ph <- phases_list[[pname]]
    nrep <- as.integer(ph$n_trials)
    chosen <- ph$trialtypes
    if(is.null(chosen) || length(chosen)==0) next
    block_types <- types_df %>% filter(name %in% chosen)
    if(nrow(block_types)==0) next
    if(isTRUE(ph$randomize)){
      idx <- sample(seq_len(nrow(block_types)), nrep, replace = TRUE)
    } else {
      idx <- rep(seq_len(nrow(block_types)), length.out = nrep)
    }
    tmp <- block_types[idx, , drop = FALSE]
    tmp$phase <- pname
    tmp$trial <- counter:(counter + nrep - 1)
    out[[length(out)+1]] <- tmp
    counter <- counter + nrep
  }
  if(length(out)==0) return(NULL)
  trials <- bind_rows(out)
  # ensure cs columns exist
  for(cs in cs_names) if(!(cs %in% names(trials))) trials[[cs]] <- 0
  existing_cs <- intersect(cs_names, names(trials))
  trials_df <- trials %>% select(trial, phase, lambda, all_of(existing_cs))
  # rename existing_cs -> <cs>_present
  rename_map <- setNames(existing_cs, paste0(existing_cs, "_present"))
  trials_df <- trials_df %>% rename(!!!rename_map)
  return(trials_df)
}

# simulate function
simulate_rw_trials <- function(trials_df, cs_list, beta = 0.5){
  cs_list <- normalize_cs_list(cs_list)
  cs_names <- names(cs_list)
  V <- map_dbl(cs_list, ~ .x$V0)
  n_trials <- nrow(trials_df)
  V_hist <- matrix(NA, n_trials + 1, length(cs_names))
  colnames(V_hist) <- cs_names
  V_hist[1,] <- V
  get_alpha <- function(cs, phase){
    cs_alpha <- cs_list[[cs]]$alpha
    if(is.numeric(cs_alpha)) return(cs_alpha)
    if(phase %in% names(cs_alpha)) return(cs_alpha[[phase]])
    if('general' %in% names(cs_alpha)) return(cs_alpha$general)
    return(cs_alpha[[1]])
  }
  for(t in 1:n_trials){
    lambda_t <- trials_df$lambda[t]
    phase_t <- trials_df$phase[t]
    present <- sapply(cs_names, function(cs){
      val <- trials_df[[paste0(cs, '_present')]][t]
      if(is.null(val) || is.na(val)) return(0)
      return(as.numeric(val))
    })
    alpha_t <- sapply(cs_names, get_alpha, phase = phase_t)
    V_pred <- sum(V * present)
    delta <- lambda_t - V_pred
    V <- V + alpha_t * beta * delta * present
    V_hist[t+1, ] <- V
  }
  # add trial0 row (all absent)
  present_cols <- grep('_present$', names(trials_df), value = TRUE)
  trial0 <- trials_df[1, , drop = FALSE]
  trial0$trial <- 0
  trial0[present_cols] <- 0
  trials_df2 <- bind_rows(trial0, trials_df)
  df_out <- bind_cols(trials_df2, as_tibble(V_hist))
  df_out <- df_out %>% rowwise() %>% mutate(net_V = sum(c_across(all_of(cs_names)), na.rm = TRUE)) %>% ungroup()
  df_long <- df_out %>% pivot_longer(all_of(cs_names), names_to = 'CS', values_to = 'V')
  df_out$phase <- factor(df_out$phase, levels = unique(trials_df2$phase))
  df_long$phase <- factor(df_long$phase, levels = unique(trials_df2$phase))
  presence_long <- df_out %>% select(trial, phase, ends_with('_present')) %>% pivot_longer(ends_with('_present'), names_to = 'CS', values_to = 'present') %>% mutate(CS = gsub('_present','',CS)) %>% filter(present==1)
  p <- ggplot(df_long, aes(x=trial, y=V, color = CS)) +
    geom_line(size=1) +
    facet_wrap(~phase, scales = 'free_x', space = 'free_x') +
    geom_line(aes(y = net_V), color='black', linetype=2) +
    geom_point(data = presence_long, aes(x=trial, y=0, color = CS), inherit.aes = FALSE, position = position_dodge2(width=0.6, padding = 0.1), size = 2) +
    theme_minimal() + labs(title='Rescorla-Wagner', y='V', x='Trial')
  return(list(plot = p, outtable = df_out))
}

# ------------------ UI ------------------
ui <- fluidPage(
  titlePanel('Rescorla-Wagner Simulator (flexible)'),
  sidebarLayout(
    sidebarPanel(width = 3,
                 h4('Global settings'),
                 numericInput('beta','beta (US learning rate)', value = 0.5, min=0, step=0.01),
                 textInput('cs_input','CS names (comma separated)', value = 'A,X'),
                 actionButton('set_cs','Create CS'),
                 hr(),
                 h4('Trial type editor'),
                 actionButton('add_type','Add trial type'),
                 br(), br(),
                 fileInput('load_experiment','Load experiment (.rds)', accept = '.rds'),
                 downloadButton('save_experiment','Save experiment (.rds)'),
                 hr(),
                 downloadButton('export_xlsx','Export results (.xlsx)')
    ),

    mainPanel(
      tabsetPanel(
        tabPanel('CS editor',
                 DTOutput('cs_table'),
                 helpText('Edit per-CS alpha for phases as JSON-like: e.g. {"learning":0.8,"extinction":0.1} or a single number 0.5 for general')
        ),
        tabPanel('Trial Types',
                 fluidRow(
                   column(4, textInput('new_tt_name','New trial type name')),
                   column(2, numericInput('new_lambda','λ', value = 1, min = 0, max = 1)),
                   column(2, actionButton('add_trialtype','Add Trial Type'))
                 ),
                 uiOutput('trialtype_list_ui')
        ),
        tabPanel('Phases',
                 fluidRow(
                   column(4, textInput('new_phase','New phase name')),
                   column(2, numericInput('new_phase_trials','# trials', 10, min=1)),
                   column(2, checkboxInput('new_phase_randomize','Randomize order', TRUE)),
                   column(2, actionButton('add_phase','Add Phase'))
                 ),
                 uiOutput('phase_list_ui')
        ),
        tabPanel('Preview Trials',
                 actionButton('gen_trials','Generate trials'),
                 DTOutput('trials_preview')
        ),
        tabPanel('Run',
                 actionButton('run_sim','Run simulation'),
                 plotOutput('rw_plot', height = '500px'),
                 DTOutput('out_table')
        )
      )
    )
  )
)

# ------------------ Server ------------------
server <- function(input, output, session){
  cs_names_r <- reactiveVal(NULL)
  rv <- reactiveValues(cs_df = NULL, trialtypes = list(), phases = list(), phases_df = NULL, trials_df = NULL, sim_result = NULL)

  observeEvent(input$set_cs, {
    names <- strsplit(input$cs_input, ',')[[1]] %>% trimws()
    cs_names_r(names)
    df_cs <- data.frame(name = names, alpha = rep('0.5', length(names)), V0 = rep(0, length(names)), stringsAsFactors = FALSE)
    rv$cs_df <- df_cs
  })

  output$cs_table <- renderDT({
    if(is.null(rv$cs_df)) return(datatable(data.frame()))
    datatable(rv$cs_df, editable = TRUE)
  })

  observeEvent(input$cs_table_cell_edit, {
    info <- input$cs_table_cell_edit
    i <- info$row; j <- info$col; v <- info$value
    rv$cs_df[i, j] <<- v
    replaceData(proxy = dataTableProxy('cs_table'), rv$cs_df, resetPaging = FALSE)
  })

  # trial types (checkbox UI)
  observeEvent(input$add_trialtype, {
    tt <- rv$trialtypes
    names_cs <- cs_names_r()
    nm <- input$new_tt_name
    if(is.null(nm) || nm == '' || nm %in% names(tt)) { showNotification('Invalid or duplicate trial type name', type = 'error'); return() }
    new <- list(lambda = as.numeric(input$new_lambda), cs_present = setNames(as.list(rep(FALSE, length(names_cs))), names_cs))
    tt[[nm]] <- new
    rv$trialtypes <- tt
  })

  output$trialtype_list_ui <- renderUI({
    tt <- rv$trialtypes
    cs_names <- cs_names_r()
    if(length(tt) == 0) return(tags$em('No trial types defined yet'))
    tagList(lapply(names(tt), function(nm){
      fluidRow(
        column(12, style='border:1px solid #ccc; padding:10px; margin-bottom:8px;',
               strong(nm),
               checkboxInput(paste0('lambda_', nm), 'λ present?', value = isTRUE(tt[[nm]]$lambda == 1)),
               lapply(cs_names, function(cs) checkboxInput(paste0('cs_', nm, '_', cs), label = cs, value = isTRUE(tt[[nm]]$cs_present[[cs]]))),
               actionButton(paste0('del_tt_', nm), 'Delete', class = 'btn-danger')
        )
      )
    }))
  })

  # sync UI inputs back into reactive list
  observe({
    tt <- rv$trialtypes
    cs_names <- cs_names_r()
    if(length(tt) == 0) return()
    for(nm in names(tt)){
      lam <- input[[paste0('lambda_', nm)]]
      if(!is.null(lam)) tt[[nm]]$lambda <- ifelse(lam, 1, 0)
      for(cs in cs_names){
        chk <- input[[paste0('cs_', nm, '_', cs)]]
        if(!is.null(chk)) tt[[nm]]$cs_present[[cs]] <- chk
      }
    }
    rv$trialtypes <- tt
  })

  # delete trial type
  observe({
    tt <- rv$trialtypes
    for(nm in names(tt)){
      btn <- input[[paste0('del_tt_', nm)]]
      if(!is.null(btn) && btn > 0){ tt[[nm]] <- NULL; rv$trialtypes <- tt }
    }
  })

  # phases
  observeEvent(input$add_phase, {
    ph <- rv$phases
    nm <- input$new_phase
    if(is.null(nm) || nm == '' || nm %in% names(ph)) { showNotification('Invalid or duplicate phase name', type = 'error'); return() }
    ph[[nm]] <- list(n_trials = as.integer(input$new_phase_trials), trialtypes = character(0), randomize = isTRUE(input$new_phase_randomize))
    rv$phases <- ph
  })

  output$phase_list_ui <- renderUI({
    ph <- rv$phases
    tt_names <- names(rv$trialtypes)
    if(length(ph) == 0) return(tags$em('No phases defined.'))
    tagList(lapply(names(ph), function(nm){
      fluidRow(
        column(12, style='border:1px solid #ccc; padding:10px; margin-bottom:8px;',
               strong(nm),
               numericInput(paste0('ntrial_', nm), '# trials', value = ph[[nm]]$n_trials, min = 1),
               checkboxInput(paste0('rand_', nm), 'Randomize order', value = isTRUE(ph[[nm]]$randomize)),
               checkboxGroupInput(paste0('ph_tt_', nm), 'Use trial types:', choices = tt_names, selected = ph[[nm]]$trialtypes, inline = TRUE),
               actionButton(paste0('del_phase_', nm), 'Delete', class = 'btn-danger')
        )
      )
    }))
  })

  observe({
    ph <- rv$phases
    tt_names <- names(rv$trialtypes)
    if(length(ph) == 0) return()
    for(nm in names(ph)){
      nt <- input[[paste0('ntrial_', nm)]]
      if(!is.null(nt)) ph[[nm]]$n_trials <- as.integer(nt)
      sel <- input[[paste0('ph_tt_', nm)]]
      if(!is.null(sel)) ph[[nm]]$trialtypes <- sel
      rnd <- input[[paste0('rand_', nm)]]
      if(!is.null(rnd)) ph[[nm]]$randomize <- rnd
      del <- input[[paste0('del_phase_', nm)]]
      if(!is.null(del) && del > 0){ ph[[nm]] <- NULL; rv$phases <- ph; return() }
    }
    rv$phases <- ph
  })

  # keep a small legacy types_df for saving/loading compatibility
  observe({
    if(is.null(rv$types_df) && length(rv$trialtypes) > 0){
      csn <- cs_names_r()
      tt <- rv$trialtypes
      rows <- lapply(names(tt), function(nm){
        r <- data.frame(name = nm, lambda = as.numeric(tt[[nm]]$lambda), stringsAsFactors = FALSE)
        for(cs in csn) r[[cs]] <- as.numeric(tt[[nm]]$cs_present[[cs]])
        r
      })
      rv$types_df <- bind_rows(rows)
    }
  })

  # preview trials: build types dataframe and phases dataframe then generate
  observeEvent(input$gen_trials, {
    cs_names <- cs_names_r()
    if(is.null(cs_names) || length(cs_names)==0){ showNotification('Define CS first', type='error'); return() }
    # build types df from rv$trialtypes
    types2 <- build_types_df_from_list(rv$trialtypes, cs_names)
    if(is.null(types2) || nrow(types2)==0){ showNotification('No trial types defined', type='error'); return() }
    # build phases_list from rv$phases
    ph_list <- rv$phases
    if(length(ph_list)==0){ showNotification('No phases defined', type='error'); return() }
    trials_df <- generate_trials(ph_list, types2, cs_names)
    if(is.null(trials_df)){ showNotification('No trials generated (check phases & selections)', type='warning'); return() }
    rv$trials_df <- trials_df
    output$trials_preview <- renderDT(datatable(trials_df, options = list(pageLength = 25)))
  })

  # run simulation
  observeEvent(input$run_sim, {
    if(is.null(rv$trials_df)){ showNotification('Generate trials first', type='error'); return() }
    if(is.null(rv$cs_df)){ showNotification('Define CS in CS editor', type='error'); return() }
    # build cs_list from cs_df
    csdf <- rv$cs_df
    cs_names <- cs_names_r()
    cs_list <- list()
    for(nm in cs_names){
      row <- csdf[csdf$name == nm, , drop = FALSE]
      if(nrow(row)==0){ cs_list[[nm]] <- list(alpha = list(general = 0.5), V0 = 0) }
      else {
        parsed <- parse_alpha(row$alpha[1])
        cs_list[[nm]] <- list(alpha = parsed, V0 = as.numeric(row$V0[1]))
      }
    }
    res <- simulate_rw_trials(rv$trials_df, cs_list, beta = input$beta)
    rv$sim_result <- res
    output$rw_plot <- renderPlot(res$plot)
    output$out_table <- renderDT(datatable(res$outtable, options = list(pageLength = 25)))
  })

  # save/load
  output$save_experiment <- downloadHandler(filename = function(){ paste0('rwm_experiment_', Sys.Date(), '.rds') },
                                            content = function(file){ saveRDS(list(cs = rv$cs_df, trialtypes = rv$trialtypes, phases = rv$phases), file) })
  observeEvent(input$load_experiment, {
    req(input$load_experiment)
    obj <- readRDS(input$load_experiment$datapath)
    rv$cs_df <- obj$cs
    rv$trialtypes <- obj$trialtypes
    rv$phases <- obj$phases
    # update cs table proxy
    replaceData(dataTableProxy('cs_table'), rv$cs_df, resetPaging = FALSE)
    showNotification('Experiment loaded', type = 'message')
  })

  # export xlsx
  output$export_xlsx <- downloadHandler(filename = function(){ paste0('rwm_results_', Sys.Date(), '.xlsx') },
                                        content = function(file){ req(rv$sim_result); wb <- createWorkbook(); addWorksheet(wb,'trials'); writeData(wb,'trials', rv$sim_result$outtable); addWorksheet(wb,'meta'); writeData(wb,'meta', data.frame(beta = input$beta, date = Sys.time())); saveWorkbook(wb,file, overwrite = TRUE) })
}

shinyApp(ui, server)
