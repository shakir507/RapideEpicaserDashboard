# Run with shiny::runApp("epicaser_dashboard")
library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(epicaser)

options(shiny.maxRequestSize = 50 * 1024^2)

# Local copy: provide the empty datetime constructor used by epicaser 0.3.
# The installed namespace and original engine function are unchanged.
generate_pathways_fixed <- epicaser::generate_patient_pathways_NL2
pathway_environment <- new.env(parent = environment(generate_pathways_fixed))
pathway_environment$POSIXct <- function() {
  as.POSIXct(character(), tz = "UTC")
}
environment(generate_pathways_fixed) <- pathway_environment

age_labels <- c("00-19", "20-39", "40-59", "60-79", "80-99", "100-110")
cohort_summary <- function(data, date_column) {
  # Use actual patient ages when available, rather than guessing cohort codes.
  if ("Age" %in% names(data) && is.numeric(data$Age)) {
    data$CohortGroup <- cut(data$Age, breaks = c(0, 20, 40, 60, 80, 100, 111),
                            labels = age_labels, right = FALSE)
    if (anyNA(data$CohortGroup)) stop("Patient ages fall outside the supported 0–110 age ranges.")
  } else {
    data$CohortGroup <- factor(data$CohortGroup)
  }
  count(data, .data[[date_column]], CohortGroup, name = "Cases")
}

check_columns <- function(x, required, stage) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(stage, ": missing columns: ", paste(missing, collapse = ", "))
  x
}
check_dates <- function(x, column) {
  x[[column]] <- as.Date(as.character(x[[column]]), format = "%Y-%m-%d")
  if (anyNA(x[[column]])) stop(column, " must contain valid YYYY-MM-DD dates.")
  x
}
curve <- function(data, x, y, group = NULL, title) {
  if (!is.null(group)) {
    values <- as.character(data[[group]])
    levels <- if (all(values %in% age_labels)) age_labels else unique(values)
    data[[group]] <- factor(values, levels = levels)
  }
  p <- if (is.null(group)) {
    ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
      geom_line(colour = "#176b87") + geom_point(colour = "#176b87", size = 1)
  } else {
    ggplot(data, aes(x = .data[[x]], y = .data[[y]], fill = .data[[group]])) +
      geom_area() + scale_fill_discrete(drop = FALSE)
  }
  # Explicit text styles avoid inheriting invisible axis/legend text.
  p + scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b", minor_breaks = NULL) +
    scale_y_continuous(labels = function(v) format(v, big.mark = ",", scientific = FALSE, trim = TRUE)) +
    labs(title = title,
         x = if (x == "DateAdmitted") "Admission date" else "Date",
         y = if (x == "DateAdmitted") "Hospital admissions per day" else "Reported cases per day",
         fill = "Age group (years)") +
    theme_minimal(base_size = 13) +
    theme(
      text = element_text(colour = "#222222", size = 13),
      axis.text.x = element_text(colour = "#222222", size = 11, angle = 45, hjust = 1),
      axis.text.y = element_text(colour = "#222222", size = 11),
      axis.title.x = element_text(colour = "#222222", size = 12, margin = margin(t = 10)),
      axis.title.y = element_text(colour = "#222222", size = 12, margin = margin(r = 10)),
      legend.text = element_text(colour = "#222222", size = 11),
      legend.title = element_text(colour = "#222222", size = 12),
      legend.position = "bottom",
      plot.title = element_text(colour = "#222222", size = 14, face = "bold"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      panel.grid.minor = element_blank(),
      plot.margin = margin(12, 16, 12, 16)
    ) + guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}

ui <- page_sidebar(
  title = "EpiCASER · Epidemic and hospital dashboard",
  theme = bs_theme(version = 5, bootswatch = "flatly"),
  sidebar = sidebar(width = 320,
                    radioButtons(
                      "source",
                      "Starting point",
                      c(
                        "Run model + generate stochastic cases" = "model",
                        "Upload deterministic CSV + generate stochastic cases" = "csv",
                        "Upload existing stochastic cases CSV" = "stochastic"
                      )
                    ),

                    conditionalPanel(
                      "input.source == 'model'",

                      numericInput(
                        "N", "Population (N)",
                        value = 100000, min = 1, step = 1000
                      ),

                      numericInput(
                        "I0", "Initial infected (I0)",
                        value = 10, min = 0, step = 1
                      ),
                      numericInput(
                        "end_time", "Last day of simulation",
                        value = 100, min = 10, step = 10
                      ),
                      dateInput(
                        inputId = "start_day",
                        label = "First date of simulation",
                        value = "2023-01-01"
                      ),
                      sliderInput(
                        "RF", "RF (Reporting Fraction)",
                        min = 0, max = 1, value = 0.3, step = 0.01
                      ),
                      sliderInput(
                        "R0", "R0 (Basic Reproduction number)",
                        min = 1, max = 10, value = 1.5, step = 0.1
                      )
                    ),

                    conditionalPanel(
                      "input.source == 'csv' || input.source == 'stochastic'",

                      fileInput(
                        "csv", "Upload CSV",
                        accept = ".csv"
                      ),

                      helpText(
                        "Deterministic: Day and Cases columns; keep all engine output columns.",
                        "Stochastic: Date and Cases columns with daily, non-negative integer counts.",
                        "Only stochastic uploads require dates in YYYY-MM-DD format."
                      )
                    ),

                    conditionalPanel(
                      "input.source == 'model' || input.source == 'csv'",

                      numericInput(
                        "n_runs", "Number of stochastic case-generation runs",
                        value = 1, min = 1, step = 1
                      )
                    ),
                    conditionalPanel("input.source != 'stochastic'",
                                     selectInput("measurement", "Measurement model", c("Negative binomial" = "nb", "Poisson" = "poisson"))
                    ),
                    numericInput("seed", "Random seed", 123, min = 0, max = 2147483647, step = 1),
                    checkboxInput("hospital", "Generate hospital arrivals and NL pathways", TRUE),
                    actionButton("run", "Run pipeline", class = "btn-primary"),
                    hr(),
                    helpText("Results update only when you click Run pipeline. Age profiles, hospital probabilities and pathway settings use your installed EpiCASER defaults."),
                    textOutput("status"),
                    selectInput("view_run", "View completed run", choices = c("Run 1" = "1"))
  ),
  layout_columns(
    value_box("Synthetic cases", textOutput("total")),
    value_box("Peak daily cases", textOutput("peak")),
    value_box("Hospital arrivals", textOutput("admissions"))
  ),
  navset_card_tab(
    nav_panel("Cases", plotOutput("cases_plot", height = "480px"), DTOutput("cases_table")),
    nav_panel("Age cohorts", plotOutput("age_plot", height = "480px"), DTOutput("age_table")),
    nav_panel("Case list check", plotOutput("list_plot", height = "480px"), DTOutput("list_table")),
    nav_panel("Hospital arrivals", plotOutput("hospital_plot", height = "480px"), DTOutput("hospital_table")),
    nav_panel(
      "Patient pathways",
      navset_tab(
        nav_panel(
          "Pathway stages",
          DTOutput("pathway_table")
        ),
        nav_panel(
          "Patient details + pathways",
          DTOutput("joined_table")
        )
      )
    ),    nav_panel("Export / run details", verbatimTextOutput("details"),
                    selectInput("export", "Dataset", c("NL pathways" = "pathways", "Joined hospital dataset" = "joined", "Deterministic output" = "deterministic", "Stochastic cases" = "cases", "Age cohorts" = "age", "Epidemiological case list" = "epi", "Hospital arrivals" = "arrivals")),
                    downloadButton("download", "Download CSV"),
                    helpText("The joined dataset follows your left_join(arrivals, pathways_NL). It contains hospital arrivals and their pathways, not all community cases.")
    )
  ),

  fluidRow(
    style = "
      text-align: center;
      margin-top: 30px;
      padding-top: 10px;
      border-top: 1px solid #ddd;
      color: #555;
      font-size: 13px;
    ",

    column(
      width = 6,
      tags$img(
        src = "RAPIDE_logo.png",
        height = "60px"
      ),
      div("RAPIDE Project"),
      div("WP3 Simulation Dashboard Supporting WP6 TTX")
    ),

    column(
      width = 6,
      tags$img(
        src = "EUfunding.png",
        height = "60px"
      ),
      div("Funded by the European Union"),
      div(
        style = "max-width: 700px; margin: auto;",
        "Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Health and Digital Executive Agency (HADEA). Neither the European Union nor the granting authority can be held responsible for them."
      )
    )
  )
)

server <- function(input, output, session) {
  result <- reactiveVal(NULL)
  completed_runs <- reactiveVal(list())
  status <- reactiveVal("Choose a starting point, then run the pipeline.")
  observeEvent(input$run, {
    # Clear previous results so an error cannot display stale outputs.
    result(NULL)
    completed_runs(list())
    status("Running…")
    tryCatch({
      if (!is.finite(input$seed) || input$seed != floor(input$seed)) stop("Use an integer random seed.")
      set.seed(as.integer(input$seed))
      withProgress(message = "EpiCASER pipeline", value = 0, {
        deterministic <- NULL
        uploaded_cases <- NULL
        if (input$source %in% c("csv", "stochastic")) {
          if (is.null(input$csv)) stop("Upload a CSV first.")
          uploaded <- read.csv(input$csv$datapath, check.names = FALSE,
                               stringsAsFactors = FALSE)
          message("USING UPLOADED CSV: ", input$csv$name,
                  " | Rows: ", nrow(uploaded), " | Mode: ", input$source)
          if (!nrow(uploaded)) stop("The CSV is empty.")
          if (input$source == "csv") {
            check_columns(uploaded, c("Day", "Cases"), "Deterministic CSV")
            if (!is.numeric(uploaded$Day) || any(!is.finite(uploaded$Day)) ||
                anyDuplicated(uploaded$Day)) {
              stop("Day must contain unique, finite numeric values.")
            }
            if (!is.numeric(uploaded$Cases) || any(!is.finite(uploaded$Cases)) ||
                any(uploaded$Cases < 0)) {
              stop("Deterministic Cases must contain finite, non-negative numeric values; decimals are allowed.")
            }
            # Preserve Day and every other engine column; no Date conversion here.
            deterministic <- uploaded[order(uploaded$Day), , drop = FALSE]
          } else {
            check_columns(uploaded, c("Date", "Cases"), "Stochastic CSV")
            uploaded <- check_dates(uploaded, "Date")
            if (anyDuplicated(uploaded$Date)) stop("Use one row per date.")
            uploaded_cases <- uploaded[order(uploaded$Date), , drop = FALSE]
          }
        } else if (input$source == "model") {
          vals <- c(input$N, input$I0, input$RF)
          if (any(!is.finite(vals)) || input$N < 1 || input$I0 < 0 || input$I0 > input$N ||
              input$N != floor(input$N) || input$I0 != floor(input$I0)) {
            stop("N and I0 must be integers with 0 <= I0 <= N and N >= 1.")
          }
          message("RUNNING DETERMINISTIC ENGINE")
          deterministic <- generate_epi_cases_deterministic(R0 = input$R0, RF = input$RF, N = input$N, I0 = input$I0,
                                                            start_time = 0, end_time = input$end_time)
        } else {
          stop("Choose a valid starting point.")
        }

        n_runs <- if (input$source == "stochastic") 1 else input$n_runs
        if (length(n_runs) != 1 || !is.finite(n_runs) || n_runs < 1 ||
            n_runs != floor(n_runs) || n_runs > .Machine$integer.max) {
          stop("Number of runs must be a positive integer.")
        }
        if (input$seed < 0 || input$seed + n_runs - 1 > .Machine$integer.max) {
          stop("Seed plus number of runs exceeds the supported seed range.")
        }
        run_results <- vector("list", as.integer(n_runs))
        for (run_index in seq_len(as.integer(n_runs))) {
          run_seed <- as.integer(input$seed + run_index - 1)
          set.seed(run_seed)
          setProgress(value = (run_index - 1) / n_runs,
                      detail = paste("Run", run_index, "of", n_runs))
          if (input$source == "stochastic") {
            cases <- uploaded_cases
            message("USING EXISTING STOCHASTIC CASES; stochastic generation skipped")
          } else {
            cases <- generate_epi_cases_stochastic(deterministic,
                                                   Poisson = input$measurement == "poisson",start_day =input$start_day)
          }
          check_columns(cases, c("Date", "Cases"), "Stochastic output")
          cases <- check_dates(cases, "Date")
          if (!nrow(cases) || !is.numeric(cases$Cases) || anyNA(cases$Cases) || any(!is.finite(cases$Cases)) || any(cases$Cases < 0) || any(cases$Cases != floor(cases$Cases))) stop("Stochastic output must contain non-negative integer Cases.")
          incProgress(0.2 / n_runs, detail = "Age cohorts")
          age <- generate_cohort_epi_cases(cases$Date, cases$Cases)
          check_columns(age, c("Date", "00-19", "100-110"), "Cohort output")
          age <- check_dates(age, "Date")
          age_tidy <- pivot_longer(age, cols = `00-19`:`100-110`, names_to = "Age", values_to = "Incidence")
          incProgress(0.2 / n_runs, detail = "Individual case list")
          epi <- generate_epi_case_list(age)
          check_columns(epi, c("Date", "CohortGroup"), "Case list")
          epi <- check_dates(epi, "Date")
          epi_sum <- cohort_summary(epi, "Date")
          arrivals <- pathways <- joined <- arr_sum <- NULL
          hospital_error <- NULL
          if (isTRUE(input$hospital)) {
            # Preserve successful community outputs if a hospital stage fails.
            tryCatch({
              incProgress(0.1 / n_runs, detail = "Hospital arrivals")
              arrivals <- generate_hospital_arrivals(epi)
              check_columns(arrivals, c("DateAdmitted", "CohortGroup", "CaseID"), "Hospital arrivals")
              arrivals <- check_dates(arrivals, "DateAdmitted")
              arr_sum <- cohort_summary(arrivals, "DateAdmitted")
              if (nrow(arrivals) > 0) {
                incProgress(0.2 / n_runs, detail = "NL patient pathways")
                pathways <- generate_pathways_fixed(arrivals, set_seed = FALSE)
                check_columns(pathways, c("admission", "startdatetime", "enddatetime"), "NL pathways")
                # Use the engine's timestamps directly, exactly as in your script.
                pathways$Duration <- round(as.numeric(difftime(pathways$enddatetime, pathways$startdatetime, units = "days")), 2)
                joined <- left_join(arrivals, pathways, by = c("CaseID" = "admission"), keep = TRUE)
              }
            }, error = function(e) hospital_error <<- conditionMessage(e))
          }
          run_results[[run_index]] <- list(deterministic = deterministic, cases = cases, age = age, age_tidy = age_tidy,
                                           epi = epi, epi_sum = epi_sum, arrivals = arrivals, arr_sum = arr_sum,
                                           pathways = pathways, joined = joined,
                                           settings = list(source = input$source, file = if (input$source %in% c("csv", "stochastic")) input$csv$name else NULL,
                                                           model_parameters = if (input$source == "model") list(N = input$N, I0 = input$I0, RF = input$RF) else NULL,
                                                           measurement = if (input$source == "stochastic") "Already supplied in CSV" else input$measurement,
                                                           seed = run_seed, run_index = run_index, n_runs = n_runs,
                                                           hospital = input$hospital,
                                                           engine_version = as.character(packageVersion("epicaser")),
                                                           run_time = as.character(Sys.time()), hospital_error = hospital_error))
          setProgress(value = run_index / n_runs)
        }
        completed_runs(run_results)
        result(run_results[[1]])
        updateSelectInput(session, "view_run",
                          choices = setNames(as.character(seq_len(n_runs)), paste("Run", seq_len(n_runs))),
                          selected = "1")
        failed <- vapply(run_results, function(r) !is.null(r$settings$hospital_error), logical(1))
        status(paste("Completed", n_runs, "run(s).",
                     if (any(failed)) "Some hospital stages failed; see each run's details." else "Select a run to inspect its results."))
        if (any(failed)) showNotification("Some hospital stages failed; see run details.", type = "error", duration = NULL)
      })
    }, error = function(e) {
      status(paste("Run failed:", conditionMessage(e)))
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })
  observeEvent(input$view_run, {
    runs <- completed_runs()
    index <- suppressWarnings(as.integer(input$view_run))
    if (length(runs) && length(index) == 1 && !is.na(index) && index >= 1 && index <= length(runs)) {
      result(runs[[index]])
    }
  })
  observeEvent(result(), {
    r <- result()
    if (is.null(r)) return()
    available <- c("NL pathways" = "pathways",
                   "Joined hospital dataset" = "joined",
                   "Deterministic output" = "deterministic", "Stochastic cases" = "cases",
                   "Age cohorts" = "age", "Epidemiological case list" = "epi",
                   "Hospital arrivals" = "arrivals")
    available <- available[vapply(unname(available), function(key) !is.null(r[[key]]), logical(1))]
    chosen <- if (input$export %in% available) input$export else "cases"
    updateSelectInput(session, "export", choices = available, selected = chosen)
  })
  output$status <- renderText(status())
  output$total <- renderText({ if (is.null(result())) "—" else format(sum(result()$cases$Cases), big.mark = ",") })
  output$peak <- renderText({ if (is.null(result())) "—" else format(max(result()$cases$Cases), big.mark = ",") })
  output$admissions <- renderText({ if (is.null(result()) || is.null(result()$arrivals)) "—" else format(nrow(result()$arrivals), big.mark = ",") })
  output$cases_plot <- renderPlot({ req(result()); curve(result()$cases, "Date", "Cases", title = "Synthetic daily cases") })
  output$age_plot <- renderPlot({ req(result()); curve(result()$age_tidy, "Date", "Incidence", "Age", "Cases by age cohort") })
  output$list_plot <- renderPlot({ req(result()); curve(result()$epi_sum, "Date", "Cases", "CohortGroup", "Case list summary") })
  output$hospital_plot <- renderPlot({
    req(result())
    validate(need(!is.null(result()$arr_sum), "Enable hospital generation and run the pipeline."),
             need(nrow(result()$arr_sum) > 0, "No hospital arrivals in this run."))
    curve(result()$arr_sum, "DateAdmitted", "Cases", "CohortGroup", "Daily hospital arrivals")
  })
  table_map <- c(cases_table = "cases", age_table = "age", list_table = "epi_sum", hospital_table = "arrivals", pathway_table = "pathways", joined_table = "joined")
  for (id in names(table_map)) local({
    output_id <- id
    dataset <- table_map[[id]]

    output[[output_id]] <- renderDT({
      req(result())
      data <- result()[[dataset]]

      validate(need(
        !is.null(data),
        "This dataset has not been generated. Check run status."
      ))

      # Format datetime columns for display only
      for (column in names(data)) {
        if (inherits(data[[column]], "POSIXt")) {
          data[[column]] <- format(
            data[[column]],
            format = "%Y-%m-%d %H:%M:%S",
            tz = "UTC"
          )
        }
      }

      datatable(
        data,
        rownames = FALSE,
        options = list(pageLength = 4, scrollX = TRUE)
      )
    }, server = TRUE)
  })
  output$details <- renderPrint({ req(result()); result()$settings })
  output$download <- downloadHandler(
    filename = function() paste0("epicaser_", input$export, "_run_", result()$settings$run_index, ".csv"),
    content = function(file) {
      req(result())
      data <- result()[[input$export]]
      if (is.null(data)) stop("Selected dataset is unavailable. Check run status or select another dataset.")
      write.csv(data, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
