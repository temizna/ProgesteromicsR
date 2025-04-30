#' Logger Module Server
#'
#' @param id Shiny module id
#' @param input_all The full Shiny input object (passed from main server)
#' @return A reactive expression with log entries
#' @importFrom shiny reactiveValuesToList moduleServer onStop req observeEvent renderText reactiveValues reactive
#' @export
mod_logger_server <- function(id, input_all) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    user_log <- reactiveValues(entries = character())
    
    observe({
      req(input_all)
      all_inputs <- reactiveValuesToList(input_all)
      
      for (name in names(all_inputs)) {
        observeEvent(input_all[[name]], {
          value <- input_all[[name]]
          if (!is.null(value)) {
            log_entry <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", name, ": ", paste(value, collapse = ", "))
            user_log$entries <- c(user_log$entries, log_entry)
          }
        }, ignoreInit = TRUE)
      }
    })
    
    output$log_text <- renderText({
      paste(user_log$entries, collapse = "\n")
    })
    
    output$download_log <- downloadHandler(
      filename = function() { paste0("user_session_log_", Sys.Date(), ".txt") },
      content = function(file) {
        writeLines(user_log$entries, con = file)
      }
    )
    
    session$onSessionEnded(function() {
      log_dir <- file.path(getwd(), "logs")
      if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
      log_file <- file.path(log_dir, paste0("user_session_log_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".txt"))
      writeLines(user_log$entries, log_file)
      message("[Logger] Saved session log to: ", log_file)
    })
    
    # Return reactive expression for use in other modules
    return(reactive(user_log$entries))
  })
}

#' Logger Module UI
#'
#' @param id Shiny module id
#' @return A Shiny UI element to display the log and a download button
#' @importFrom shiny NS tagList h3 verbatimTextOutput div downloadButton
#' @export
mod_logger_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("User Session Log"),
    div(
      style = "overflow-y: auto; height: 500px; border: 1px solid #ccc; padding: 10px; white-space: pre-wrap; background-color: #f9f9f9; font-family: monospace; font-size: 12px;",
      verbatimTextOutput(ns("log_text"))
    ),
    downloadButton(ns("download_log"), "Download Log")
  )
}
