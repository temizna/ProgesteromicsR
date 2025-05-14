# === Module: mod_sample_select ===

#' Sample Selection Module
#'
#' This module enables selection of samples from uploaded or demo data.
#' Users can view the design table and choose samples for downstream analyses.
#'
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param loaded_data_rv list storing raw_counts samples normailized counts and species
#' @param filtered_data_rv a reactive data with filtered tables.
#' @param filtered_dds_rv a reactive filtered dds object
#' @param dds_rv reactive value full dds object
#' @return Updates reactive values: `filtered_data_rv` (filtered data) and `filtered_dds_rv` (filtered DESeq2 object)
#' @importFrom DT renderDT datatable
#' @importFrom shinythemes shinytheme 
#' @export
mod_sample_select <- function(input, output, session, dds_rv, loaded_data_rv, filtered_data_rv, filtered_dds_rv) {
  observe({
    req(loaded_data_rv()$samples)
    samples<-loaded_data_rv()$samples
    updateSelectInput(session, "filter_dim", choices = unique(samples$dimension))
    updateSelectInput(session, "filter_treatment", choices = unique(samples$treatment))
    updateSelectInput(session, "filter_batch", choices = unique(samples$batch))
    updateSelectInput(session, "filter_cellline", choices = unique(samples$cellline))
    updateSelectInput(session, "filter_PR", choices = unique(samples$PR))
    updateSelectInput(session, "filter_ER", choices = unique(samples$ER))
    updateSelectInput(session, "sample_select", choices = rownames(samples),selected = NULL)
    #print(colnames(samples))
  })
  

  # Reactive expression for filtered data based on user input
  
  observeEvent(input$run_filter, {
    # No req() check for filters, just let them be empty or NULL
    req(loaded_data_rv(), dds_rv())
    
    # Print the head of loaded data
    print("Loaded Data:")
   # print(head(loaded_data_rv()))
    
    # Start with the full dataset
    filtering_data <- loaded_data_rv()
    
    # Print to verify the dataset being filtered
    #print("Filtering Data:")
    #print(head(filtering_data))
    
    print("Input filter values:")
    print(paste("filter_cellline: ", input$filter_cellline))
    print(paste("filter_treatment: ", input$filter_treatment))
    print(paste("filter_batch: ", input$filter_batch))
    print(paste("filter_dim: ", input$filter_dim))
    print(paste("filter_ER: ", input$filter_ER))
    print(paste("filter_PR: ", input$filter_PR))
    
    # Apply filters based on user input, but skip filters that are NULL or empty
    if (!is.null(input$filter_cellline) && length(input$filter_cellline) > 0) {
      print("Filtering based on cellline...")
      filtering_data$samples <- filtering_data$samples[filtering_data$samples$cellline %in% input$filter_cellline, ]
    }
    
    if (!is.null(input$filter_treatment) && length(input$filter_treatment) > 0) {
      print("Filtering based on treatment...")
      filtering_data$samples <- filtering_data$samples[filtering_data$samples$treatment %in% input$filter_treatment, ]
    }
    
    if (!is.null(input$filter_batch) && length(input$filter_batch) > 0) {
      print("Filtering based on batch...")
      filtering_data$samples <- filtering_data$samples[filtering_data$samples$batch %in% input$filter_batch, ]
    }
    
    if (!is.null(input$filter_dim) && length(input$filter_dim) > 0) {
      print("Filtering based on dimension...")
      filtering_data$samples <- filtering_data$samples[filtering_data$samples$dimension %in% input$filter_dim, ]
    }
    
    if (!is.null(input$filter_ER) && length(input$filter_ER) > 0) {
      print("Filtering based on ER...")
      filtering_data$samples <- filtering_data$samples[filtering_data$samples$ER %in% input$filter_ER, ]
    }
    
    if (!is.null(input$filter_PR) && length(input$filter_PR) > 0) {
      print("Filtering based on PR...")
      filtering_data$samples <- filtering_data$samples[filtering_data$samples$PR %in% input$filter_PR, ]
    }
    
      # Ensure we don't lose other information (counts, norm_counts, etc.)
      filtering_data$counts <- filtering_data$counts[, rownames(filtering_data$samples)]
      filtering_data$norm_counts <- filtering_data$norm_counts[, rownames(filtering_data$samples)]
      filtering_data$species <- loaded_data_rv()$species
      
      # Update the filtered dataset in filtered_data_rv
      filtered_data_rv(filtering_data)
      #print(head(filtered_data_rv()))
      # Now, create the filtered DESeq2 object based on the filtered samples
      if (!is.null(dds_rv())) {
        filtered_dds <- dds_rv()
        filtered_dds <- filtered_dds[, rownames(filtering_data$samples)]  # Subset dds to filtered samples
        filtered_dds_rv(filtered_dds)  # Store the filtered DESeq2 object
        #print(head(filtered_dds))
      }
    }
  )
  
  # Observe "Deselect All" button
  observeEvent(input$deselect_all, {
    req(loaded_data_rv(), filtered_data_rv(), dds_rv(), filtered_dds_rv())
    # Deselect all samples in "Select Samples" dropdown
    samples<-loaded_data_rv()$samples
    # Optionally, reset filtered data to the full dataset
    filtered_data_rv(loaded_data_rv())  # Reset to full data
    # Optionally, reset the filtered DESeq2 object to the full dds_rv
    filtered_dds_rv(dds_rv())  # Reset to full DESeq2 object
    updateSelectInput(session, "filter_dim", choices = unique(samples$dimension))
    updateSelectInput(session, "filter_treatment", choices = unique(samples$treatment))
    updateSelectInput(session, "filter_batch", choices = unique(samples$batch))
    updateSelectInput(session, "filter_cellline", choices = unique(samples$cellline))
    updateSelectInput(session, "filter_PR", choices = unique(samples$PR))
    updateSelectInput(session, "filter_ER", choices = unique(samples$ER))
    updateSelectInput(session, "sample_select", choices = rownames(samples),selected = NULL)
  })
  
  observeEvent(input$select_all, {
    req(loaded_data_rv(), dds_rv())
    # Deselect all samples in "Select Samples" dropdown
    samples<-loaded_data_rv()$samples
    # Optionally, reset filtered data to the full dataset
    filtered_data_rv(loaded_data_rv())  # Reset to full data
    # Optionally, reset the filtered DESeq2 object to the full dds_rv
    filtered_dds_rv(dds_rv())  # Reset to full DESeq2 object
  })
  # Example: Render filtered data table
  output$filteredDataTable <- renderDT({
    req(filtered_data_rv())
    datatable(filtered_data_rv()$samples)
  })
}


