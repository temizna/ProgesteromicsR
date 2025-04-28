# === Module: mod_tf_enrichment_analysis ===

#' Transcription Factor Enrichment Analysis Module
#'
#' This module allows users to perform transcription factor enrichment analysis
#' based on selected datasets stored locally (gzipped) in inst/extdata.
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param geneList_rv Reactive value containing the gene list (log2FC vector) for enrichment
#' @param tf_enrichment_result reactive data containing the tf enrichment results
#' @importFrom gson read.gmt
#' @importFrom vroom vroom
#' @export
mod_tf_enrichment_analysis <- function(input, output, session, geneList_rv, tf_enrichment_result) {
  
  load_tf_data <- function(tf_data_source) {
    tf_data_file_map <- list(
      "TRANSFAC and JASPAR PWMs" = "TRANSFAC_and_JASPAR_PWMs.h.gmt.gz",
      "ENCODE and ChEA Consensus TFs from ChIP-X" = "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt.gz",
      "TRRUST_Transcription_Factors_2019" = "TRRUST_Transcription_Factors_2019h.gmt.gz",
      "TF_Perturbations_Followed_by_Expression" = "TF_Perturbations_Followed_by_Expression_Human.gmt.gz",
      "hTFtarget" = "TF-Target-information.txt.gz",
      "TFLink" = "TFLink_Homo_sapiens_interactions_All_GMT_proteinName_v1.0.gmt.gz"
    )
    
    tf_data_file <- tf_data_file_map[[tf_data_source]]
    if (is.null(tf_data_file)) stop("Invalid TF data source selected.")
    
    tf_file_path <- system.file("extdata", tf_data_file, package = "ProgesteromicsR")
    if (tf_file_path == "") stop("TF data file not found.")
    
    if (grepl(".gmt", tf_data_file)) {
      tf_data <- gson::read.gmt(tf_file_path)
    } else {
      tf_data <- vroom::vroom(tf_file_path, delim = "\t", show_col_types = FALSE)
    }
    
    return(tf_data)
  }
  
  # Reactive value to store enrichment results
  tf_enrichment_results_rv <- reactiveVal()
  processed_tf_data <- reactiveVal()
  
  # Observe button click
  observeEvent(input$run_tf_enrichment, {
    req(geneList_rv())
    print("Running TF Enrichment Analysis...")
    tf_data_source <- input$tf_data_source
    print(paste("Selected TF data source:", tf_data_source))
    
    tryCatch({
      # Load TF data using helper function
      tf_data <- load_tf_data(tf_data_source)
      
      print("TF data loaded:")
      print(head(tf_data))
      
      # Process TF data: always map target genes to ENTREZ IDs
      tf_data <- tf_data %>%
        rename(TF_name = 1, target_gene = 2)
      
      tf_data_entrez <- bitr(
        tf_data$target_gene,
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db
      )
      
      tf_data_gmt <- tf_data %>%
        left_join(tf_data_entrez, by = c("target_gene" = "SYMBOL")) %>%
        select(TF_name, ENTREZID) %>%
        filter(!is.na(ENTREZID))
      
      print("Processed TERM2GENE:")
      print(head(tf_data_gmt))
      
      print("Processed GENELIST:")
      print(head(geneList_rv()))
      print("TF QVaL:")
      print(input$tf.qval)
      # Perform enrichment analysis
      tf_result <- clusterProfiler::enricher(
        gene = names(geneList_rv()),
        TERM2GENE = tf_data_gmt,
        pvalueCutoff = 0.05,
        qvalueCutoff = input$tf.qval
      )
      
      if (is.null(tf_result) || nrow(as.data.frame(tf_result)) < 1) {
        showNotification("No enriched terms found in TF enrichment.", type = "warning")
        return()
      }
      
      tf_enrichment_results_rv(tf_result)
      
      output$tf_dotplot <- renderPlot({
        enrichplot::dotplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
      })
      
      output$download_tf_dotplot <- downloadHandler(
        filename = function() paste0("TF_Enrichment_", input$tf_data_source, "_dotplot.pdf"),
        content = function(file) {
          pdf(file)
          print(enrichplot::dotplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
          dev.off()
        }
      )
      
      # output$tf_treeplot <- renderPlot({
      #   enrichplot::treeplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
      # })
      # 
      # output$download_tf_treeplot <- downloadHandler(
      #   filename = function() paste0("TF_Enrichment_", input$tf_data_source, "_treeplot.pdf"),
      #   content = function(file) {
      #     pdf(file)
      #     print(enrichplot::treeplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
      #     dev.off()
      #   }
      # )
      # 
      # output$tf_heatmap <- renderPlot({
      #   enrichplot::heatplot(tf_result, foldChange = geneList_rv(), showCategory = 10)
      # })
      # 
      # output$download_tf_heatmap <- downloadHandler(
      #   filename = function() paste0("TF_Enrichment_", input$tf_data_source, "_heatmap.pdf"),
      #   content = function(file) {
      #     pdf(file)
      #     print(enrichplot::heatplot(tf_result, foldChange = geneList_rv(), showCategory = 10))
      #     dev.off()
      #   }
      # )
      # 
      # output$tf_upsetplot <- renderPlot({
      #   enrichplot::upsetplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold"))
      # })
      # 
      # output$download_tf_upsetplot <- downloadHandler(
      #   filename = function() paste0("TF_Enrichment_", input$tf_data_source, "_upsetplot.pdf"),
      #   content = function(file) {
      #     pdf(file)
      #     print(enrichplot::upsetplot(tf_result) + theme(axis.text.y = element_text(size = 6, face = "bold")))
      #     dev.off()
      #   }
      # )
      
      output$tf_results_table <- renderDT({
        as.data.frame(tf_result@result)
      })
      
      output$download_tf_results_table <- downloadHandler(
        filename = function() paste0("TF_Enrichment_", input$tf_data_source, "_results.csv"),
        content = function(file) {
          write.csv(as.data.frame(tf_result@result), file, row.names = FALSE)
        }
      )
      
    }, error = function(e) {
      showNotification(paste("Error loading or processing TF data:", e$message), type = "error")
    })
  })
}
