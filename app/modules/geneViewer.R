# modules/geneViewer.R

# =====================
# UI
# =====================
geneViewerUI <- function(id){
  
  ns <- NS(id)
  
  div(
    class = "content",
    
    h2("🧬 | GENE VIEWER"),
    
    # ===== FILTER =====
    div(
      class = "filter-box",
      
      textInput(
        ns("gene"),
        "Filter by gene:",
        placeholder = "e.g. DPM1"
      )
    ),
    
    # ===== SUMMARY =====
    div(
      class = "table-box",
      
      withSpinner(
        uiOutput(ns("gene_summary")),
        type = 4,
        color = "#8b1e5b"
      )
    ),
    
    br(),
    
    # ===== INFO =====
    div(
      class = "table-box",
      
      withSpinner(
        uiOutput(ns("gene_info")),
        type = 4,
        color = "#8b1e5b"
      )
    ),
    
    br(),
    
    # ===== LINKS + UCSC =====
    div(
      style = "display:flex; gap:20px; align-items:stretch; flex-wrap:wrap;",
      
      div(
        class = "table-box",
        style = "flex:1; min-width:300px;",
        
        withSpinner(
          uiOutput(ns("external_links")),
          type = 4,
          color = "#8b1e5b"
        )
      ),
      
      div(
        class = "table-box",
        style = "flex:1; min-width:300px;",
        
        withSpinner(
          uiOutput(ns("ucsc_button")),
          type = 4,
          color = "#8b1e5b"
        )
      )
    )
  )
}


# =====================
# SERVER
# =====================
geneViewerServer <- function(id, pool, selected_gene){
  
  moduleServer(id, function(input, output, session){
    
    `%||%` <- function(a,b){
      if(is.null(a)) b else a
    }
    
    # =====================
    # SYNC INPUT
    # =====================
    observeEvent(selected_gene(),{
      
      gene <- selected_gene()
      if(
        is.null(gene) ||
        length(gene)==0 ||
        is.na(gene) ||
        trimws(gene)==""
      ){
        return()
      }
      gene <- as.character(gene)[1]
      current_input <- input$gene %||% ""
      
      if(!identical(gene,current_input)){
        updateTextInput(
          session,
          "gene",
          value=gene
        )
      }
    }, ignoreInit=TRUE)
    
    # =====================
    # CURRENT GENE
    # =====================
    current_gene <- reactive({
      gene <- selected_gene()
      
      if(
        !is.null(gene) &&
        length(gene) > 0 &&
        gene != ""
      ){
        return(trimws(gene))
      }
      
      if(
        !is.null(input$gene) &&
        trimws(input$gene) != ""
      ){
        return(trimws(input$gene))
      }
      
      NULL
    })
    
    # =====================
    # VARIANTS
    # =====================
    gene_variants <- reactive({
      
      req(current_gene())
      
      tryCatch(
        {
          get_variants_by_gene(
            pool,
            current_gene()
          )
        },
        error=function(e){
          print(e)
          data.frame()
        }
      )
    })
    
    
    # =====================
    # GENE INFO
    # =====================
    gene_info_filtered <- reactive({
      
      req(current_gene())
      
      tryCatch(
        {
          get_gene_info_by_gene(
            pool,
            current_gene()
          )
        },
        error=function(e){
          print(e)
          data.frame()
        }
      )
    })
    
    
    # =====================
    # SUMMARY
    # =====================
    output$gene_summary <- renderUI({
      
      df <- gene_variants()
      
      if(nrow(df)==0){
        return(tags$p("No gene found"))
      }
      
      gene <- unique(df$`Gene name`)[1]
      gene_id <- unique(df$`Gene ID`)[1]
      n_var <- nrow(df)
      
      tags$div(
        style="display:flex; justify-content:space-between;",
        
        tags$div(
          tags$h3(
            gene,
            style="color:#8b1e5b;"
          ),
          
          tags$p(
            paste("Gene ID:",gene_id),
            style="margin:0;color:#555;"
          )
        ),
        
        tags$div(
          style="text-align:right;",
          
          tags$div(
            "Variants",
            style="font-size:12px;color:#777;"
          ),
          
          tags$div(
            n_var,
            style="font-size:22px;font-weight:700;color:#8b1e5b;"
          )
        )
      )
    })
    
    
    # =====================
    # GENE INFO
    # =====================
    output$gene_info <- renderUI({
      
      df <- gene_info_filtered()
      
      if(nrow(df)==0){
        return(NULL)
      }
      
      row <- df[1,]
      
      tags$div(
        
        tags$h4(
          "Gene information",
          style="color:#8b1e5b;"
        ),
        
        tags$table(
          class="table table-sm gene-info-table",
          
          tags$tr(tags$th("HGNC ID"),tags$td(row$HGNC_ID)),
          tags$tr(tags$th("Biotype"),tags$td(row$BIOTYPE)),
          tags$tr(tags$th("Gene phenotype"),tags$td(row$GENE_PHENO)),
          tags$tr(tags$th("Function"),tags$td(row$Function_description)),
          tags$tr(tags$th("Disease"),tags$td(row$Disease_description)),
          tags$tr(tags$th("HPO ID"),tags$td(row$HPO_id)),
          tags$tr(tags$th("HPO name"),tags$td(row$HPO_name))
        )
      )
    })
    
    
    # =====================
    # EXTERNAL LINKS
    # =====================
    output$external_links <- renderUI({
      
      df_info <- gene_info_filtered()
      df_var <- gene_variants()
      
      if(nrow(df_var)==0){
        return(NULL)
      }
      
      gene <- unique(df_var$`Gene name`)[1]
      
      genecards_url <- paste0(
        "https://www.genecards.org/cgi-bin/carddisp.pl?gene=",
        gene
      )
      
      gtex_url <- paste0(
        "https://gtexportal.org/home/multiGeneQueryPage/",
        gene
      )
      
      omim <- NULL
      
      if(
        !is.null(df_info) &&
        "OMIM_id" %in% colnames(df_info)
      ){
        
        omim <- unique(
          na.omit(
            trimws(df_info$OMIM_id)
          )
        )[1]
      }
      
      valid_omim <- (
        !is.null(omim) &&
          length(omim)>0 &&
          omim!="" &&
          !is.na(omim)
      )
      
      
      tags$div(
        
        tags$h4(
          "External resources",
          style="color:#8b1e5b;"
        ),
        
        tags$div(
          style="display:flex;gap:10px;flex-wrap:wrap;",
          
          tags$a(
            "GeneCards",
            href=genecards_url,
            onclick=paste0(
              "Shiny.setInputValue('",
              session$ns("open_external_url"),
              "',this.href,{priority:'event'});return false;"
            ),
            class="btn-download"
          ),
          
          tags$a(
            "GTEx Portal",
            href=gtex_url,
            onclick=paste0(
              "Shiny.setInputValue('",
              session$ns("open_external_url"),
              "',this.href,{priority:'event'});return false;"
            ),
            class="btn-download"
          ),
          
          if(valid_omim){
            
            tags$a(
              "OMIM",
              href=paste0(
                "https://www.omim.org/entry/",
                omim
              ),
              onclick=paste0(
                "Shiny.setInputValue('",
                session$ns("open_external_url"),
                "',this.href,{priority:'event'});return false;"
              ),
              class="btn-download"
            )
            
          } else {
            
            tags$button(
              "OMIM",
              class="btn-download",
              disabled=TRUE,
              style="
              opacity:0.5;
              cursor:not-allowed;"
            )
          }
        )
      )
    })
    
    
    # =====================
    # UCSC
    # =====================
    output$ucsc_button <- renderUI({
      
      req(current_gene())
      
      df <- gene_variants()
      
      if(nrow(df)==0){
        return(NULL)
      }
      
      chr <- unique(df$CHROM)[1]
      start <- min(df$POS,na.rm=TRUE)
      end <- max(df$POS,na.rm=TRUE)
      
      if(
        is.na(chr) ||
        is.na(start) ||
        is.na(end)
      ){
        return(NULL)
      }
      
      tags$div(
        
        tags$h4(
          "Genome browser",
          style="color:#8b1e5b;"
        ),
        
        tags$a(
          "Open in UCSC Genome Browser",
          
          href=paste0(
            "https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=",
            chr,":",
            start-1000,
            "-",
            end+1000
          ),
          
          onclick=paste0(
            "Shiny.setInputValue('",
            session$ns("open_external_url"),
            "',this.href,{priority:'event'});return false;"
          ),
          
          class="btn-download"
        )
      )
    })
    
    
    # =====================
    # OPEN URL
    # =====================
    observeEvent(input$open_external_url,{
      
      req(input$open_external_url)
      
      utils::browseURL(
        input$open_external_url
      )
      
    })
    
  })
}
