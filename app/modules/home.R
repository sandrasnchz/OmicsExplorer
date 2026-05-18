# modules/home.R

# =====================
# UI
# =====================
homeUI <- function(id){
  
  ns <- NS(id)
  
  div(
    class="content",
    
    h2("🏠 | HOME"),
    
    # =========================
    # HERO
    # =========================
    
    div(
      class="home-hero",
      
      h1(
        class="home-title",
        "OmicsExplorer"
      ),
      
      p(
        class="home-subtitle",
        "Multi-omics integration for clinical variant interpretation"
      ),
      
      div(
        class="filter-card home-description",
        tags$span(
          " OmicsExplorer is an interactive Shiny application designed to support the analysis and prioritization of genetic variants through the integration of genomic (WES/WGS) and transcriptomic (RNA-seq) data."
        )
      )
    ),
    
    # =========================
    # CARDS
    # =========================
    
    fluidRow(
      
      column(
        6,
        
        div(
          class="filter-card home-card",
          
          h4(
            icon("flask"),
            " Key Features"
          ),
          
          tags$ul(
            class="home-list",
            
            tags$li("Integration of genomic and transcriptomic data"),
            tags$li("Variant filtering and functional annotation"),
            tags$li("Gene-level exploration and visualization"),
            tags$li("Interactive plots and data export")
          )
        )
      ),
      
      column(
        6,
        
        div(
          class="filter-card home-card",
          
          h4(
            icon("dna"),
            " Use Case"
          ),
          
          p(
            "Designed for researchers and clinical geneticists working on rare diseases, OmicsExplorer facilitates the identification and prioritization of candidate variants through an intuitive and unified analysis environment."
          )
        )
      )
    ),
    
    br(),
    
    div(
      class="footer home-footer",
      "Developed as part of a Master's Thesis in Bioinformatics and Computational Biology (UAM)"
    )
    
  )
}

# =====================
# SERVER
# =====================
homeServer <- function(id){
  moduleServer(id, function(input, output, session){})
}
