# modules/introduction.R

# =====================
# UI
# =====================
introUI <- function(id){
  
  ns <- NS(id)
  
  div(
    class="content",
    
    h2("📄 | INTRODUCTION"),
    
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
      )
    ),
    
    # =========================
    # MAIN INTRO CARD
    # =========================
    
    div(
      
      class="filter-card intro-main",
      
      h4(
        icon("info-circle"),
        " About OmicsExplorer"
      ),
      
      p(
        strong("OmicsExplorer"),
        " is an interactive application developed for the integration and clinical interpretation of genomic and transcriptomic data, with a particular focus on patients with rare diseases without a definitive diagnosis."
      ),
      
      p(
        "This tool has been developed as part of a Master’s Thesis within the Master’s Degree in Bioinformatics and Computational Biology at the Universidad Autónoma de Madrid, in collaboration with the Instituto de Salud Carlos III, and is framed within the SpainUDP (Spain Undiagnosed Rare Diseases Program)."
      )
    ),
    
    br(),
    
    # =========================
    # TWO CARDS
    # =========================
    
    fluidRow(
      
      column(
        6,
        
        div(
          
          class="filter-card intro-card",
          
          h4(
            icon("dna"),
            " Background"
          ),
          
          p(
            "Advances in Next-Generation Sequencing technologies (NGS), such as Whole Exome Sequencing (WES), Whole Genome Sequencing (WGS), and RNA-seq, have enabled the identification of genetic variants and transcriptomic alterations associated with disease."
          ),
          
          p(
            "However, interpretation of these datasets remains one of the main challenges in clinical genomics."
          )
        )
      ),
      
      column(
        6,
        
        div(
          
          class="filter-card intro-card",
          
          h4(
            icon("chart-line"),
            " OmicsExplorer goal"
          ),
          
          p(
            "OmicsExplorer provides a unified environment for exploring and prioritizing candidate variants through the integration of genomic annotations and transcriptomic evidence."
          ),
          
          p(
            "The platform incorporates filtering strategies, functional annotations and interactive visualizations to support clinical interpretation."
          )
        )
      )
    ),
    
    br(),
    
    div(
      
      class="filter-card intro-highlight",
      
      h4(
        icon("bullseye"),
        " Final objective"
      ),
      
      p(
        "By combining multi-omics data with intuitive analytical tools, OmicsExplorer aims to facilitate data interpretation and support clinical decision-making in the context of precision medicine."
      )
    )
  )
}

# =====================
# SERVER
# =====================
introServer <- function(id){
  moduleServer(id, function(input, output, session){})
}
