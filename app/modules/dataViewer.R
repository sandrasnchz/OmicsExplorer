# modules/data_viewer.R

# =====================
# UI
# =====================
dataViewerUI <- function(id){
  
  ns <- NS(id)
  
  div(
    class="content",
    
    h2("🔎 | DATA VIEWER"),
    
    div(
      class="filter-box",
      
      textInput(
        ns("gene"),
        "Filter by gene:",
        placeholder="e.g. DPM1"
      )
    ),
    
    tabsetPanel(
      
      id=ns("main_tabs"),
      
      # ==========================================================
      # VARIANTS
      # ==========================================================
      tabPanel(
        
        "Variants (WES + WGS)",
        
        div(
          
          class="filter-box",
          
          # =====================
          # QUALITY FILTERS
          # =====================
          
          div(
            
            class="filter-card",
            
            h4(icon("sliders-h"), " Quality filters"),
            
            fluidRow(
              
              column(
                4,
                
                sliderInput(
                  ns("af"),
                  "Max AF",
                  min=0,
                  max=1,
                  value=0.05,
                  step=0.01
                )
              ),
              
              column(
                4,
                
                div(
                  
                  numericInput(
                    ns("dp"),
                    "Min DP",
                    value=90,
                    min=0,
                    step=1,
                    width="80px"
                  ),
                  
                  actionButton(
                    ns("apply_dp"),
                    "Apply",
                    class="btn-download"
                  )
                )
              ),
              
              column(
                4,
                
                sliderInput(
                  ns("gq"),
                  "GQ range",
                  min=0,
                  max=100,
                  value=c(90,100),
                  step=1
                )
              )
              
            )
          ),
          
          # =====================
          # VARIANT FILTERS
          # =====================
          
          div(
            
            class="filter-card",
            
            h4(icon("dna"), " Variant filters"),
            
            fluidRow(
              
              column(
                6,
                
                prettyCheckboxGroup(
                  inputId=ns("impact"),
                  label="Impact",
                  choices=c(
                    "HIGH",
                    "MODERATE",
                    "LOW",
                    "MODIFIER"
                  ),
                  selected=c(
                    "HIGH",
                    "MODERATE",
                    "LOW",
                    "MODIFIER"
                  ),
                  inline=TRUE,
                  shape="square",
                  icon=icon("check")
                )
              ),
              
              column(
                3,
                
                prettyCheckboxGroup(
                  inputId=ns("source"),
                  label="Source",
                  choices=c(
                    "WES",
                    "WGS",
                    "BOTH"
                  ),
                  selected=c(
                    "WES",
                    "WGS",
                    "BOTH"
                  ),
                  inline=FALSE,
                  shape="square",
                  icon=icon("check")
                )
              ),
              
              column(
                3,
                
                prettyCheckboxGroup(
                  inputId=ns("variant_class"),
                  label="Variant class",
                  choices=c(
                    "SNV",
                    "insertion",
                    "deletion"
                  ),
                  selected=c(
                    "SNV",
                    "insertion",
                    "deletion"
                  ),
                  inline=FALSE,
                  shape="square",
                  icon=icon("check")
                )
              )
              
            )
          ),
          
          # =====================
          # INHERITANCE
          # =====================
          
          div(
            
            class="filter-card",
            
            h4(icon("project-diagram"), " Inheritance & gene filters"),
            
            fluidRow(
              
              column(
                
                5,
                
                prettyCheckboxGroup(
                  inputId=ns("inheritance_source"),
                  label="Inheritance source",
                  choices=c(
                    "AR",
                    "AD",
                    "XR",
                    "XD",
                    "MT"
                  ),
                  selected=c(
                    "AR",
                    "AD",
                    "XR",
                    "XD",
                    "MT"
                  ),
                  inline=TRUE,
                  shape="square",
                  icon=icon("check")
                ),
                
                br(),
                
                prettyCheckboxGroup(
                  inputId=ns("de_novo_source"),
                  label=" ",
                  choices=c(
                    "AD - de_novo"="AD",
                    "XD - de_novo"="XD"
                  ),
                  selected=c(
                    "AD",
                    "XD"
                  ),
                  inline=TRUE,
                  shape="square",
                  icon=icon("check")
                )
                
              ),
              
              column(
                4,
                
                br(),
                
                div(
                  class="omics-switch",
                  
                  prettySwitch(
                    inputId=ns("hide_intergenic"),
                    label="Hide intergenic variants",
                    value=TRUE,
                    fill=TRUE
                  )
                ),
                
                div(
                  class="omics-switch",
                  
                  prettySwitch(
                    inputId=ns("only_omim"),
                    label="Only OMIM genes",
                    value=FALSE,
                    fill=TRUE
                  )
                ),
                
                div(
                  class="omics-switch",
                  
                  prettySwitch(
                    inputId=ns("only_canonical"),
                    label="Only canonical transcripts",
                    value=TRUE,
                    fill=TRUE
                  )
                )
              )
              
            )
          )
          
        ),
        
        br(),
        
        withSpinner(
          DTOutput(ns("variants")),
          type=4,
          color="#8b1e5b"
        ),
        
        div(
          class="download-box",
          
          radioButtons(
            ns("variants_download_format"),
            "Download format:",
            choices=c(
              "TSV"="tsv",
              "TXT"="txt"
            ),
            selected="tsv",
            inline=TRUE
          ),
          
          downloadButton(
            ns("download_variants"),
            "Download variants",
            class="omics-download-btn"
          )
        )
      ),
      
      # ==========================================================
      # RNA
      # ==========================================================
      
      tabPanel(
        
        "RNA Data",
        
        div(
          
          class="filter-box",
          
          div(
            class="filter-card",
            
            h4(icon("dna"), " RNA filters"),
            
            fluidRow(
              
              column(
                6,
                
                prettyCheckboxGroup(
                  ns("drop_filters"),
                  "DROP filters",
                  
                  choices=c(
                    "Aberrant Expression"="expr",
                    "Aberrant Splicing"="splicing",
                    "MAE"="mae"
                  ),
                  
                  inline=TRUE,
                  shape="curve",
                  icon=icon("check")
                )
              ),
              
              column(
                6,
                
                prettyCheckboxGroup(
                  ns("tpm_filters"),
                  "TPM outliers",
                  
                  choices=c(
                    "Overexpressed"="over",
                    "Infraexpressed"="under"
                  ),
                  
                  inline=TRUE,
                  shape="curve",
                  icon=icon("check")
                )
              )
            )
          )
        ),
        
        withSpinner(
          DTOutput(ns("rna")),
          type=4,
          color="#8b1e5b"
        ),
        
        div(
          
          class="download-box",
          
          radioButtons(
            ns("rna_download_format"),
            "Download format:",
            choices=c(
              "TSV"="tsv",
              "TXT"="txt"
            ),
            selected="tsv",
            inline=TRUE
          ),
          
          downloadButton(
            ns("download_rna"),
            "Download RNA table",
            class="omics-download-btn"
          )
        )
        
      )
      
    )
  )
}

# =====================
# SERVER
# =====================
dataViewerServer <- function(id, pool, selected_gene){
  moduleServer(id, function(input, output, session){
    
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    all_impacts <- c("HIGH","MODERATE","LOW","MODIFIER")
    all_sources <- c("WES","WGS","BOTH")
    all_inheritance_sources <- c("AR","AD","XR","XD","MT")
    all_de_novo_sources <- c("AD","XD")
    all_variant_classes <- c("SNV","insertion","deletion")
    
    navigating <- reactiveVal(FALSE)
    
    dp_filter <- reactiveVal(90)
    
    observeEvent(input$apply_dp,{
      
      dp_filter(
        input$dp
      )
      
    })
    
    # =====================
    # RESET CONTROLADO
    # =====================
    observeEvent(input$main_tabs,{
      
      if(!navigating()){
        
        selected_gene(NULL)
        
        updateTextInput(
          session,
          "gene",
          value=""
        )
      }
      navigating(FALSE)
    })
    
    # =====================
    # SINCRONIZAR INPUT GENE
    # =====================
    observeEvent(selected_gene(),{
      
      if(
        !is.null(selected_gene()) &&
        selected_gene() != ""
      ){
        
        updateTextInput(
          session,
          "gene",
          value = selected_gene()
        )
      }
    }, ignoreInit=TRUE)
    
    
    # =====================
    # FUNCIÓN PARA OBTENER FLAGS DE DROP (PARA RESUMEN EN RNA)
    # =====================
    get_drop_flags <- function(pool){
      
      expr <- get_drop_expr(pool)
      spl  <- get_drop_splicing(pool)
      mae  <- get_drop_mae(pool)
      
      # =====================
      # EXPRESSION
      # =====================
      
      if(nrow(expr) > 0){
        
        gene_col_expr <- if("gene_name" %in% colnames(expr)) "gene_name" else "hgncSymbol"
        
        expr <- expr %>%
          mutate(drop_expr = TRUE) %>%
          distinct(.data[[gene_col_expr]], .keep_all = TRUE) %>%
          select(
            gene = all_of(gene_col_expr),
            drop_expr
          )
        
      } else {
        
        expr <- data.frame(
          gene = character(),
          drop_expr = logical()
        )
      }
      
      # =====================
      # SPLICING
      # =====================
      
      if(nrow(spl) > 0){
        
        gene_col_spl <- if("gene_name" %in% colnames(spl)) "gene_name" else "hgncSymbol"
        
        spl <- spl %>%
          mutate(drop_splicing = TRUE) %>%
          distinct(.data[[gene_col_spl]], .keep_all = TRUE) %>%
          select(
            gene = all_of(gene_col_spl),
            drop_splicing
          )
        
      } else {
        
        spl <- data.frame(
          gene = character(),
          drop_splicing = logical()
        )
      }
      
      # =====================
      # MAE
      # =====================
      
      if(nrow(mae) > 0){
        
        gene_col_mae <- if("gene_name" %in% colnames(mae)) "gene_name" else "hgncSymbol"
        
        mae <- mae %>%
          mutate(drop_mae = TRUE) %>%
          distinct(.data[[gene_col_mae]], .keep_all = TRUE) %>%
          select(
            gene = all_of(gene_col_mae),
            drop_mae
          )
        
      } else {
        
        mae <- data.frame(
          gene = character(),
          drop_mae = logical()
        )
      }
      
      # =====================
      # JOIN
      # =====================
      
      df <- full_join(expr, spl, by = "gene") %>%
        full_join(mae, by = "gene") %>%
        mutate(
          across(starts_with("drop"), ~replace_na(., FALSE))
        )
      
      return(df)
    }
    
    get_value <- function(row, field){
      if(!field %in% names(row)) return("")
      value <- row[[field]][1]
      if(is.null(value) || is.na(value) || identical(value, "")) return("-")
      htmltools::htmlEscape(as.character(value))
    }
    
    variant_detail_html <- function(row){
      paste0(
        "<div class='variant-detail-tabs-scroll'>",
        "<ul class='nav nav-tabs variant-detail-tabs'>",
        "<li class='nav-item'><a class='nav-link active freq-tab' data-target='#geno-detail'>Trio genotypes</a></li>",
        "<li class='nav-item'><a class='nav-link freq-tab' data-target='#pop-detail'>Population frequencies</a></li>",
        "<li class='nav-item'><a class='nav-link freq-tab' data-target='#pred-detail'>Predictors</a></li>",
        "<li class='nav-item'><a class='nav-link freq-tab' data-target='#trans-detail'>Transcript / Protein</a></li>",
        "<li class='nav-item'><a class='nav-link freq-tab' data-target='#annot-detail'>Transcript annotations</a></li>",
        "<li class='nav-item'><a class='nav-link freq-tab' data-target='#hgvs-detail'>HGVS info</a></li>",
        "<li class='nav-item'><a class='nav-link freq-tab' data-target='#gene-detail'>Gene information</a></li>",
        "<li class='nav-item'><a class='nav-link freq-tab' data-target='#motif-detail'>Regulation and motifs</a></li>",
        "</ul>",
        "</div>",
        
        "<div class='tab-content' style='margin-top:10px'>",
        
        "<div class='tab-pane active' id='geno-detail'>",
        "<table class='table table-sm table-bordered genotype-table'>",
        "<tr><th></th><th>GT</th><th>DP</th><th>AD</th><th>GQ</th></tr>",
        "<tr><td><b>P1</b></td><td>", get_value(row, "PARENT1_GT"), "</td><td>", get_value(row, "PARENT1_DP"), "</td><td>", get_value(row, "PARENT1_AD"), "</td><td>", get_value(row, "PARENT1_GQ"), "</td></tr>",
        "<tr><td><b>P2</b></td><td>", get_value(row, "PARENT2_GT"), "</td><td>", get_value(row, "PARENT2_DP"), "</td><td>", get_value(row, "PARENT2_AD"), "</td><td>", get_value(row, "PARENT2_GQ"), "</td></tr>",
        "<tr><td><b>Child</b></td><td>", get_value(row, "CHILD_GT"), "</td><td>", get_value(row, "CHILD_DP"), "</td><td>", get_value(row, "CHILD_AD"), "</td><td>", get_value(row, "CHILD_GQ"), "</td></tr>",
        "</table>",
        "</div>",
        
        "<div class='tab-pane' id='pop-detail'>",
        "<ul class='nav nav-tabs'>",
        "<li class='nav-item'><a class='nav-link active subtab' data-target='#general-detail'>General population</a></li>",
        "<li class='nav-item'><a class='nav-link subtab' data-target='#gnom-detail'>gnomAD frequencies</a></li>",
        "</ul>",
        "<div class='tab-content' style='margin-top:10px'>",
        
        "<div class='tab-pane active' id='general-detail'>",
        "AF: ", get_value(row, "AF"), "<br>",
        "AFR: ", get_value(row, "AFR_AF"), "<br>",
        "AMR: ", get_value(row, "AMR_AF"), "<br>",
        "EAS: ", get_value(row, "EAS_AF"), "<br>",
        "EUR: ", get_value(row, "EUR_AF"), "<br>",
        "SAS: ", get_value(row, "SAS_AF"), "<br>",
        "AA: ", get_value(row, "AA_AF"), "<br>",
        "EA: ", get_value(row, "EA_AF"),
        "</div>",
        

        "<div class='tab-pane' id='gnom-detail'>",
        "AF: ", get_value(row, "gnomAD_AF"), "<br>",
        "AFR: ", get_value(row, "gnomAD_AFR_AF"), "<br>",
        "AMR: ", get_value(row, "gnomAD_AMR_AF"), "<br>",
        "ASJ: ", get_value(row, "gnomAD_ASJ_AF"), "<br>",
        "EAS: ", get_value(row, "gnomAD_EAS_AF"), "<br>",
        "FIN: ", get_value(row, "gnomAD_FIN_AF"), "<br>",
        "NFE: ", get_value(row, "gnomAD_NFE_AF"), "<br>",
        "OTH: ", get_value(row, "gnomAD_OTH_AF"),
        "</div>",
        "</div>",
        "</div>",
        
        "<div class='tab-pane' id='pred-detail'>",
        "<ul class='nav nav-tabs'>",
        "<li class='nav-item'><a class='nav-link active subtab' data-target='#sift-detail'>SIFT</a></li>",
        "<li class='nav-item'><a class='nav-link subtab' data-target='#poly-detail'>PolyPhen</a></li>",
        "<li class='nav-item'><a class='nav-link subtab' data-target='#mut-detail'>Mutation</a></li>",
        "<li class='nav-item'><a class='nav-link subtab' data-target='#meta-detail'>Meta</a></li>",
        "<li class='nav-item'><a class='nav-link subtab' data-target='#cadd-detail'>CADD</a></li>",
        "</ul>",
        "<div class='tab-content' style='margin-top:10px'>",
        "<div class='tab-pane active' id='sift-detail'>Score: ", get_value(row, "SIFT_score"), "<br>Pred: ", get_value(row, "SIFT_pred"), "</div>",
        "<div class='tab-pane' id='poly-detail'>HDIV: ", get_value(row, "Polyphen2_HDIV_score"), " (", get_value(row, "Polyphen2_HDIV_pred"), ")<br>",
        "HVAR: ", get_value(row, "Polyphen2_HVAR_score"), " (", get_value(row, "Polyphen2_HVAR_pred"), ")</div>",
        "<div class='tab-pane' id='mut-detail'>Taster: ", get_value(row, "MutationTaster_score"), " (", get_value(row, "MutationTaster_pred"), ")<br>",
        "Assessor: ", get_value(row, "MutationAssessor_score"), " (", get_value(row, "MutationAssessor_pred"), ")</div>",
        "<div class='tab-pane' id='meta-detail'>MetaSVM: ", get_value(row, "MetaSVM_score"), " (", get_value(row, "MetaSVM_pred"), ")<br>",
        "MetaLR: ", get_value(row, "MetaLR_score"), " (", get_value(row, "MetaLR_pred"), ")</div>",
        "<div class='tab-pane' id='cadd-detail'>Raw: ", get_value(row, "CADD_raw"), "<br>Phred: ", get_value(row, "CADD_phred"), "</div>",
        "</div>",
        "</div>",
        
        "<div class='tab-pane' id='trans-detail'>",
        "<b>Transcript position</b><br><br>",
        "cDNA position: ", get_value(row, "cDNA_position"), "<br>",
        "CDS position: ", get_value(row, "CDS_position"), "<br>",
        "Protein position: ", get_value(row, "Protein_position"), "<br><br>",
        "<b>Protein change</b><br><br>",
        "Amino acids: ", get_value(row, "Amino_acids"), "<br>",
        "Codons: ", get_value(row, "Codons"), "<br><br>",
        "<b>Genomic location</b><br><br>",
        "Exon: ", get_value(row, "EXON"), "<br>",
        "Intron: ", get_value(row, "INTRON"),
        "</div>",
        
        "<div class='tab-pane' id='annot-detail'>",
        "<b>Transcript quality</b><br><br>",
        "Canonical: ", get_value(row, "CANONICAL"), "<br>",
        "MANE Select: ", get_value(row, "MANE_sel"), "<br>",
        "MANE Plus Clinical: ", get_value(row, "MANE_plus"), "<br>",
        "TSL: ", get_value(row, "TSL"), "<br>",
        "APPRIS: ", get_value(row, "APPRIS"), "<br><br>",
        "<b>Transcript identifiers</b><br><br>",
        "CCDS: ", get_value(row, "CCDS"), "<br>",
        "ENSP: ", get_value(row, "ENSP"), "<br>",
        "SwissProt: ", get_value(row, "SWISSPROT"), "<br>",
        "TrEMBL: ", get_value(row, "TREMBL"), "<br>",
        "UniParc: ", get_value(row, "UNIPARC"),
        "</div>",
        
        "<div class='tab-pane' id='hgvs-detail'>",
        "<b>VEP HGVS</b><br><br>",
        "HGVSc: ", get_value(row, "HGVSc"), "<br>",
        "HGVSp: ", get_value(row, "HGVSp"), "<br>",
        "HGVS offset: ", get_value(row, "HGVS_OFFSET"), "<br><br>",
        "<b>snpEff HGVS</b><br><br>",
        "HGVSc snpEff: ", get_value(row, "HGVSc_snpEff"), "<br>",
        "HGVSp snpEff: ", get_value(row, "HGVSp_snpEff"),
        "</div>",
        
        "<div class='tab-pane' id='gene-detail'>",
        "<b>Gene information</b><br><br>",
        " ", get_value(row, "Function_description"),"<br><br>",
        "OMIM ID: ", get_value(row, "OMIM_id"), "</div>",
        
        "<div class='tab-pane' id='motif-detail'>",
        "<b>Regulation y motifs</b><br><br>",
        "Motif name: ", get_value(row, "MOTIF_NAME"), "<br>",
        "Motif position: ", get_value(row, "MOTIF_POS"), "<br>",
        "High information position: ", get_value(row, "HIGH_INF_POS"), "<br>",
        "Motif score change: ", get_value(row, "MOTIF_SCORE_CHANGE"), "<br>",
        "Transcription factors: ", get_value(row, "TRANSCRIPTION_FACTORS"),
        "</div>",
        
        "</div>",
        
        "<hr>",
        "<div class='explore-box'>",
        "<div class='explore-title'>Explore</div>",
        "<div class='explore-desc'>Navigate to gene viewer or RNA data table</div>",
        "<div style='display:flex; gap:10px; margin-top:8px;'>",
        "<button class='go-gene' data-gene='", get_value(row, "SYMBOL"), "'>View gene info</button>",
        "<button class='go-rna' data-gene='", get_value(row, "SYMBOL"), "'>Go to RNA Data table</button>",
        "</div>",
        "</div>"
      )
    }
    
    variants_files <- reactivePoll(
      intervalMillis = 1000,
      session = session,
      checkFunc = function(){
        files <- Sys.glob("../data/variants/*.parquet")
        if(length(files) == 0) return("")
        
        info <- file.info(files)
        paste(files, info$mtime, info$size, collapse = "|")
      },
      valueFunc = function(){
        Sys.glob("../data/variants/*.parquet")
      }
    )
    
    variants_all <- reactive({
      
      req(pool)
      
      # fuerza invalidación cuando cambian archivos
      files <- variants_files()
      
      validate(
        need(length(files) > 0,
             "No variant files loaded")
      )
      
      # limpiar posibles resultados anteriores
      gc()
      
      isolate({
        
        get_variants_with_inheritance(pool) %>%
          mutate(
            .variant_row_key = row_number()
          )
        
      })
    })
    
    format_inheritance_display <- function(x){
      x <- as.character(x)
      p1 <- sub(".*(P1:[^|\\)]+).*", "\\1", x)
      p2 <- sub(".*(P2:[^\\)]+).*", "\\1", x)
      
      ifelse(
        grepl("P1:", p1) & grepl("P2:", p2),
        paste0(trimws(p1), "<br>", trimws(p2)),
        x
      )
    }
    
    filtered_variants_data <- reactive({
      df <- variants_all()
      
      if(nrow(df) == 0){
        return(df)
      }
      
      if(!is.null(selected_gene()) && selected_gene() != ""){
        df <- df %>% filter(toupper(SYMBOL) == toupper(selected_gene()))
      }
      
      if(nzchar(input$gene)){
        df <- df %>% filter(toupper(SYMBOL) == toupper(input$gene))
      }
      
      if(!is.null(input$af) && input$af < 1){
        df <- df %>%
          mutate(MAX_AF_FILTER = suppressWarnings(as.numeric(MAX_AF))) %>%
          filter(is.na(MAX_AF_FILTER) | MAX_AF_FILTER <= input$af) %>%
          select(-MAX_AF_FILTER)
      }
      
      if(!is.null(dp_filter())){
        
        df <- df %>%
          mutate(
            CHILD_DP_FILTER=
              suppressWarnings(
                as.numeric(CHILD_DP)
              )
          ) %>%
          filter(
            is.na(CHILD_DP_FILTER) |
              CHILD_DP_FILTER >= dp_filter()
          ) %>%
          select(-CHILD_DP_FILTER)
      }
      
      if(!is.null(input$gq) && !identical(input$gq, c(0, 100))){
        df <- df %>%
          mutate(CHILD_GQ_FILTER = suppressWarnings(as.numeric(CHILD_GQ))) %>%
          filter(
            is.na(CHILD_GQ_FILTER) |
              (CHILD_GQ_FILTER >= input$gq[1] &
                 CHILD_GQ_FILTER <= input$gq[2])
          ) %>%
          select(-CHILD_GQ_FILTER)
      }
      
      if(length(input$impact) > 0 && !setequal(input$impact, all_impacts)){
        df <- df %>% filter(IMPACT %in% input$impact)
      }
      
      if(length(input$source) > 0 && !setequal(input$source, all_sources)){
        df <- df %>% filter(source %in% input$source)
      }
      
      # =====================
      # INHERITANCE FILTER
      # =====================
      
      selected_main <- input$inheritance_source
      selected_denovo <- input$de_novo_source
      
      df <- df %>%
        filter(
          (inheritance_source %in% selected_main &
             inheritance_type != "de_novo") |
            
            (inheritance_source %in% selected_denovo &
               inheritance_type == "de_novo") |
            
            (inheritance_source %in% intersect(c("AD","XD"), selected_main) &
               inheritance_type == "de_novo")
        )
      
      if(nrow(df)==0) df <- df[0,]
      
      if(input$hide_intergenic){
        df <- df %>%
          filter(
            !grepl("intergenic_variant", Consequence, ignore.case = TRUE)
          )
      }
      
      if(input$only_omim){
        df <- df %>%
          filter(
            !is.na(OMIM_id) &
              OMIM_id != "" &
              OMIM_id != "." &
              OMIM_id != "-"
          )
      }
      
      if(input$only_canonical && "CANONICAL" %in% colnames(df)){
        df <- df %>%
          filter(toupper(CANONICAL) == "YES")
      }
      
      if(length(input$variant_class) > 0 && !setequal(input$variant_class, all_variant_classes)){
        df <- df %>%
          filter(VARIANT_CLASS %in% input$variant_class)
      }
      
      df
    })
    
    write_download_table <- function(df, file){
      if(".variant_row_key" %in% colnames(df)){
        df <- df %>% select(-.variant_row_key)
      }
      
      write.table(
        df,
        file = file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
      )
    }
    
    format_impact_badge <- function(x){
      impact <- as.character(x)
      class <- paste0("impact-badge impact-", tolower(impact))
      ifelse(
        is.na(impact) | impact == "",
        "",
        paste0("<span class='", class, "'>", htmltools::htmlEscape(impact), "</span>")
      )
    }
    
    
    # =====================
    # VARIANTS
    # =====================
    output$variants <- renderDT({
      
      req(pool)
      
      tryCatch({
        
        df <- variants_all()
        
        if(nrow(df) == 0){
          return(datatable(data.frame(Message="No variants loaded")))
        }
        
        # ===== FILTROS =====
        if(!is.null(selected_gene()) && selected_gene() != ""){
          df <- df %>% filter(toupper(SYMBOL) == toupper(selected_gene()))
        }
        
        if(nzchar(input$gene)){
          df <- df %>% filter(toupper(SYMBOL) == toupper(input$gene))
        }
        
        if(!is.null(input$af) && input$af < 1){
          df <- df %>%
            mutate(MAX_AF_FILTER = suppressWarnings(as.numeric(MAX_AF))) %>%
            filter(is.na(MAX_AF_FILTER) | MAX_AF_FILTER <= input$af) %>%
            select(-MAX_AF_FILTER)
        }
        
        if(!is.null(dp_filter())){
          
          df <- df %>%
            mutate(
              CHILD_DP_FILTER=
                suppressWarnings(
                  as.numeric(CHILD_DP)
                )
            ) %>%
            filter(
              is.na(CHILD_DP_FILTER) |
                CHILD_DP_FILTER >= dp_filter()
            ) %>%
            select(-CHILD_DP_FILTER)
        }
        
        if(!is.null(input$gq) && !identical(input$gq, c(0, 100))){
          df <- df %>%
            mutate(CHILD_GQ_FILTER = suppressWarnings(as.numeric(CHILD_GQ))) %>%
            filter(
              is.na(CHILD_GQ_FILTER) |
                (CHILD_GQ_FILTER >= input$gq[1] &
                   CHILD_GQ_FILTER <= input$gq[2])
            ) %>%
            select(-CHILD_GQ_FILTER)
        }
        
        if(length(input$impact) > 0 && !setequal(input$impact, all_impacts)){
          df <- df %>% filter(IMPACT %in% input$impact)
        }
        
        if(length(input$source) > 0 && !setequal(input$source, all_sources)){
          df <- df %>% filter(source %in% input$source)
        }
        
        # =====================
        # INHERITANCE FILTER
        # =====================
        
        selected_main <- input$inheritance_source
        selected_denovo <- input$de_novo_source
        
        df <- df %>%
          filter(
            (inheritance_source %in% selected_main &
               inheritance_type != "de_novo") |
              
              (inheritance_source %in% selected_denovo &
                 inheritance_type == "de_novo") |
              
              (inheritance_source %in% intersect(c("AD","XD"), selected_main) &
                 inheritance_type == "de_novo")
          )
        
        if(nrow(df)==0) df <- df[0,]
        
        if(input$hide_intergenic){
          df <- df %>%
            filter(
              !grepl("intergenic_variant", Consequence, ignore.case = TRUE)
            )
        }
        
        if(input$only_omim){
          df <- df %>%
            filter(
              !is.na(OMIM_id) &
                OMIM_id != "" &
                OMIM_id != "." &
                OMIM_id != "-"
            )
        }
        
        if(input$only_canonical && "CANONICAL" %in% colnames(df)){
          df <- df %>%
            filter(toupper(CANONICAL) == "YES")
        }
        
        if(length(input$variant_class) > 0 && !setequal(input$variant_class, all_variant_classes)){
          df <- df %>%
            filter(VARIANT_CLASS %in% input$variant_class)
        }
        
        if(nrow(df) == 0){
          return(datatable(data.frame(Message="No variants match filters")))
        }
        
        # ===== COLUMNAS =====
        cols_to_show <- c(
          ".variant_row_key",
          "ID","CHROM","POS","REF","ALT","FILTER",
          "CHILD_GT_N","CHILD_DP","CHILD_AD","CHILD_GQ",
          "SYMBOL","Gene","ENSP",
          "Consequence","IMPACT","MAX_AF","VARIANT_CLASS",
          "clinvar_trait",
          "clinvar_id",
          "SIFT_pred", "Polyphen2_HVAR_pred",
          "source","inheritance_source","inheritance_type","inheritance"
        )
        
        cols_to_show <- intersect(cols_to_show, colnames(df))
        df <- df[, cols_to_show]
        
        # ===== RENOMBRE =====
        df <- df %>%
          rename(
            GT = CHILD_GT_N,
            DP = CHILD_DP,
            AD = CHILD_AD,
            GQ = CHILD_GQ,
            gene_name = SYMBOL,
            gene_id   = Gene,
            transcript_id = ENSP,
            consequence = Consequence,
            impact      = IMPACT,
            max_AF      = MAX_AF,
            variant_class = VARIANT_CLASS
          ) %>%
          mutate(
            impact_sort = factor(
              impact,
              levels = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
              ordered = TRUE
            ),
            inheritance = format_inheritance_display(inheritance)
          )
        
        df <- df %>%
          arrange(impact_sort) %>%
          mutate(impact = format_impact_badge(impact)) %>%
          select(-impact_sort)
        
        showNotification(
          "By default, variants are ordered by impact severity: HIGH → MODERATE → LOW → MODIFIER.",
          type = "default",
          duration = 5
        )
        
        # =====================
        # DATATABLE
        # =====================
        datatable(
          df,
          rownames = FALSE,
          escape = FALSE,
          selection = "none",
          extensions = 'FixedHeader',
          options = list(
            scrollX = TRUE,
            scrollCollapse = TRUE,
            autoWidth = TRUE,
            pageLength = 10,
            fixedHeader = TRUE,
            initComplete = JS("function(settings, json) { this.api().columns.adjust(); }"),
            drawCallback = JS("function(settings) { this.api().columns.adjust(); }"),
            columnDefs = list(
              list(
                targets = which(colnames(df) == ".variant_row_key") - 1,
                visible = FALSE,
                searchable = FALSE
              ),
              list(
                targets = which(colnames(df) == "inheritance") - 1,
                width = "100px",
                className = "inheritance-col"
              )
            )
          ),
          
          callback = JS(paste0("

var format = function(rowData){

  // ===== convertir a objeto =====
  var data = {};
  table.columns().every(function(i){
    data[table.column(i).header().innerText] = rowData[i];
  });

  var gene = data['gene_name'];
  var uid = Math.random().toString(36).substring(2,9);

  return '<div class=\"variant-detail\" id=\"variant-detail-body-'+uid+'\" data-uid=\"'+uid+'\" style=\"padding:10px\">Loading variant details...</div>';

  return '<div class=\"variant-detail\" style=\"padding:10px\">' +

    '<div class=\"variant-detail-tabs-scroll\">' +
    '<ul class=\"nav nav-tabs variant-detail-tabs\">' +
      '<li class=\"nav-item\"><a class=\"nav-link active freq-tab\" data-target=\"#geno-'+uid+'\">Trio genotypes</a></li>' +
      '<li class=\"nav-item\"><a class=\"nav-link freq-tab\" data-target=\"#pop-'+uid+'\">Population frequencies</a></li>' +
      '<li class=\"nav-item\"><a class=\"nav-link freq-tab\" data-target=\"#gnom-'+uid+'\">gnomAD frequencies</a></li>' +
      '<li class=\"nav-item\"><a class=\"nav-link freq-tab\" data-target=\"#pred-'+uid+'\">Predictors</a></li>' +
      '<li class=\"nav-item\"><a class=\"nav-link freq-tab\" data-target=\"#trans-'+uid+'\">Transcript / Protein </a></li>' +
      '<li class=\"nav-item\"><a class=\"nav-link freq-tab\" data-target=\"#annot-'+uid+'\">Transcript annotations</a></li>' +
      '<li class=\"nav-item\"><a class=\"nav-link freq-tab\" data-target=\"#hgvs-'+uid+'\">HGVS info</a></li>' +
      '<li class=\"nav-item\"><a class=\"nav-link freq-tab\" data-target=\"#motif-'+uid+'\">Regulación y motifs</a></li>' +
    '</ul>' +
    '</div>' +

    '<div class=\"tab-content\" style=\"margin-top:10px\">' +
    
    '<div class=\"tab-pane active\" id=\"geno-'+uid+'\">' +

  '<table class=\"table table-sm table-bordered genotype-table\">' +

    '<tr>' +
      '<th></th>' +
      '<th>GT</th>' +
      '<th>DP</th>' +
      '<th>AD</th>' +
      '<th>GQ</th>' +
    '</tr>' +

    '<tr>' +
      '<td><b>P1</b></td>' +
      '<td>'+data['PARENT1_GT']+'</td>' +
      '<td>'+data['PARENT1_DP']+'</td>' +
      '<td>'+data['PARENT1_AD']+'</td>' +
      '<td>'+data['PARENT1_GQ']+'</td>' +
    '</tr>' +

    '<tr>' +
      '<td><b>P2</b></td>' +
      '<td>'+data['PARENT2_GT']+'</td>' +
      '<td>'+data['PARENT2_DP']+'</td>' +
      '<td>'+data['PARENT2_AD']+'</td>' +
      '<td>'+data['PARENT2_GQ']+'</td>' +
    '</tr>' +

    '<tr>' +
      '<td><b>Child</b></td>' +
      '<td>'+data['CHILD_GT']+'</td>' +
      '<td>'+data['DP']+'</td>' +
      '<td>'+data['AD']+'</td>' +
      '<td>'+data['GQ']+'</td>' +
    '</tr>' +

   '</table>' +

    '</div>' +
    

      '<div class=\"tab-pane\" id=\"pop-'+uid+'\">' +
        'AF: '+data['AF']+'<br>' +
        'AFR: '+data['AFR_AF']+'<br>' +
        'AMR: '+data['AMR_AF']+'<br>' +
        'EAS: '+data['EAS_AF']+'<br>' +
        'EUR: '+data['EUR_AF']+'<br>' +
        'SAS: '+data['SAS_AF']+'<br>' +
        'AA: '+data['AA_AF']+'<br>' +
        'EA: '+data['EA_AF'] +
      '</div>' +

      '<div class=\"tab-pane\" id=\"gnom-'+uid+'\">' +
        'AF: '+data['gnomAD_AF']+'<br>' +
        'AFR: '+data['gnomAD_AFR_AF']+'<br>' +
        'AMR: '+data['gnomAD_AMR_AF']+'<br>' +
        'ASJ: '+data['gnomAD_ASJ_AF']+'<br>' +
        'EAS: '+data['gnomAD_EAS_AF']+'<br>' +
        'FIN: '+data['gnomAD_FIN_AF']+'<br>' +
        'NFE: '+data['gnomAD_NFE_AF']+'<br>' +
        'OTH: '+data['gnomAD_OTH_AF'] +
      '</div>' +

      '<div class=\"tab-pane\" id=\"pred-'+uid+'\">' +

        '<ul class=\"nav nav-tabs\">' +
          '<li class=\"nav-item\"><a class=\"nav-link active subtab\" data-target=\"#sift-'+uid+'\">SIFT</a></li>' +
          '<li class=\"nav-item\"><a class=\"nav-link subtab\" data-target=\"#poly-'+uid+'\">PolyPhen</a></li>' +
          '<li class=\"nav-item\"><a class=\"nav-link subtab\" data-target=\"#mut-'+uid+'\">Mutation</a></li>' +
          '<li class=\"nav-item\"><a class=\"nav-link subtab\" data-target=\"#meta-'+uid+'\">Meta</a></li>' +
          '<li class=\"nav-item\"><a class=\"nav-link subtab\" data-target=\"#cadd-'+uid+'\">CADD</a></li>' +
        '</ul>' +

        '<div class=\"tab-content\" style=\"margin-top:10px\">' +

          '<div class=\"tab-pane active\" id=\"sift-'+uid+'\">Score: '+data['SIFT_score']+'<br>Pred: '+data['SIFT_pred']+'</div>' +

          '<div class=\"tab-pane\" id=\"poly-'+uid+'\">HDIV: '+data['Polyphen2_HDIV_score']+' ('+data['Polyphen2_HDIV_pred']+')<br>HVAR: '+data['Polyphen2_HVAR_score']+' ('+data['Polyphen2_HVAR_pred']+')</div>' +

          '<div class=\"tab-pane\" id=\"mut-'+uid+'\">Taster: '+data['MutationTaster_score']+' ('+data['MutationTaster_pred']+')<br>Assessor: '+data['MutationAssessor_score']+' ('+data['MutationAssessor_pred']+')</div>' +

          '<div class=\"tab-pane\" id=\"meta-'+uid+'\">MetaSVM: '+data['MetaSVM_score']+' ('+data['MetaSVM_pred']+')<br>MetaLR: '+data['MetaLR_score']+' ('+data['MetaLR_pred']+')</div>' +

          '<div class=\"tab-pane\" id=\"cadd-'+uid+'\">Raw: '+data['CADD_raw']+'<br>Phred: '+data['CADD_phred']+'</div>' +

        '</div>' +
        '</div>' +
        
  '<div class=\"tab-pane\" id=\"trans-'+uid+'\">' +

  '<b>Transcript position</b><br><br>' +

  'cDNA position: '+data['cDNA_position']+'<br>' +
  'CDS position: '+data['CDS_position']+'<br>' +
  'Protein position: '+data['Protein_position']+'<br><br>' +

  '<b>Protein change</b><br><br>' +

  'Amino acids: '+data['Amino_acids']+'<br>' +
  'Codons: '+data['Codons']+'<br><br>' +

  '<b>Genomic location</b><br><br>' +

  'Exon: '+data['EXON']+'<br>' +
  'Intron: '+data['INTRON'] +

'</div>' +

'<div class=\"tab-pane\" id=\"annot-'+uid+'\">' +

  '<b>Transcript quality</b><br><br>' +

  'Canonical: '+data['CANONICAL']+'<br>' +
  'MANE Select: '+data['MANE_sel']+'<br>' +
  'MANE Plus Clinical: '+data['MANE_plus']+'<br>' +
  'TSL: '+data['TSL']+'<br>' +
  'APPRIS: '+data['APPRIS']+'<br><br>' +

  '<b>Transcript identifiers</b><br><br>' +

  'CCDS: '+data['CCDS']+'<br>' +
  'ENSP: '+data['ENSP']+'<br>' +
  'SwissProt: '+data['SWISSPROT']+'<br>' +
  'TrEMBL: '+data['TREMBL']+'<br>' +
  'UniParc: '+data['UNIPARC'] +

'</div>' +

'<div class=\"tab-pane\" id=\"hgvs-'+uid+'\">' +

  '<b>VEP HGVS</b><br><br>' +

  'HGVSc: '+data['HGVSc']+'<br>' +
  'HGVSp: '+data['HGVSp']+'<br>' +
  'HGVS offset: '+data['HGVS_OFFSET']+'<br><br>' +

  '<b>snpEff HGVS</b><br><br>' +

  'HGVSc snpEff: '+data['HGVSc_snpEff']+'<br>' +
  'HGVSp snpEff: '+data['HGVSp_snpEff'] +

'</div>' +

'<div class=\"tab-pane\" id=\"motif-'+uid+'\">' +

  '<b>Regulación y motifs</b><br><br>' +

  'Motif name: '+data['MOTIF_NAME']+'<br>' +
  'Motif position: '+data['MOTIF_POS']+'<br>' +
  'High information position: '+data['HIGH_INF_POS']+'<br>' +
  'Motif score change: '+data['MOTIF_SCORE_CHANGE']+'<br>' +
  'Transcription factors: '+data['TRANSCRIPTION_FACTORS'] +

'</div>' +


      '</div>' +

    '</div>' +

    '<hr>' +

    '<div class=\"explore-box\">' +
      '<div class=\"explore-title\">Explore</div>' +
      '<div class=\"explore-desc\">Navigate to gene viewer or RNA data table</div>' +
      '<div style=\"display:flex; gap:10px; margin-top:8px;\">' +
      '<button class=\"go-gene\" data-gene=\"'+gene+'\">View gene info</button><br>' +
      '<button class=\"go-rna\" data-gene=\"'+gene+'\">Go to RNA Data table</button>' +
      '</div>' +
    '</div>' +

  '</div>';
};

// expand
table.on('click','tr',function(){
  var tr=$(this);
  var row=table.row(tr);
  if(row.child.isShown()){
    row.child.hide();
    tr.removeClass('shown');
  } else {
    row.child(format(row.data())).show();
    tr.addClass('shown');
    var detail = tr.next('tr').find('.variant-detail');
    Shiny.setInputValue('", session$ns("variant_detail_request"), "', {
      row_key: row.data()[0],
      uid: detail.data('uid')
    }, {priority:'event'});
  }
});

// tabs
table.on('click','.freq-tab',function(e){
  e.stopPropagation();
  var container = $(this).closest('.variant-detail');
  container.find('.variant-detail-tabs .nav-link').removeClass('active');
  $(this).addClass('active');
  var target = $(this).data('target');
  container.children('.tab-content').children('.tab-pane').removeClass('active');
  container.find(target).addClass('active');
});

table.on('click','.subtab',function(e){
  e.stopPropagation();
  var container = $(this).closest('div');
  container.find('.subtab').removeClass('active');
  $(this).addClass('active');
  var target = $(this).data('target');
  container.find('.tab-pane').removeClass('active');
  container.find(target).addClass('active');
});

// navegación
table.on('click','.go-gene',function(e){
  e.stopPropagation();
  var gene=$(this).data('gene');
  Shiny.setInputValue('", session$ns("nav_click"), "',{gene:gene,tab:'gene'},{priority:'event'});
  Shiny.setInputValue('menu','gene',{priority:'event'});
});

table.on('click','.go-rna',function(e){
  e.stopPropagation();
  var gene=$(this).data('gene');
  Shiny.setInputValue('", session$ns("nav_click"), "',{gene:gene,tab:'rna'},{priority:'event'});
});

Shiny.addCustomMessageHandler('variant_detail_render', function(msg) {
  $('#variant-detail-body-' + msg.uid).html(msg.html);
  setTimeout(function(){
    table.columns.adjust();
  }, 100);
});

"
          ))
          
          
        )}, error = function(e){
          print(e)
          datatable(data.frame(Message="Error loading variants"))
        })
    })
    
    output$download_variants <- downloadHandler(
      filename = function(){
        ext <- input$variants_download_format %||% "tsv"
        paste0("variants_filtered_", Sys.Date(), ".", ext)
      },
      content = function(file){
        df <- filtered_variants_data()
        write_download_table(df, file)
      }
    )
    
    observeEvent(input$variant_detail_request, {
      req(input$variant_detail_request$row_key, input$variant_detail_request$uid)
      
      row <- variants_all() %>%
        filter(.variant_row_key == as.integer(input$variant_detail_request$row_key)) %>%
        slice(1)
      
      html <- if(nrow(row) == 0){
        "<i>No detail available for this variant</i>"
      } else {
        variant_detail_html(row)
      }
      
      session$sendCustomMessage(
        "variant_detail_render",
        list(
          uid = input$variant_detail_request$uid,
          html = html
        )
      )
    })
    
    # =====================
    # NAVEGACIÓN CONTROLADA DESDE BOTONES EN VARIANTES Y RNA
    # =====================
    observeEvent(input$nav_click, {
      
      navigating(TRUE)
      selected_gene(input$nav_click$gene)
      
      if(input$nav_click$tab == "rna"){
        updateTabsetPanel(session, "main_tabs", selected = "RNA Data")
        
      } else if(input$nav_click$tab == "drop"){
        updateTabsetPanel(session, "main_tabs", selected = "DROP")
        updateSelectInput(session, "drop_type", selected = input$nav_click$type)
        
      } else if(input$nav_click$tab == "variants"){
        updateTabsetPanel(session, "main_tabs", selected = "Variants (WES + WGS)")
      }
    })
    
    rna_files <- reactivePoll(
      
      intervalMillis = 1000,
      session = session,
      
      checkFunc = function(){
        
        files <- Sys.glob("../data/rnaseq/*.parquet")
        
        if(length(files)==0){
          return("")
        }
        
        info <- file.info(files)
        
        paste(
          files,
          info$mtime,
          info$size,
          collapse="|"
        )
      },
      
      valueFunc=function(){
        Sys.glob("../data/rnaseq/*.parquet")
      }
      
    )
    
    rna_filtered_data <- reactive({
      
      rna_files()
      
      df <- get_rna(pool)
      
      # ===== TPM RENOMBRE =====
      tpm_col <- grep("^tpm$", names(df), value=TRUE)
      
      if(length(tpm_col)==1){
        
        files_sample <- list.files(
          "../data/rnaseq",
          pattern="sample_",
          full.names=TRUE
        )
        
        sample_name <- tools::file_path_sans_ext(
          basename(files_sample[1])
        )
        
        sample_name <- sub(
          "_[0-9]{8}_[0-9]{6}$",
          "",
          sample_name
        )
        
        names(df)[names(df)=="tpm"] <- paste0(
          sample_name,
          "_tpm"
        )
      }
      
      # ===== DROP =====
      drop_flags <- get_drop_flags(pool)
      
      df <- df %>%
        left_join(
          drop_flags,
          by=c("gene_name"="gene")
        ) %>%
        mutate(
          drop_expr=replace_na(drop_expr,FALSE),
          drop_splicing=replace_na(drop_splicing,FALSE),
          drop_mae=replace_na(drop_mae,FALSE),
          
          DROP_status=paste0(
            "expr: ",drop_expr,"<br>",
            "splicing: ",drop_splicing,"<br>",
            "mae: ",drop_mae
          )
        )
      
      # ===== GENE FILTER =====
      if(!is.null(selected_gene()) && selected_gene()!=""){
        df <- df %>%
          filter(
            toupper(gene_name)==toupper(selected_gene())
          )
      }
      
      if(nzchar(input$gene)){
        df <- df %>%
          filter(
            toupper(gene_name)==toupper(input$gene)
          )
      }
      
      # ===== DROP FILTERS =====
      if("expr" %in% input$drop_filters){
        df <- df %>% filter(drop_expr)
      }
      
      if("splicing" %in% input$drop_filters){
        df <- df %>% filter(drop_splicing)
      }
      
      if("mae" %in% input$drop_filters){
        df <- df %>% filter(drop_mae)
      }
      
      # ===== TPM FILTERS =====
      
      sample_col <- grep(
        "_tpm$",
        names(df),
        value=TRUE
      )
      
      if(length(sample_col)==1){
        
        sample_col <- sample_col[1]
        
        df <- df %>%
          mutate(
            sample=as.numeric(.data[[sample_col]]),
            minv=as.numeric(.data[["min_TPM"]]),
            maxv=as.numeric(.data[["max_TPM"]])
          )
        
        if(identical(input$tpm_filters,"over")){
          df <- df %>%
            filter(sample>maxv)
        }
        
        if(identical(input$tpm_filters,"under")){
          df <- df %>%
            filter(sample<minv)
        }
        
        if(length(input$tpm_filters)==2){
          df <- df %>%
            filter(
              sample>maxv |
                sample<minv
            )
        }
        
        df <- df %>%
          select(
            -sample,
            -minv,
            -maxv
          )
      }
      
      df
    })
    
    # =====================
    # RNA DATA
    # =====================
    output$rna <- renderDT({
      
      rna_files()
      
      df <- rna_filtered_data()
      
      if(nrow(df)==0){
        return(datatable(data.frame(Message="No RNA data")))
      }
      
      # ===== RENOMBRAR TPM DINÁMICAMENTE =====
      
      tpm_col <- grep("^tpm$", colnames(df), value = TRUE)
      
      if(length(tpm_col) == 1){
        
        files_sample <- list.files("../data/rnaseq", pattern="sample_", full.names=TRUE)
        
        sample_name <- tools::file_path_sans_ext(basename(files_sample[1]))
        sample_name <- sub("_[0-9]{8}_[0-9]{6}$", "", sample_name)
        
        new_name <- paste0(sample_name, "_tpm")
        
        colnames(df)[colnames(df) == "tpm"] <- new_name
      }
      
 
      
      # ===== FILTROS =====
      if(!is.null(selected_gene()) && selected_gene()!=""){
        df <- df %>% filter(toupper(gene_name) == toupper(selected_gene()))
      }
      
      if(nzchar(input$gene)){
        df <- df %>% filter(toupper(gene_name) == toupper(input$gene))
      }
      
  
      
      # ===== HIDDEN COLUMN =====
      df$gene_hidden <- df$gene_name
      
      df <- df %>% select(-drop_expr, -drop_splicing, -drop_mae)
      
      datatable(
        df,
        rownames = FALSE,
        escape = FALSE,
        selection = "none",
        extensions = 'FixedHeader',
        options=list(
          scrollX=TRUE,
          scrollCollapse=TRUE,
          pageLength=10,
          fixedHeader=TRUE,
          autoWidth=TRUE,
          
          initComplete=JS("function(settings,json){this.api().columns.adjust();}"),
          
          drawCallback=JS("function(settings){this.api().columns.adjust();}"),
          
          columnDefs=list(
            list(
              targets=which(colnames(df)=="gene_hidden")-1,
              visible=FALSE
            )
          )
        ),
        
        callback = JS(sprintf("
      
      var format=function(rowData){

        var gene=rowData[rowData.length-1];
        var uid = Math.random().toString(36).substring(2,9);

        return '<div style=\"padding:10px\">' +

          '<b>DROP summary</b><br><br>' +

          '<ul class=\"nav nav-tabs\" role=\"tablist\">' +

            '<li class=\"nav-item\">' +
              '<a class=\"nav-link active drop-tab\" data-gene=\"'+gene+'\" data-type=\"expr\" data-uid=\"'+uid+'\" href=\"#\">Expression</a>' +
            '</li>' +

            '<li class=\"nav-item\">' +
              '<a class=\"nav-link drop-tab\" data-gene=\"'+gene+'\" data-type=\"splicing\" data-uid=\"'+uid+'\" href=\"#\">Splicing</a>' +
            '</li>' +

            '<li class=\"nav-item\">' +
              '<a class=\"nav-link drop-tab\" data-gene=\"'+gene+'\" data-type=\"mae\" data-uid=\"'+uid+'\" href=\"#\">MAE</a>' +
            '</li>' +

          '</ul>' +

          '<div class=\"tab-content\" style=\"margin-top:10px\">' +
            '<div class=\"tab-pane active\" id=\"expr-'+uid+'\">Click tab</div>' +
            '<div class=\"tab-pane\" id=\"splicing-'+uid+'\"></div>' +
            '<div class=\"tab-pane\" id=\"mae-'+uid+'\"></div>' +
          '</div>' +

          '<hr>' +

          '<div class=\"explore-box\">' +
            '<div class=\"explore-title\">Explore variants</div>' +
            '<div class=\"explore-desc\">Go to variants table filtered by this gene</div>' +
            '<button class=\"go-var\" data-gene=\"'+gene+'\">View variants</button>' +
          '</div>' +

        '</div>';
      };

      table.on('click','tr',function(){
        var tr=$(this); var row=table.row(tr);
        if(row.child.isShown()){
          row.child.hide();tr.removeClass('shown');
        } else {
          row.child(format(row.data())).show();tr.addClass('shown');
        }
      });

      table.on('click','.drop-tab',function(e){
        e.preventDefault();
        e.stopPropagation();
        
        var gene=$(this).data('gene');
        var type=$(this).data('type');
        var uid=$(this).data('uid');
        
        var target = '#'+type+'-'+uid;
        var container = $(target);
        
        $(this).closest('ul').find('.nav-link').removeClass('active');
        $(this).addClass('active');
        
        container.siblings().removeClass('active');
        container.addClass('active');
        
        container.html('Loading...');
        
        Shiny.setInputValue('%s',{
          gene:gene,
          type:type,
          uid:uid
        },{priority:'event'});
      });

      table.on('click','.go-var',function(e){
        e.stopPropagation();
        var gene=$(this).data('gene');
        Shiny.setInputValue('%s',{gene:gene,tab:'variants'},{priority:'event'});
      });

      Shiny.addCustomMessageHandler('drop_render', function(msg) {
        $(msg.target).html(msg.html);
      });

      setTimeout(function(){
        table.columns.adjust();
      }, 300);

    ",
                              session$ns("drop_request"),
                              session$ns("nav_click")
        ))
      )
    })
    
    output$download_rna <- downloadHandler(
      filename = function(){
        ext <- input$rna_download_format %||% "tsv"
        paste0("rna_table_filtered_", Sys.Date(), ".", ext)
      },
      content = function(file){
        df <- rna_filtered_data()
        if("gene_hidden" %in% colnames(df)){
          df <- df %>% select(-gene_hidden)
        }
        write_download_table(df, file)
      }
    )
    
    # =====================
    # DROP LOADER (FUERA)
    # =====================
    observeEvent(input$drop_request, {
      
      gene <- input$drop_request$gene
      type  <- input$drop_request$type
      uid   <- input$drop_request$uid
      
      target <- paste0("#", type, "-", uid)
      
      df <- switch(
        type,
        expr = get_drop_expr(pool),
        splicing = get_drop_splicing(pool),
        mae = get_drop_mae(pool)
      )
      
      # =====================
      # NO HAY DATOS DROP
      # =====================
      
      if(nrow(df) == 0){
        
        html <- "<i>No DROP data loaded</i>"
        
        session$sendCustomMessage(
          "drop_render",
          list(
            html = html,
            target = target
          )
        )
        
        return()
      }
      
      # =====================
      # DETECTAR COLUMNA GEN
      # =====================
      
      gene_col <- NULL
      
      if("gene_name" %in% colnames(df)){
        gene_col <- "gene_name"
      }
      
      if("hgncSymbol" %in% colnames(df)){
        gene_col <- "hgncSymbol"
      }
      
      if(is.null(gene_col)){
        
        html <- "<i>Gene column not found</i>"
        
        session$sendCustomMessage(
          "drop_render",
          list(
            html = html,
            target = target
          )
        )
        
        return()
      }
      
      # =====================
      # FILTRAR GEN
      # =====================
      
      df <- df %>%
        filter(.data[[gene_col]] == gene)
      
      # =====================
      # NO HAY RESULTADOS
      # =====================
      
      if(nrow(df) == 0){
        
        html <- "<i>No data for this gene</i>"
        
      } else {
        
        html <- paste0(
          
          "<div style='width:100%;overflow-x:auto;'>",
          
          "<div style='width:2500px;'>",
          
          "<table style='font-size:12px;
             border-collapse:collapse;
             white-space:nowrap;'>",
          
          "<thead><tr>",
          
          paste0(
            "<th style='padding:6px 10px;border:1px solid #ddd'>",
            colnames(df),
            "</th>",
            collapse=""
          ),
          
          "</tr></thead>",
          
          "<tbody>",
          
          paste0(
            apply(df,1,function(row){
              
              paste0(
                "<tr>",
                
                paste0(
                  "<td style='padding:6px 10px;border:1px solid #ddd'>",
                  row,
                  "</td>",
                  collapse=""
                ),
                
                "</tr>"
              )
              
            }),
            collapse=""
          ),
          
          "</tbody>",
          "</table>",
          
          "</div>",
          "</div>"
        )
      }
      
      session$sendCustomMessage(
        "drop_render",
        list(
          html = html,
          target = target
        )
      )
    })
    
  })
}
