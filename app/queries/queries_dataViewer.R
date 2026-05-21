# queries_dataViewer.R

# =====================
# VARIANTS BASE
# =====================
get_variants_with_inheritance <- function(pool){
  
  files <- list.files("../data/variants", full.names = TRUE)
  
  if(length(files) == 0){
    stop("No variant files found")
  }
  
  has_wes <- any(grepl("wes", files, ignore.case = TRUE))
  has_wgs <- any(grepl("wgs", files, ignore.case = TRUE))

  single_source_query <- function(pattern, source){
    sprintf("
      WITH normalized AS (
        SELECT *,

          CASE
            WHEN REPLACE(CHILD_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(CHILD_GT, '|', '/')
          END AS CHILD_GT_N,

          CASE
            WHEN REPLACE(PARENT1_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(PARENT1_GT, '|', '/')
          END AS P1_GT_N,

          CASE
            WHEN REPLACE(PARENT2_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(PARENT2_GT, '|', '/')
          END AS P2_GT_N

        FROM read_parquet('%s')
      )

      SELECT
        *,
        '%s' AS source,

       CASE

        -- AR heredada
        WHEN inheritance_source='AR'
             AND CHILD_GT_N='1/1'
             AND P1_GT_N='0/1'
             AND P2_GT_N='0/1'
        THEN 'recessive'
      
        -- AD/XD de novo
        WHEN inheritance_source IN ('AD','XD')
             AND CHILD_GT_N='0/1'
             AND P1_GT_N='0/0'
             AND P2_GT_N='0/0'
        THEN 'de_novo'
      
        -- AD/XD heredada
        WHEN inheritance_source IN ('AD','XD')
             AND CHILD_GT_N='0/1'
             AND (
                  P1_GT_N='0/1'
                  OR P2_GT_N='0/1'
             )
        THEN 'dominant'
      
        -- XR
        WHEN inheritance_source='XR'
        THEN 'x_recessive'
      
        -- MT
        WHEN inheritance_source='MT'
        THEN 'mitochondrial'
      
        ELSE ''
      
      END AS inheritance_type,

        CONCAT(
          CASE
            WHEN inheritance_source IN ('AD','XD') AND CHILD_GT_N = '0/1' AND P1_GT_N = '0/0' AND P2_GT_N = '0/0' THEN 'de_novo'
            WHEN inheritance_source = 'AR' AND CHILD_GT_N = '1/1' AND P1_GT_N = '0/1' AND P2_GT_N = '0/1' THEN 'recessive'
            WHEN inheritance_source != 'AR' AND CHILD_GT_N = '0/1' AND (P1_GT_N = '0/1' OR P2_GT_N = '0/1') THEN 'dominant'
            ELSE ''
          END,
          ' (P1:', P1_GT_N, ' AD:', PARENT1_AD,
          ' | P2:', P2_GT_N, ' AD:', PARENT2_AD, ')'
        ) AS inheritance

      FROM (
        SELECT DISTINCT *
        FROM normalized
      )
    ", pattern, source)
  }
  
  # =====================
  # WES + WGS
  # =====================
  if(has_wes & has_wgs){
    
    query <- "
      WITH wes AS (
        SELECT *, 'WES' AS src FROM read_parquet('../data/variants/wes*.parquet')
      ),
      wgs AS (
        SELECT *, 'WGS' AS src FROM read_parquet('../data/variants/wgs*.parquet')
      ),
      combined AS (
        SELECT * FROM wes
        UNION ALL
        SELECT * FROM wgs
      ),
      
      normalized AS (
        SELECT *,
        
          CASE 
            WHEN REPLACE(CHILD_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(CHILD_GT, '|', '/')
          END AS CHILD_GT_N,
          
          CASE 
            WHEN REPLACE(PARENT1_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(PARENT1_GT, '|', '/')
          END AS P1_GT_N,
          
          CASE 
            WHEN REPLACE(PARENT2_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(PARENT2_GT, '|', '/')
          END AS P2_GT_N
          
        FROM combined
      )
      
      SELECT 
        *,
      
        -- SOURCE
        CASE
          WHEN COUNT(DISTINCT src) OVER (PARTITION BY ID) = 2 THEN 'BOTH'
          ELSE MAX(src) OVER (PARTITION BY ID)
        END AS source,
        
        -- INHERITANCE TYPE
        CASE

        -- AR heredada
        WHEN inheritance_source='AR'
             AND CHILD_GT_N='1/1'
             AND P1_GT_N='0/1'
             AND P2_GT_N='0/1'
        THEN 'recessive'
      
        -- AD/XD de novo
        WHEN inheritance_source IN ('AD','XD')
             AND CHILD_GT_N='0/1'
             AND P1_GT_N='0/0'
             AND P2_GT_N='0/0'
        THEN 'de_novo'
      
        -- AD/XD heredada
        WHEN inheritance_source IN ('AD','XD')
             AND CHILD_GT_N='0/1'
             AND (
                  P1_GT_N='0/1'
                  OR P2_GT_N='0/1'
             )
        THEN 'dominant'
      
        -- XR
        WHEN inheritance_source='XR'
        THEN 'x_recessive'
      
        -- MT
        WHEN inheritance_source='MT'
        THEN 'mitochondrial'
      
        ELSE ''
      
      END AS inheritance_type,
        
        -- INHERITANCE TEXT
        CONCAT(
          CASE
            WHEN inheritance_source IN ('AD','XD') AND CHILD_GT_N = '0/1' AND P1_GT_N = '0/0' AND P2_GT_N = '0/0' THEN 'de_novo'
            WHEN inheritance_source = 'AR' AND CHILD_GT_N = '1/1' AND P1_GT_N = '0/1' AND P2_GT_N = '0/1' THEN 'recessive'
            WHEN inheritance_source != 'AR' AND CHILD_GT_N = '0/1' AND (P1_GT_N = '0/1' OR P2_GT_N = '0/1') THEN 'dominant'
            ELSE ''
          END,
          ' (P1:', P1_GT_N, ' AD:', PARENT1_AD,
          ' | P2:', P2_GT_N, ' AD:', PARENT2_AD, ')'
        ) AS inheritance
        
      FROM (
        SELECT DISTINCT *
        FROM normalized
      )
    "
    
  } else if(has_wes){
    query <- single_source_query("../data/variants/wes*.parquet", "WES")
    
  } else if(has_wgs){
    query <- single_source_query("../data/variants/wgs*.parquet", "WGS")
  }
  
  df <- dbGetQuery(pool, query)
  
  return(df)
}


get_variants_sql_base <- function(){
  
  files <- list.files("../data/variants", full.names = TRUE)
  
  has_wes <- any(grepl("wes", files, ignore.case = TRUE))
  has_wgs <- any(grepl("wgs", files, ignore.case = TRUE))
  
  if(has_wes & has_wgs){
    
    return("
      WITH wes AS (
        SELECT *, 'WES' AS src FROM read_parquet('../data/variants/wes*.parquet')
      ),
      wgs AS (
        SELECT *, 'WGS' AS src FROM read_parquet('../data/variants/wgs*.parquet')
      ),
      combined AS (
        SELECT * FROM wes
        UNION ALL
        SELECT * FROM wgs
      ),
      
      normalized AS (
        SELECT *,
        
          CASE 
            WHEN REPLACE(CHILD_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(CHILD_GT, '|', '/')
          END AS CHILD_GT_N,
          
          CASE 
            WHEN REPLACE(PARENT1_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(PARENT1_GT, '|', '/')
          END AS P1_GT_N,
          
          CASE 
            WHEN REPLACE(PARENT2_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
            ELSE REPLACE(PARENT2_GT, '|', '/')
          END AS P2_GT_N
          
        FROM combined
      )
      
      SELECT *
      FROM normalized
    ")
    
  } else if(has_wes){
    
    return("
      SELECT *
      FROM read_parquet('../data/variants/wes*.parquet')
    ")
    
  } else {
    
    return("
      SELECT *
      FROM read_parquet('../data/variants/wgs*.parquet')
    ")
    
  }
}

get_variants_sql <- function(){
  
  return("
    WITH wes AS (
      SELECT *, 'WES' AS src FROM read_parquet('../data/variants/wes*.parquet')
    ),
    wgs AS (
      SELECT *, 'WGS' AS src FROM read_parquet('../data/variants/wgs*.parquet')
    ),
    combined AS (
      SELECT * FROM wes
      UNION ALL
      SELECT * FROM wgs
    ),
    
    normalized AS (
      SELECT *,
      
        CASE 
          WHEN REPLACE(CHILD_GT, '|', '/') IN ('1/0','0/1') THEN '0/1'
          ELSE REPLACE(CHILD_GT, '|', '/')
        END AS CHILD_GT_N
        
      FROM combined
    )
    
    SELECT *
    FROM normalized
  ")
}

# =====================
# RNA
# =====================
get_rna <- function(pool){
  
  files_sample <- list.files("../data/rnaseq", pattern="sample_", full.names=TRUE)
  files_ctrl   <- list.files("../data/rnaseq", pattern="controls", full.names=TRUE)
  
  if(length(files_sample) == 0){
    stop("No RNA sample data")
  }
  
  if(length(files_ctrl) == 0){
    
    df <- dbGetQuery(pool, "
      SELECT 
        gene_id,
        gene_name,
        tpm AS \"gene tpm\"
      FROM read_parquet('../data/rnaseq/sample_*.parquet')
    ")
    
  } else {
    
    df <- dbGetQuery(pool, "
      SELECT 
        r.gene_id,
        r.gene_name,
        r.tpm AS \"tpm\",
        c.max_TPM,
        c.min_TPM,
        c.mean_TPM,
        c.median_TPM
      FROM read_parquet('../data/rnaseq/sample_*.parquet') r
      LEFT JOIN read_parquet('../data/rnaseq/controls*.parquet') c
      USING(gene_id)
    ")
  }
  
  return(df)
}


# =====================
# DROP
# =====================

get_drop_expr <- function(pool){
  
  files <- Sys.glob("../data/drop/expression*.parquet")
  
  if(length(files) == 0){
    return(data.frame())
  }
  
  dbGetQuery(
    pool,
    "SELECT * FROM read_parquet('../data/drop/expression*.parquet')"
  )
}

get_drop_splicing <- function(pool){
  
  files <- Sys.glob("../data/drop/splicing*.parquet")
  
  if(length(files) == 0){
    return(data.frame())
  }
  
  dbGetQuery(
    pool,
    "SELECT * FROM read_parquet('../data/drop/splicing*.parquet')"
  )
}

get_drop_mae <- function(pool){
  
  files <- Sys.glob("../data/drop/mae*.parquet")
  
  if(length(files) == 0){
    return(data.frame())
  }
  
  dbGetQuery(
    pool,
    "SELECT * FROM read_parquet('../data/drop/mae*.parquet')"
  )
}
