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
          CONCAT(
            COALESCE(CAST(CHROM AS VARCHAR), ''), ':',
            COALESCE(CAST(POS AS VARCHAR), ''), ':',
            COALESCE(CAST(REF AS VARCHAR), ''), ':',
            COALESCE(CAST(ALT AS VARCHAR), ''), ':',
            COALESCE(CAST(Gene AS VARCHAR), ''), ':',
            COALESCE(CAST(ENSP AS VARCHAR), '')
          ) AS variant_source_key,
        
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
      ),

      ranked AS (
        SELECT *,
          ROW_NUMBER() OVER (
            PARTITION BY variant_source_key
            ORDER BY CASE WHEN src = 'WGS' THEN 1 ELSE 2 END
          ) AS source_rank,
          COUNT(DISTINCT src) OVER (PARTITION BY variant_source_key) AS source_count,
          MAX(CASE WHEN src = 'WES' THEN CHILD_GT END) OVER (PARTITION BY variant_source_key) AS WES_CHILD_GT,
          MAX(CASE WHEN src = 'WES' THEN CHILD_DP END) OVER (PARTITION BY variant_source_key) AS WES_CHILD_DP,
          MAX(CASE WHEN src = 'WES' THEN CHILD_AD END) OVER (PARTITION BY variant_source_key) AS WES_CHILD_AD,
          MAX(CASE WHEN src = 'WES' THEN CHILD_GQ END) OVER (PARTITION BY variant_source_key) AS WES_CHILD_GQ,
          MAX(CASE WHEN src = 'WES' THEN PARENT1_GT END) OVER (PARTITION BY variant_source_key) AS WES_PARENT1_GT,
          MAX(CASE WHEN src = 'WES' THEN PARENT1_DP END) OVER (PARTITION BY variant_source_key) AS WES_PARENT1_DP,
          MAX(CASE WHEN src = 'WES' THEN PARENT1_AD END) OVER (PARTITION BY variant_source_key) AS WES_PARENT1_AD,
          MAX(CASE WHEN src = 'WES' THEN PARENT1_GQ END) OVER (PARTITION BY variant_source_key) AS WES_PARENT1_GQ,
          MAX(CASE WHEN src = 'WES' THEN PARENT2_GT END) OVER (PARTITION BY variant_source_key) AS WES_PARENT2_GT,
          MAX(CASE WHEN src = 'WES' THEN PARENT2_DP END) OVER (PARTITION BY variant_source_key) AS WES_PARENT2_DP,
          MAX(CASE WHEN src = 'WES' THEN PARENT2_AD END) OVER (PARTITION BY variant_source_key) AS WES_PARENT2_AD,
          MAX(CASE WHEN src = 'WES' THEN PARENT2_GQ END) OVER (PARTITION BY variant_source_key) AS WES_PARENT2_GQ,
          MAX(CASE WHEN src = 'WGS' THEN CHILD_GT END) OVER (PARTITION BY variant_source_key) AS WGS_CHILD_GT,
          MAX(CASE WHEN src = 'WGS' THEN CHILD_DP END) OVER (PARTITION BY variant_source_key) AS WGS_CHILD_DP,
          MAX(CASE WHEN src = 'WGS' THEN CHILD_AD END) OVER (PARTITION BY variant_source_key) AS WGS_CHILD_AD,
          MAX(CASE WHEN src = 'WGS' THEN CHILD_GQ END) OVER (PARTITION BY variant_source_key) AS WGS_CHILD_GQ,
          MAX(CASE WHEN src = 'WGS' THEN PARENT1_GT END) OVER (PARTITION BY variant_source_key) AS WGS_PARENT1_GT,
          MAX(CASE WHEN src = 'WGS' THEN PARENT1_DP END) OVER (PARTITION BY variant_source_key) AS WGS_PARENT1_DP,
          MAX(CASE WHEN src = 'WGS' THEN PARENT1_AD END) OVER (PARTITION BY variant_source_key) AS WGS_PARENT1_AD,
          MAX(CASE WHEN src = 'WGS' THEN PARENT1_GQ END) OVER (PARTITION BY variant_source_key) AS WGS_PARENT1_GQ,
          MAX(CASE WHEN src = 'WGS' THEN PARENT2_GT END) OVER (PARTITION BY variant_source_key) AS WGS_PARENT2_GT,
          MAX(CASE WHEN src = 'WGS' THEN PARENT2_DP END) OVER (PARTITION BY variant_source_key) AS WGS_PARENT2_DP,
          MAX(CASE WHEN src = 'WGS' THEN PARENT2_AD END) OVER (PARTITION BY variant_source_key) AS WGS_PARENT2_AD,
          MAX(CASE WHEN src = 'WGS' THEN PARENT2_GQ END) OVER (PARTITION BY variant_source_key) AS WGS_PARENT2_GQ
        FROM (
          SELECT DISTINCT *
          FROM normalized
        )
      )
      
      SELECT 
        * EXCLUDE (variant_source_key, source_rank, source_count, CHILD_DP, CHILD_AD, CHILD_GQ),

        -- MAIN TABLE METRICS
        CASE
          WHEN source_count = 2 THEN GREATEST(
            TRY_CAST(WES_CHILD_DP AS DOUBLE),
            TRY_CAST(WGS_CHILD_DP AS DOUBLE)
          )
          ELSE TRY_CAST(CHILD_DP AS DOUBLE)
        END AS CHILD_DP,

        CASE
          WHEN source_count = 2
               AND TRY_CAST(WES_CHILD_AD AS DOUBLE) IS NOT NULL
               AND TRY_CAST(WGS_CHILD_AD AS DOUBLE) IS NOT NULL
          THEN CAST(GREATEST(
            TRY_CAST(WES_CHILD_AD AS DOUBLE),
            TRY_CAST(WGS_CHILD_AD AS DOUBLE)
          ) AS VARCHAR)
          WHEN source_count = 2
               AND COALESCE(TRY_CAST(WES_CHILD_DP AS DOUBLE), -1) >=
                   COALESCE(TRY_CAST(WGS_CHILD_DP AS DOUBLE), -1)
          THEN WES_CHILD_AD
          WHEN source_count = 2 THEN WGS_CHILD_AD
          ELSE CHILD_AD
        END AS CHILD_AD,

        CASE
          WHEN source_count = 2 THEN GREATEST(
            TRY_CAST(WES_CHILD_GQ AS DOUBLE),
            TRY_CAST(WGS_CHILD_GQ AS DOUBLE)
          )
          ELSE TRY_CAST(CHILD_GQ AS DOUBLE)
        END AS CHILD_GQ,
      
        -- SOURCE
        CASE
          WHEN source_count = 2 THEN 'BOTH'
          ELSE src
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
        
      FROM ranked
      WHERE source_rank = 1
    "
    
  } else if(has_wes){
    query <- single_source_query("../data/variants/wes*.parquet", "WES")
    
  } else if(has_wgs){
    query <- single_source_query("../data/variants/wgs*.parquet", "WGS")
  }
  
  df <- dbGetQuery(pool, query)
  
  return(df)
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
