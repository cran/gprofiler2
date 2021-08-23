## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
#  install.packages("gprofiler2")

## ----setup--------------------------------------------------------------------
library(gprofiler2)

## -----------------------------------------------------------------------------
gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

## -----------------------------------------------------------------------------
names(gostres)

## -----------------------------------------------------------------------------
head(gostres$result, 3)

## -----------------------------------------------------------------------------
names(gostres$meta)

## -----------------------------------------------------------------------------
gostres2 <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL)

## -----------------------------------------------------------------------------
head(gostres2$result, 3)

## ----eval = FALSE-------------------------------------------------------------
#  gostres_link <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"),
#                  as_short_link = TRUE)

## -----------------------------------------------------------------------------
multi_gostres1 <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340", 
                                                 "GO:0005005", "ENSG00000156103", "NLRP1"),
                             "chromY" = c("Y:1:10000000", "rs17396340", 
                                          "GO:0005005", "ENSG00000156103", "NLRP1")), 
                       multi_query = FALSE)

## -----------------------------------------------------------------------------
head(multi_gostres1$result, 3)

## -----------------------------------------------------------------------------
multi_gostres2 <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340",
                                                 "GO:0005005", "ENSG00000156103", "NLRP1"),
                             "chromY" = c("Y:1:10000000", "rs17396340", 
                                          "GO:0005005", "ENSG00000156103", "NLRP1")), 
                       multi_query = TRUE)

## -----------------------------------------------------------------------------
head(multi_gostres2$result, 3)

## ----fig.width = 9.5----------------------------------------------------------
gostplot(gostres, capped = TRUE, interactive = TRUE)

## ----fig.width = 9.5, fig.height = 4------------------------------------------
p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p

## ----fig.width = 9.5----------------------------------------------------------
pp <- publish_gostplot(p, highlight_terms = c("GO:0048013", "REAC:R-HSA-3928663"), 
                       width = NA, height = NA, filename = NULL )

## ----fig.width = 9.5, fig.height = 3------------------------------------------
publish_gosttable(gostres, highlight_terms = gostres$result[c(1:2,10,120),],
                        use_colors = TRUE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_size"),
                        filename = NULL)

## ----fig.width = 9.5, fig.height = 7, warning = F-----------------------------
gostplot(multi_gostres2, capped = TRUE, interactive = TRUE)

## ----fig.width = 10, fig.height = 2, warning = F------------------------------
publish_gosttable(multi_gostres1, 
                         highlight_terms = multi_gostres1$result[c(1, 82, 176),],
                        use_colors = TRUE, 
                        show_columns = c("source", "term_name", "term_size"),
                        filename = NULL)

## -----------------------------------------------------------------------------
get_version_info(organism = "hsapiens")

## ----eval = F-----------------------------------------------------------------
#  download.file(url = "http://software.broadinstitute.org/gsea/resources/msigdb/7.0/c2.cp.biocarta.v7.0.symbols.gmt", destfile = "extdata/biocarta.gmt")

## ----eval = F-----------------------------------------------------------------
#  upload_GMT_file(gmtfile = "extdata/biocarta.gmt")

## -----------------------------------------------------------------------------
custom_gostres <- gost(query = c("MAPK3",	"PIK3C2G", "HRAS", "PIK3R1", "MAP2K1", 
                                 "RAF1", "PLCG1",	"GNAQ",	"MAPK1", "PRKCB",	"CRK", "BCAR1", "NFKB1"),
                       organism = "gp__TEXF_hZLM_d18")
head(custom_gostres$result, 3)

## -----------------------------------------------------------------------------
gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                evcodes = TRUE, multi_query = FALSE, 
                sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"))

gem <- gostres$result[,c("term_id", "term_name", "p_value", "intersection")]
colnames(gem) <- c("GO.ID", "Description", "p.Val", "Genes")
gem$FDR <- gem$p.Val
gem$Phenotype = "+1"
gem <- gem[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]
head(gem, 3)

## ---- eval=F------------------------------------------------------------------
#  write.table(gem, file = "extdata/gProfiler_gem.txt", sep = "\t", quote = F, row.names = F)

## ---- eval=F------------------------------------------------------------------
#  # enrichment for two input gene lists
#  multi_gostres <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340",
#                                                   "GO:0005005", "ENSG00000156103", "NLRP1"),
#                               "chromY" = c("Y:1:10000000", "rs17396340",
#                                            "GO:0005005", "ENSG00000156103", "NLRP1")),
#                        evcodes = TRUE, multi_query = FALSE,
#                        sources = c("GO", "REAC", "MIRNA", "CORUM", "HP", "HPA", "WP"))
#  
#  # format to GEM
#  gem <- multi_gostres$result[,c("query", "term_id", "term_name", "p_value", "intersection")]
#  colnames(gem) <- c("query", "GO.ID", "Description", "p.Val", "Genes")
#  gem$FDR <- gem$p.Val
#  gem$Phenotype = "+1"
#  
#  # write separate files for queries
#  
#  # install.packages("dplyr")
#  library(dplyr)
#  
#  gem %>% group_by(query) %>%
#    group_walk(~
#      write.table(data.frame(.x[,c("GO.ID", "Description", "p.Val", "FDR", "Phenotype", "Genes")]),
#                  file = paste0("gProfiler_", unique(.y$query), "_gem.txt"),
#                  sep = "\t", quote = F, row.names = F))

## -----------------------------------------------------------------------------
gconvert(query = c("GO:0005030", "rs17396340", "NLRP1"), organism = "hsapiens", 
         target="ENSG", mthreshold = Inf, filter_na = TRUE)

## -----------------------------------------------------------------------------
gorth(query = c("Klf4", "Sox2", "71950"), source_organism = "mmusculus", 
      target_organism = "hsapiens", mthreshold = Inf, filter_na = TRUE,
      numeric_ns = "ENTREZGENE_ACC")

## -----------------------------------------------------------------------------
gsnpense(query = c("rs11734132", "rs4305276", "rs17396340", "rs3184504"), 
         filter_na = TRUE)

## -----------------------------------------------------------------------------
gsnpense(query = c("rs11734132", "rs4305276", "rs17396340", "rs3184504"), 
         filter_na = FALSE)

## -----------------------------------------------------------------------------
set_base_url("http://biit.cs.ut.ee/gprofiler_beta")

## -----------------------------------------------------------------------------
get_base_url()

## -----------------------------------------------------------------------------
set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e95_eg42_p13")

