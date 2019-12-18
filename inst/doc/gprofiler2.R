## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE--------------------------------------------------------
#  install.packages("gprofiler2")

## ----setup---------------------------------------------------------------
library(gprofiler2)

## ------------------------------------------------------------------------
gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

## ------------------------------------------------------------------------
names(gostres)

## ------------------------------------------------------------------------
head(gostres$result)

## ------------------------------------------------------------------------
names(gostres$meta)

## ------------------------------------------------------------------------
gostres2 <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL)

## ------------------------------------------------------------------------
head(gostres2$result)

## ----eval = FALSE--------------------------------------------------------
#  gostres_link <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"),
#                  as_short_link = TRUE)

## ------------------------------------------------------------------------
multi_gostres1 <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340", 
                                                 "GO:0005005", "ENSG00000156103", "NLRP1"),
                             "chromY" = c("Y:1:10000000", "rs17396340", 
                                          "GO:0005005", "ENSG00000156103", "NLRP1")), 
                       multi_query = FALSE)

## ------------------------------------------------------------------------
head(multi_gostres1$result)

## ------------------------------------------------------------------------
multi_gostres2 <- gost(query = list("chromX" = c("X:1000:1000000", "rs17396340",
                                                 "GO:0005005", "ENSG00000156103", "NLRP1"),
                             "chromY" = c("Y:1:10000000", "rs17396340", 
                                          "GO:0005005", "ENSG00000156103", "NLRP1")), 
                       multi_query = TRUE)

## ------------------------------------------------------------------------
head(multi_gostres2$result)

## ----fig.width = 9.5-----------------------------------------------------
gostplot(gostres, capped = TRUE, interactive = TRUE)

## ----fig.width = 9.5, fig.height = 4-------------------------------------
p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p

## ----fig.width = 9.5-----------------------------------------------------
pp <- publish_gostplot(p, highlight_terms = c("GO:0048013", "REAC:R-HSA-3928663"), 
                       width = NA, height = NA, filename = NULL )

## ----fig.width = 9.5, fig.height = 3-------------------------------------
pt <- publish_gosttable(gostres, highlight_terms = gostres$result[c(1:2,10,100:102,120,124,125),],
                        use_colors = TRUE, 
                        show_columns = c("source", "term_name", "term_size", "intersection_size"),
                        filename = NULL)

## ----fig.width = 9.5, fig.height = 7, warning = F------------------------
gostplot(multi_gostres2, capped = TRUE, interactive = TRUE)

## ----fig.width = 10, fig.height = 2, warning = F-------------------------
pt2 <- publish_gosttable(multi_gostres1, 
                         highlight_terms = multi_gostres1$result[c(1, 24, 82, 176, 204, 234),],
                        use_colors = TRUE, 
                        show_columns = c("source", "term_name", "term_size"),
                        filename = NULL)

## ----eval = F------------------------------------------------------------
#  download.file(url = "http://software.broadinstitute.org/gsea/resources/msigdb/7.0/c2.cp.biocarta.v7.0.symbols.gmt", destfile = "extdata/biocarta.gmt")

## ----eval = F------------------------------------------------------------
#  upload_GMT_file(gmtfile = "extdata/biocarta.gmt")

## ------------------------------------------------------------------------
custom_gostres <- gost(query = c("MAPK3",	"PIK3C2G", "HRAS", "PIK3R1", "MAP2K1", 
                                 "RAF1", "PLCG1",	"GNAQ",	"MAPK1", "PRKCB",	"CRK", "BCAR1", "NFKB1"),
                       organism = "gp__TEXF_hZLM_d18")
head(custom_gostres$result)

## ------------------------------------------------------------------------
gconvert(query = c("REAC:R-HSA-3928664", "rs17396340", "NLRP1"), organism = "hsapiens", 
         target="ENSG", mthreshold = Inf, filter_na = TRUE)

## ------------------------------------------------------------------------
gorth(query = c("Klf4", "Sox2", "71950"), source_organism = "mmusculus", 
      target_organism = "hsapiens", mthreshold = Inf, filter_na = TRUE,
      numeric_ns = "ENTREZGENE_ACC")

## ------------------------------------------------------------------------
gsnpense(query = c("rs11734132", "rs4305276", "rs17396340", "rs3184504"), 
         filter_na = TRUE)

## ------------------------------------------------------------------------
gsnpense(query = c("rs11734132", "rs4305276", "rs17396340", "rs3184504"), 
         filter_na = FALSE)

## ------------------------------------------------------------------------
set_base_url("http://biit.cs.ut.ee/gprofiler_beta")

## ------------------------------------------------------------------------
get_base_url()

## ------------------------------------------------------------------------
set_base_url("http://biit.cs.ut.ee/gprofiler_archive3/e95_eg42_p13")

