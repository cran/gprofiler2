gp_globals = new.env()

gp_globals$version =
  tryCatch(
    utils::packageVersion("gprofiler2"),
    error = function(e) { return("unknown_version") }
  );

# Set SSL version to TLSv1_2 with fallback to TLSv1_1
# CURL_SSLVERSION_SSLv3 is not used due to the SSLv3 vulnerability <https://access.redhat.com/articles/1232123>
# CURL_SSLVERSION_TLSv1_3 is not widespread enough to have a built-in LibreSSL support yet.
# (curl's authors may decide to change it at some point, so links to the source are provided.)
gp_globals$CURL_SSLVERSION_TLSv1_1 <- 5L # <https://github.com/curl/curl/blob/master/include/curl/curl.h#L1925>
gp_globals$CURL_SSLVERSION_TLSv1_2 <- 6L # <https://github.com/curl/curl/blob/master/include/curl/curl.h#L1926>

gp_globals$rcurl_opts =
  RCurl::curlOptions(useragent = paste("gprofiler2/", gp_globals$version, sep=""), sslversion = gp_globals$CURL_SSLVERSION_TLSv1_2)
gp_globals$base_url = "http://biit.cs.ut.ee/gprofiler"


#' Gene list functional enrichment.
#'
#' Interface to the g:Profiler tool g:GOSt (\url{https://biit.cs.ut.ee/gprofiler/gost}) for functional enrichments analysis of gene lists.
#' In case the input 'query' is a list of gene vectors, results for multiple queries will be returned in the same data frame with column 'query' indicating the corresponding query name.
#' If 'multi_query' is selected, the result is a data frame for comparing multiple input lists,
#' just as in the web tool.
#'
#' @param query vector, or a (named) list of vectors for multiple queries, that can consist of mixed types of gene IDs (proteins, transcripts, microarray IDs, etc), SNP IDs, chromosomal intervals or term IDs.
#' @param organism organism name. Organism names are constructed by concatenating the first letter of the name and the
#' family name. Example: human - 'hsapiens', mouse - 'mmusculus'.
#' @param ordered_query in case input gene lists are ranked this option may be
#'  used to get GSEA style p-values.
#' @param multi_query in case of multiple gene lists, returns comparison table of these lists.
#' If enabled, the result data frame has columns named 'p_values', 'query_sizes', 'intersection_sizes' with vectors showing values in the order of input queries.
#' To get the results in a long format set 'multi_query' to FALSE and just input query list of multiple gene vectors.
#' @param significant whether all or only statistically significant results should
#'  be returned.
#' @param exclude_iea exclude GO electronic annotations (IEA).
#' @param measure_underrepresentation measure underrepresentation.
#' @param evcodes include evidence codes to the results. Note
#'  that this can decrease performance and make the query slower.
#'  In addition, a column 'intersection' is created that contains the gene id-s that intersect between the query and term.
#'  This parameter does not work if 'multi_query' is set to TRUE.
#' @param user_threshold custom p-value threshold, results with a larger p-value are
#'  excluded.
#' @param correction_method the algorithm used for multiple testing correction, one of "gSCS" (synonyms: "analytical", "g_SCS"), "fdr" (synonyms: "false_discovery_rate"), "bonferroni".
#' @param domain_scope how to define statistical domain, one of "annotated", "known" or "custom".
#' @param custom_bg vector of gene names to use as a statistical background. If given, the domain_scope is set to 'custom'.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param sources a vector of data sources to use. Currently, these include
#'  GO (GO:BP, GO:MF, GO:CC to select a particular GO branch), KEGG, REAC, TF,
#'  MIRNA, CORUM, HP, HPA, WP. Please see the g:GOSt web tool for the comprehensive
#'  list and details on incorporated data sources.
#' @return A named list where 'result' contains data.frame with the enrichment analysis results and 'meta' contains metadata needed for Manhattan plot. If the input
#'  consisted of several lists the corresponding list is indicated with a variable
#'  'query'.
#'  When requesting a 'multi_query', either TRUE or FALSE, the columns of the resulting data frame differ.
#'  If 'evcodes' is set, the return value includes columns 'evidence_codes' and 'intersection'.
#'  The latter conveys info about the intersecting genes between the corresponding query and term.
#'
#'  The result fields are further described in \url{https://biit.cs.ut.ee/gprofiler_beta/page/apis#gost_query_results}
#' @author  Liis Kolberg <liis.kolberg@@ut.ee>, Uku Raudvere <uku.raudvere@@ut.ee>
#' @examples
#' gostres <- gost(c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"))
#'
#' @export
gost <- function(query,
                      organism = "hsapiens",
                      ordered_query = FALSE,
                      multi_query = FALSE,
                      significant = TRUE,
                      exclude_iea = FALSE,
                      measure_underrepresentation = FALSE,
                      evcodes = FALSE,
                      user_threshold = 0.05,
                      correction_method = c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS", "analytical"),
                      domain_scope = c("annotated", "known", "custom"),
                      custom_bg = NULL,
                      numeric_ns  = "",
                      sources = NULL
                    ) {

  url = paste0(file.path(gp_globals$base_url, "api", "gost", "profile"), "/")

  # Query

  if (is.null(query)) {
    stop("Missing query")
  } else if (is.list(query)) {
    if (is.data.frame(query)){
      stop("Query can't be a data.frame. Please use a vector or list of identifiers.")
    }
    # Multiple queries
    qnames = names(query)
    if (is.null(qnames)) {
      qnames = paste("query", seq(1, length(query)), sep = "_")
      names(query) = qnames
    }
    query = lapply(query, function(x) x[!is.na(x)])
  }
  else{
    query = query[!is.na(query)]
  }

  # Parameters

  ## evaluate choices
  correction_method <- match.arg(correction_method)

  if (startsWith(organism, "gp__")){
    message("Detected custom GMT source request")
    if (!is.null(sources)){
      message("Sources selection is not available for custom GMT requests. All sources in the GMT upload will be used.")
      sources <- NULL
    }
  }

  if (multi_query & evcodes){
    message("Evidence codes are not supported with multi_query option and will not be included in the results.\nYou can get evidence codes and intersection genes by setting multi_query = FALSE while keeping the input query as a list of multiple gene vectors.")
  }

  if (!is.null(custom_bg)){
    if (!is.vector(custom_bg)){
      stop("custom_bg must be a vector")
    }
    message("Detected custom background input, domain scope is set to 'custom'")
    domain_scope <- "custom"
    t <- ifelse(length(custom_bg) == 1, custom_bg <- jsonlite::unbox(custom_bg), custom_bg <- custom_bg)
  }

  domain_scope <- match.arg(domain_scope)

  body <- jsonlite::toJSON((
    list(
      organism = jsonlite::unbox(organism),
      query = query,
      sources = sources,
      user_threshold = jsonlite::unbox(user_threshold),
      all_results = jsonlite::unbox(!significant),
      ordered = jsonlite::unbox(ordered_query),
      no_evidences = jsonlite::unbox(!evcodes),
      combined = jsonlite::unbox(multi_query),
      measure_underrepresentation = jsonlite::unbox(measure_underrepresentation),
      no_iea = jsonlite::unbox(exclude_iea),
      domain_scope = jsonlite::unbox(domain_scope),
      numeric_ns = jsonlite::unbox(numeric_ns),
      significance_threshold_method = jsonlite::unbox(correction_method),
      background = custom_bg,
      output = jsonlite::unbox("json")
    )
  ),
  auto_unbox = FALSE,
  null = "null")
  # Headers

  headers <-
    list("Accept" = "application/json",
         "Content-Type" = "application/json",
         "charset" = "UTF-8")

  oldw <- getOption("warn")
  options(warn = -1)
  h1 = RCurl::basicTextGatherer(.mapUnicode = FALSE)
  h2 = RCurl::getCurlHandle() # Get info about the request

  # Request
  r = RCurl::curlPerform(
    url = url,
    postfields = body,
    httpheader = headers,
    customrequest = 'POST',
    verbose = FALSE,
    ssl.verifypeer = FALSE,
    writefunction = h1$update,
    curl = h2,
    .opts = gp_globals$rcurl_opts
  )
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2)[["response.code"]]
  txt <- h1$value()

  if (rescode != 200) {
    stop("Bad request, response code ", rescode)
  }

  res <- jsonlite::fromJSON(txt)
  df = res$result
  meta = res$meta

  if (is.null(dim(df))) {
    message("No results to show\n", "Please make sure that the organism is correct or set significant = FALSE")
    return(NULL)
  }

  # Re-order the data frame columns
  if (multi_query) {
    # add column that shows if significant
    df$significant <- lapply(df$p_values, function(x) x <= user_threshold)

    col_names <- c(
      "term_id",
      "p_values",
      "significant",
      "term_size",
      "query_sizes",
      "intersection_sizes",
      "source",
      "term_name",
      "effective_domain_size",
      "source_order",
      "parents"
    )

  } else {
    col_names <- c(
      "query",
      "significant",
      "p_value",
      "term_size",
      "query_size",
      "intersection_size",
      "precision",
      "recall",
      "term_id",
      "source",
      "term_name",
      "effective_domain_size",
      "source_order",
      "parents"
    )

    if (evcodes) {
      col_names <- append(col_names, c("evidence_codes", "intersection"))
      # Add column that lists the intersecting genes separated by comma
      df$intersection <- mapply(
        function(evcodes, query){
          genemap <- meta$genes_metadata$query[[query]]$mapping
          genes <- meta$genes_metadata$query[[query]]$ensgs[which(lengths(evcodes) > 0)]
          genes2 <- lapply(genes, function(x) ifelse(x %in% genemap, names(which(genemap == x)) , x))
          return(paste0(genes2, collapse = ","))
          },
        df$intersections,
        df$query,
        SIMPLIFY = TRUE
      )
      df$evidence_codes <- sapply(df$intersections, function(x)
          paste0(sapply(x[which(lengths(x) > 0)], paste0, collapse = " "), collapse = ","), USE.NAMES = FALSE)
    }

    # Order by query, source and p_value
    df <- df[with(df, order(query, source, p_value)), ]

  }

  # Rename native to term_id
  colnames(df)[colnames(df) == "native"] <- "term_id"
  colnames(df)[colnames(df) == "name"] <- "term_name"

  # Fix the row indexing to start from 1
  row.names(df) <- NULL
  df <- df[, col_names]
  return(list("result" = df, "meta" = meta))
}

#' Manhattan plot of functional enrichment results.
#'
#' This function creates a Manhattan plot out of the results from gprofiler2::gost().
#' The plot is very similar to the one shown in the g:GOSt web tool.
#'
#' @param gostres named list from gost() function (with names 'result' and 'meta')
#' @param capped whether the -log10(p-values) would be capped if >= 16, just as in the web options.
#' @param interactive if enabled, returns interactive plot using 'plotly'. If disabled, static 'ggplot()' object is returned.
#' @param pal values mapped to relevant colors for data sources.
#' @return The output is either a plotly object (if interactive = TRUE) or a ggplot object (if interactive = FALSE).
#' @author Liis Kolberg <liis.kolberg@@ut.ee>
#' @examples
#'  gostres <- gost(c("Klf4", "Pax5", "Sox2", "Nanog"), organism = "mmusculus")
#'  gostplot(gostres)
#' @export
#' @importFrom plotly %>%
gostplot <- function(gostres, capped = TRUE, interactive = TRUE, pal = c("GO:MF" = "#dc3912",
                                                                         "GO:BP" = "#ff9900",
                                                                         "GO:CC" = "#109618",
                                                                         "KEGG" = "#dd4477",
                                                                         "REAC" = "#3366cc",
                                                                         "WP" = "#0099c6",
                                                                         "TF" = "#5574a6",
                                                                         "MIRNA" = "#22aa99",
                                                                         "HPA" = "#6633cc",
                                                                         "CORUM" = "#66aa00",
                                                                         "HP" = "#990099")
){
  # gostres is the GOSt response list (contains results and metadata)
  # This function will plot only the sources that were asked from the query

  if( is.null(pal) ){
    pal <- c("GO:MF" = "#dc3912",
             "GO:BP" = "#ff9900",
             "GO:CC" = "#109618",
             "KEGG" = "#dd4477",
             "REAC" = "#3366cc",
             "WP" = "#0099c6",
             "TF" = "#5574a6",
             "MIRNA" = "#22aa99",
             "HPA" = "#6633cc",
             "CORUM" = "#66aa00",
             "HP" = "#990099")
  }

  if (!("result" %in% names(gostres))) stop("Name 'result' not found from the input")
  if (!("meta" %in% names(gostres))) stop("Name 'meta' not found from the input")

  source_order <- logpval <- term_id <- opacity <- NULL
  term_size <- term_name <- p_value <- term_size_scaled <- NULL

  df <- gostres$result
  # Order of data sources comes metadata
  meta <- gostres$meta

  # make sure all the essential column names are there
  essential_names <- c("source_order", "term_size", "term_name", "term_id", "source", "significant")

  if (!(all(essential_names %in% colnames(df)))) stop(paste("The following columns are missing from the result:", paste0(setdiff(essential_names, colnames(df)), collapse = ", ")))

  if (!any(grepl("p_value", colnames(df)))) stop("Column 'p_value(s)' is missing from the result")

  # nr_of_terms for every source
  widthscale <- unlist(lapply(meta$query_metadata$sources, function(x) meta$result_metadata[[x]][["number_of_terms"]]))
  names(widthscale) <- meta$query_metadata$sources # all the sources that were queried

  # Define the start positions for sources in the plot

  # start position for every term
  space <- 1000 # space between different sources
  starts <- c()
  start <- 1
  starts[1] <- start

  for(idx in 2:length(widthscale)){
    starts[idx] <- starts[idx - 1] + space + widthscale[idx - 1]
  }
  names(starts) <- names(widthscale)

  # Make sure that all the sources have colors

  if (is.null(names(pal))){
    names(pal) = meta$query_metadata$sources[1:length(pal)]
  }

  sourcediff = setdiff(meta$query_metadata$sources, names(pal))
  colors = grDevices::colors(distinct = TRUE)[grep('gr(a|e)y|white|snow|khaki|lightyellow', grDevices::colors(distinct = TRUE), invert = TRUE)]

  if (length(sourcediff) > 0){
    use_cols = sample(colors, length(sourcediff), replace = FALSE)
    pal[sourcediff] <- use_cols
  }

  # If multiquery
  if("p_values" %in% colnames(df)){
    p_values <- query <- significant <- NULL
    # spread the data frame to correct form
    df$query <- list(names(meta$query_metadata$queries))
    df <- tidyr::unnest(data = df, p_values, query, significant)
    df <- dplyr::rename(df, p_value = p_values)
  }

  # Set sizescale of circles
  logScale <- function(input, input_start = 1, input_end = 50000, output_start = 2, output_end = 10){
    m = (output_end - output_start)/(log(input_end) - log(input_start))
    b = -m*log(input_start) + output_start
    output = m * log(input) + b
    return(output)
  }

  # Scale x positions
  xScale <- function(input, input_start = 1, input_end = sum(widthscale) + (length(widthscale) - 1)*space, output_start = 2, output_end = 200){
    m = (output_end - output_start)/(input_end - input_start)
    b = -m*input_start + output_start
    output = m * input + b
    return(output)
  }

  # Add values to df needed for plotting
  # add -log10 pval
  df$logpval <- -log10(df$p_value)
  df$opacity <- ifelse(df$significant, 0.8, ifelse(df$p_value == 1, 0, 0.3))
  df$term_size_scaled = logScale(df$term_size)
  # add x axis position
  df <- df %>% dplyr::group_by(source) %>% dplyr::mutate(order = xScale(source_order, input_start = 1, input_end = widthscale[source], output_start = starts[source], output_end = starts[source] + widthscale[source]))
  df$order <- xScale(df$order)

  if (capped) {
    df$logpval[df$logpval > 16] <- 17
    ymin <- -1
    ymax <- 18.5
    ticklabels <- c("0", "2", "4", "6", "8", "10", "12", "14", ">16")
    tickvals <- c(0, 2, 4, 6, 8, 10, 12, 14, 16)
  } else {
    ymin <- -1
    ymax <- ceiling(max(df$logpval)) + 5
    ticklabels <- ggplot2::waiver()
    tickvals <- ggplot2::waiver()
  }

  if (interactive){
    # Start plotting
    sd <- crosstalk::SharedData$new(df, key = ~term_id)
  } else {
    sd <- df
  }

  p <- ggplot2::ggplot(data = sd, ggplot2::aes(x = order, y = logpval, text = paste(term_id, paste0('(', term_size,')'), '<br>', term_name, '<br>', formatC(p_value, format = "e", digits = 3)))) +
    ggplot2::geom_point(ggplot2::aes(color = source, size = term_size_scaled, alpha = opacity),
                        show.legend = FALSE) +
    ggplot2::facet_wrap(~ query, ncol = 1, scales = "free_x", shrink = FALSE) +
    ggplot2::ylab("-log10(p-adj)") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position='none',
                   panel.border=ggplot2::element_blank(),
                   strip.text=ggplot2::element_text(size=12, colour = "darkgrey"),
                   strip.background=ggplot2::element_blank(),
                   axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_text(size = 8, angle=45, hjust = 1),
                   axis.ticks.x=ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_line(color = "grey", size = 0.5),
                   axis.line.x = ggplot2::element_line(color = "grey", size = 0.1),
                   axis.line.y = ggplot2::element_line(size = 0.5, color = "grey"),
                   strip.text.x = ggplot2::element_text(angle = 0, hjust = 0, size = 10),
                   plot.margin = ggplot2::margin(t = 0, r = 5, b = 20, l = 20, unit = "pt"),
                   axis.title.y = ggplot2::element_text(size = 10, margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0))
    ) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::scale_alpha(range = c(0, 0.8), limits = c(0, 0.8)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(ymin, ymax),
                                labels = ticklabels,
                                breaks = tickvals) +
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 210), breaks = (xScale(starts) + xScale(starts + widthscale))/2, labels = names(widthscale))

  for (s in names(widthscale)){
    xstart = xScale(starts[s])
    xend = xScale(starts[s] + widthscale[s])
    p <- p + ggplot2::annotate("segment", x = xstart, xend = xend, y = -1, yend = -1,
                               size = 3, colour = pal[s])
  }

  if (capped){
    p <- p + ggplot2::annotate(geom = "text", x = 180,
                               y = 16.2, label = "values above this threshold are capped", size = 2, color = "grey") +
      ggplot2::geom_hline(yintercept = 16, linetype = "dashed", size = 0.2, color = "grey")
  }

  if (interactive){
    p <- p %>% plotly::ggplotly(tooltip = "text")
    p <- p %>% plotly::highlight(on = "plotly_click", off = "plotly_doubleclick", dynamic = FALSE, persistent = FALSE)
  }

  return(p)
}

#' Map vector of numeric values to Viridis color scale.
#'
#' @param values vector of numeric values (mostly -log10(p-values))
#' @param domain_min numeric value that corresponds to the 'yellow' in the color scale
#' @param domain_max numeric value that corresponds to the 'dark blue' in the color scale
#' @param n number of bins to generate from the color scale
#' @return The output is a corresponding vector of colors from the Viridis color scale with domain in range(domain_min, domain_max).
#' @author Liis Kolberg <liis.kolberg@@ut.ee>
#' @export
mapViridis <- function(values, domain_min = 0, domain_max = 50, n = 256){
  cols <- viridisLite::viridis(n, direction = -1)
  domain_range <- seq(from = domain_min, to = domain_max, length.out = n)
  val_colors <- unlist(lapply(values, function(x) ifelse(!is.na(x), cols[which.min(abs(x - domain_range))], "white")))
  return(val_colors)
}


#' Create and save an annotated Manhattan plot of enrichment results.
#'
#' This function allows to highlight a list of selected terms on the Manhattan plot created with the gprofiler2::gostplot() function.
#' The resulting plot is saved to a publication ready image if 'filename' is specified.
#' The plot is very similar to the one shown in the g:GOSt web tool after clicking on circles.
#'
#' @param p ggplot object from gostplot(gostres, interactive = FALSE) function
#' @param highlight_terms vector of selected term IDs from the analysis or a (subset) data.frame that has a column called 'term_id'. No annotation is added if set to NULL.
#' @param filename file name to create on disk and save the annotated plot. Filename extension should be from c("png", "pdf", "jpeg", "tiff", "bmp").
#' @param width plot width in inches. If not supplied, the size of current graphics device is used.
#' @param height plot height in inches. If not supplied, the size of current graphics device is used.
#' @return The output is a ggplot object.
#' @author Liis Kolberg <liis.kolberg@@ut.ee>
#' @examples
#'  gostres <- gost(c("Klf4", "Pax5", "Sox2", "Nanog"), organism = "mmusculus")
#'  p <- gostplot(gostres, interactive = FALSE)
#'  publish_gostplot(p, highlight_terms = c("GO:0001010", "WP:WP1763"))
#' @export
publish_gostplot <- function(p, highlight_terms = NULL, filename = NULL, width = NA, height = NA){
  # check that it is a static plot
  if(!("ggplot" %in% class(p))){
    warning("Highlighting terms in a Manhattan plot is available for a ggplot object only.\nPlease set 'interactive = F' in the gostplot() function and try again.")
    return(NULL)
  }

  term_id <- logpval <- term_size_scaled <- id <- query <- p_value <- NULL

  # add highlights
  if (!is.null(highlight_terms)) {
    if (is.data.frame(highlight_terms)){
      message("The input 'highlight_terms' is a data.frame and therefore the column 'term_id' will be used for detection.")
      if ("term_id" %in% colnames(highlight_terms)){
        highlight_terms <- highlight_terms$term_id
      }
      else{
        stop("No column named 'term_id'.")
      }
    }
    df <- p$data
    subdf <- base::subset(df, term_id %in% highlight_terms)

    if (nrow(subdf) == 0){
      message("None of the term IDs in the 'highlight_terms' was found from the results.")
      return(p)
    }

    subdf$id <- match(subdf$term_id, highlight_terms)

    p <- p + ggplot2::geom_point(data = subdf, ggplot2::aes(x = order, y = logpval, size = term_size_scaled), pch=21, colour = "black")
    p <- p + ggplot2::geom_text(data = subdf, size = 4, colour = "white", ggplot2::aes(label = as.character(id), family = "mono", fontface = "bold"), hjust = -1.2, vjust = -0.05) +
      ggplot2::geom_text(data = subdf, size = 4, colour = "black", fontface = "bold", ggplot2::aes(label = as.character(id)), hjust = -1.2, vjust = -0.05)

    # add table
    showdf <- subdf[,c("id", "source", "term_id", "term_name", "p_value", "query")]
    showdf$p_value <- formatC(showdf$p_value, format = "e", digits = 3)
    showdf <- tidyr::spread(showdf, query, p_value)

    idx <- 5:ncol(showdf)
    colours <- matrix("white", nrow(showdf), ncol(showdf))
    cols <- sapply(showdf[,idx], function(x) mapViridis(-log10(as.numeric(x))))
    colours[,idx] <- cols

    fontcolours <- matrix("black", nrow(showdf), ncol(showdf))
    fontcolours[,idx] <- "white"

    fontfaces <- matrix("plain", nrow(showdf), ncol(showdf))
    fontfaces[,idx] <- "bold"

    showdf[is.na(showdf)] <- ""
    th <- gridExtra::ttheme_default(base_size = 10,
                                    core=list(
                                      padding.h = grid::unit(c(3,3), "mm"),
                                      padding.v = grid::unit(c(2,2), "mm"),
                                      bg_params = list(fill = colours, col="black", lwd = 0.5),
                                      fg_params=list(hjust = 0, x = 0.01, col=fontcolours, fontface=fontfaces)),
                                    colhead=list(bg_params = list(fill = "gray99", lwd = 0.5, col = "black"),
                                                 fg_params=list(col="gray39", fontface="bold")),
                                    rowhead=list(fg_params=list(col="black", fontface="bold")))
    tb <- gridExtra::tableGrob(showdf, theme = th, rows = NULL)
    h <- grid::unit.c(grid::unit(1, "null"), sum(tb$heights) + grid::unit(3, "mm"))
    #w <- grid::unit.c(sum(tb$widths))
    w <- grid::unit.c(grid::unit(1, "null"))
    tg <- gridExtra::grid.arrange(p, tb, ncol = 1, heights = h, widths = w, newpage = TRUE, bottom = grid::textGrob("g:Profiler (biit.cs.ut.ee/gprofiler)", x = 0.95, hjust = 1, gp = grid::gpar(fontsize=10, font=8, col = "cornflowerblue")))
    # convert grob to ggplot object
    p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + ggplot2::geom_blank() + ggplot2::theme_void()
    #width = grid::convertWidth(sum(tg$widths), "in", TRUE) + 0.1
    #height = grid::convertHeight(sum(tg$heights), "in", TRUE) + 10.2
   }

  if (is.null(filename)){
    return(p)
  }
  else{
    imgtype <- strsplit(basename(filename), split="\\.")[[1]][-1]

    if (length(imgtype) == 0) {
      filename = paste0(filename, ".pdf")
    }

    if (tolower(imgtype) %in% c("png", "pdf", "jpeg", "tiff", "bmp")){
      if (is.na(width)){
        width = max(grDevices::dev.size()[1], 8)
      }
      if (is.na(height)){
        height = max(grDevices::dev.size()[2], 6)
      }
      ggplot2::ggsave(filename = filename, plot = p, width = width, height = height, limitsize = F)
      message("The image is saved to ", filename)
      return(p)
    } else {
      stop("The given file format is not supported.\nPlease use one of the following extensions: .png, .pdf, .jpeg, .tiff, .bmp")
    }
  }
}

#' Create and save a table with the functional enrichment analysis results.
#'
#' This function creates a table mainly for the results from gost() function.
#' However, if the input at least contains columns named 'term_id' and 'p_value' then any enrichment results data frame can be visualised in a table with this function.
#'
#' The output table is very similar to the one shown under the Manhattan plot.
#'
#' @param gostres named list from gost() function (with names 'result' and 'meta') or a data frame that has columns named "term_id" and "p_value(s)".
#' @param highlight_terms vector of selected term IDs from the analysis or a (subset) data.frame that has a column called 'term_id'. All data is shown if set to NULL.
#' @param use_colors if enabled, the p-values are highlighted in the viridis colorscale just as in g:Profiler, otherwise the table has no background colors.
#' @param show_columns names of additional columns to show besides term_id and p_value. By default the output table shows additional columns named "source", "term_name", "term_size", "intersection_size"
#' @param filename file name to create on disk and save the annotated plot. Filename extension should be from c("png", "pdf", "jpeg", "tiff", "bmp").
#' @return The output is a ggplot object.
#' @author Liis Kolberg <liis.kolberg@@ut.ee>
#' @examples
#'  gostres <- gost(c("Klf4", "Pax5", "Sox2", "Nanog"), organism = "mmusculus")
#'  publish_gosttable(gostres, highlight_terms = c("GO:0001010", "WP:WP1763"))
#' @export
publish_gosttable <- function(gostres, highlight_terms = NULL, use_colors = TRUE, show_columns = c("source", "term_name", "term_size", "intersection_size"), filename = NULL){
  # gostres is the GOSt response list (contains results and metadata) or a data frame

  term_id <- p_values <- query <- p_value <- NULL

  if (class(gostres) == "list"){
    if (!("result" %in% names(gostres))) stop("Name 'result' not found from the input")
    df <- gostres$result
  } else if (class(gostres) == "data.frame"){
    df <- gostres
  } else {
    stop("The input 'gostres' should be a data frame or a list from the gost() function.")
  }

  # make sure all the essential column names are there
  if (!"term_id" %in% colnames(df)) stop("The required column 'term_id' is missing from the input")
  if (!any(grepl("p_value", colnames(df)))) stop("Column 'p_value(s)' is missing from the input")

  # selected terms
  if (is.null(highlight_terms)) {
    # show full table if no terms given
    highlight_terms = df
  }

  if (is.data.frame(highlight_terms)){
    message("The input 'highlight_terms' is a data.frame. The column 'term_id' will be used.")
    if ("term_id" %in% colnames(highlight_terms)){
      highlight_terms <- highlight_terms$term_id
    }
    else{
      stop("No column named 'term_id'.")
    }
  }

  subdf <- base::subset(df, term_id %in% highlight_terms)

  if (nrow(subdf) == 0){
    stop("None of the term IDs in the 'highlight_terms' were found from the results.")
  }

  subdf$id <- match(subdf$term_id, highlight_terms)

  # default column names to show
  show_columns <- unique(append(show_columns, c("id", "term_id", "p_value")))
  gp_colnames <- c("id", "source", "term_id", "term_name", "term_size", "query_size", "intersection_size", "p_value")

  colnames <- gp_colnames[which(gp_colnames %in% show_columns)]

  # include non gprofiler columns
  if (length(setdiff(show_columns, gp_colnames)) > 0){
    colnames <- append(colnames, setdiff(show_columns, gp_colnames))
  }

  # If multiquery
  if ("p_values" %in% colnames(subdf)){
    if ("meta" %in% names(gostres)){
      meta <- gostres$meta
      subdf$query <- list(names(meta$query_metadata$queries))
    } else {
      qnames = paste("query", seq(1, length(subdf$p_values[[1]])), sep = "_")
      subdf$query <- list(names(qnames))
    }
    # spread the data frame to correct form
    subdf <- tidyr::unnest(data = subdf, p_values, query)
    subdf <- dplyr::rename(subdf, p_value = p_values)
    subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 3)
    showdf <- subdf[,stats::na.omit(match(c(colnames, "query"), names(subdf)))]
    showdf <- tidyr::spread(showdf, query, p_value)
    idx <- which(!is.na(match(names(showdf), unique(subdf$query))))
  } else {
    if ("query" %in% names(subdf) & length(unique(subdf$query)) > 1){
      subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 3)
      showdf <- subdf[,stats::na.omit(match(c(colnames, "query"), names(subdf)))]
      showdf <- tidyr::spread(showdf, query, p_value)
      idx <- which(!is.na(match(names(showdf), unique(subdf$query))))
    } else {
      subdf$p_value <- formatC(subdf$p_value, format = "e", digits = 3)
      showdf <- subdf[,stats::na.omit(match(colnames, names(subdf)))]
      idx <- which(!is.na(match(names(showdf), "p_value")))
    }
  }

  # Prepare table
  colours <- matrix("white", nrow(showdf), ncol(showdf))
  if (use_colors){
    cols <- sapply(showdf[,idx], function(x) mapViridis(-log10(as.numeric(x))))
  } else {
    cols <- sapply(showdf[,idx], function(x) "white")
  }

  colours[,idx] <- cols

  fontcolours <- matrix("black", nrow(showdf), ncol(showdf))
  fontcolours[,idx] <- ifelse(use_colors, "white", "black")

  fontfaces <- matrix("plain", nrow(showdf), ncol(showdf))
  fontfaces[,idx] <- "bold"

  showdf[is.na(showdf)] <- ""
  th <- gridExtra::ttheme_default(base_size = 10,
                                  padding = grid::unit(c(4, 4), "mm"),
                                  core=list(
                                    padding.h = grid::unit(c(15,15), "mm"),
                                    padding.v = grid::unit(c(15,15), "mm"),
                                    bg_params = list(fill = colours, col="black", lwd = 0.5),
                                    fg_params=list(hjust = 0, x = 0.01, col=fontcolours, fontface=fontfaces)),
                                  colhead=list(bg_params = list(fill = "gray99", lwd = 0.5, col = "black"),
                                               fg_params=list(col="gray39", fontface="bold")),
                                  rowhead=list(fg_params=list(col="black", fontface="bold")))

  tb <- gridExtra::tableGrob(showdf, theme = th, rows = NULL)
  h <- grid::unit.c(sum(tb$heights))
  w <- grid::unit.c(sum(tb$widths))
  tg <- gridExtra::grid.arrange(tb, ncol = 1, widths = w, heights = h, newpage = TRUE, bottom = grid::textGrob("g:Profiler (biit.cs.ut.ee/gprofiler)", x = 0.95, hjust = 1, gp = grid::gpar(fontsize=10, font=8, col = "cornflowerblue")))

  p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + ggplot2::geom_blank() + ggplot2::theme_void()

  if (is.null(filename)){
    return(p)
  } else{
    imgtype <- strsplit(basename(filename), split="\\.")[[1]][-1]

    if (length(imgtype) == 0) {
      filename = paste0(filename, ".pdf")
    }

    if (tolower(imgtype) %in% c("png", "pdf", "jpeg", "tiff", "bmp")){
      width = grid::convertWidth(sum(tg$widths), "in", TRUE) + 0.2
      height = grid::convertHeight(sum(tg$heights), "in", TRUE) + 0.2
      ggplot2::ggsave(filename = filename, plot = p, height = height, width = width)
      message("The image is saved to ", filename)
      return(p)
    } else {
      stop("The given file format is not supported.\nPlease use one of the following extensions: .png, .pdf, .jpeg, .tiff, .bmp")
    }
  }
}

#' Upload custom annotations for functional enrichment analysis in g:GOSt.
#'
#' Upload your own annotation data using files in the Gene Matrix Transposed file format (GMT) for functional enrichment analysis in g:GOSt.
#' The accepted file is either a single annotations file (with the extension .gmt) or a compressed directory of multiple annotation GMT files (with the extension .zip).
#' The GMT format is a tab-separated list of gene annotation sets where every line represents a separate gene set/functional term. The first column defines the function ID, second defines a short name/description of the function and the following columns
#' are the list of genes related to the specific function in that row.
#'
#' The uploaded filename is used to define 'source' name in the g:GOSt results.
#'
#' @param gmtfile the filepath of the GMT file to be uploaded. The file extension should be .gmt or .zip in case of multiple GMT files.
#' If the filepath does not contain an absolute path, the filename is relative to the current working directory.
#' @return A string that denotes the ID of the uploaded custom annotations in the g:Profiler database.
#' After the GMT file upload this unique ID can be used as a value for the argument 'organism' in the \code{gost()} function to perform
#' functional enrichment analysis based on these custom data.
#'
#' No need to repeatedly upload the same custom GMT file(s) every time you want to do the enrichment analysis.
#' The custom ID can also be used in the web tool as a token under the Custom GMT options.
#' @author Liis Kolberg <liis.kolberg@@ut.ee>
#' @examples
#' \dontrun{custom_id <- upload_GMT_file("path/to/file.gmt")}
#' @export
upload_GMT_file <- function(gmtfile){
  # check if file exists
  if (!file.exists(gmtfile)){
    stop(paste0("Can't find the input file ", gmtfile, "\nPlease check the filename and absolute path and try again."))
  }
  if (endsWith(gmtfile, ".gmt")){
    # GMT file
    gmturl = paste0(file.path(gp_globals$base_url, "api", "gost", "custom"), "/")
    gmtdata <- paste(readLines(gmtfile, skipNul = TRUE), collapse = "\n")

    body <- jsonlite::toJSON((
      list(
        gmt = gmtdata,
        name = basename(gmtfile)
      )
    ),
    auto_unbox = TRUE,
    null = "null")

    # Headers
    headers <-
      list("Accept" = "application/json",
           "Content-Type" = "application/json",
           "charset" = "UTF-8")

    oldw <- getOption("warn")
    options(warn = -1)
    h1 = RCurl::basicTextGatherer(.mapUnicode = FALSE)
    h2 = RCurl::getCurlHandle() # Get info about the request

    # Request
    r = RCurl::curlPerform(
      url = gmturl,
      postfields = body,
      httpheader = headers,
      customrequest = 'POST',
      verbose = FALSE,
      ssl.verifypeer = FALSE,
      writefunction = h1$update,
      curl = h2,
      .opts = gp_globals$rcurl_opts
    )
    options(warn = 0)
    rescode = RCurl::getCurlInfo(h2)[["response.code"]]

    if (rescode != 200) {
      stop("Bad request, response code ", rescode)
    }

    txt <- h1$value()
    res <- jsonlite::fromJSON(txt)
  }
  else if (endsWith(gmtfile, ".zip")){
    # zipped GMT files

    gmturl = paste0(file.path(gp_globals$base_url, "api", "gost", "custom", "zip"), "/")
    options(warn = -1)
    h2 = RCurl::getCurlHandle() # Get info about the request
    r = RCurl::postForm(
      uri = gmturl,
      zipfile = RCurl::fileUpload(filename = gmtfile, contentType="multipart/form-data"),
      curl = h2
    )
    options(warn = 0)
    rescode = RCurl::getCurlInfo(h2)[["response.code"]]

    if (rescode != 200) {
      stop("Bad request, response code ", rescode)
    }
    res <- jsonlite::fromJSON(r)
  }
  else {
    stop("Custom GMT file should have extension .gmt or .zip")
  }
  custom_id <- res$organism
  message(paste("Your custom annotations ID is", custom_id), "\nYou can use this ID as an 'organism' name in all the related enrichment tests against this custom source.")
  message(paste0("Just use: gost(my_genes, organism = '", custom_id, "')"))

  return(custom_id)
}


#' Generate a random gene list.
#'
#' This function returns a vector of randomly selected genes from the selected organism.
#'
#' @param organism organism name. Organism names are constructed by concatenating the first letter of the name and the
#' family name. Example: human - 'hsapiens', mouse - 'mmusculus'.
#' @return a character vector containing randomly selected gene IDs from the selected organism.
#' @author Liis Kolberg <liis.kolberg@@ut.ee>
#' @examples
#' random_genes <- random_query()
#' @export
random_query <- function(organism = "hsapiens"){
  url = paste0(file.path(gp_globals$base_url, "api", "gost", "random_query"), "/")
  body <- jsonlite::toJSON((
    list(organism = jsonlite::unbox(organism))),
    auto_unbox = FALSE,
    null = "null")
  # Headers

  headers <- list("Accept" = "application/json",
                  "Content-Type" = "application/json",
                  "charset" = "UTF-8")

  oldw <- getOption("warn")
  options(warn = -1)
  h1 = RCurl::basicTextGatherer(.mapUnicode = FALSE)
  h2 = RCurl::getCurlHandle() # Get info about the request

  # Request
  r = RCurl::curlPerform(
    url = url,
    postfields = body,
    httpheader = headers,
    customrequest = 'POST',
    verbose = FALSE,
    ssl.verifypeer = FALSE,
    writefunction = h1$update,
    curl = h2,
    .opts = gp_globals$rcurl_opts
  )
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2)[["response.code"]]
  txt <- h1$value()

  if (rescode != 200) {
    stop("Bad request, response code ", rescode)
  }
  res <- jsonlite::fromJSON(txt)
  return(res)
}

#' Gene ID conversion.
#'
#' Interface to the g:Profiler tool g:Convert (\url{https://biit.cs.ut.ee/gprofiler/convert}) that uses the information in Ensembl databases to handle hundreds of types of identifiers for genes, proteins, transcripts, microarray probesets, etc, for many species,
#' experimental platforms and biological databases.
#' The input is flexible: it accepts a mixed list of IDs and recognises their types automatically.
#' It can also serve as a service to get all genes belonging to a particular functional category.
#'
#' @param query vector that can consist of mixed types of gene IDs (proteins, transcripts, microarray IDs, etc), SNP IDs, chromosomal intervals or term IDs.
#' @param organism organism name. Organism names are constructed by
#' concatenating the first letter of the name and the family name. Example: human
#' - 'hsapiens', mouse - 'mmusculus'.
#' @param target target namespace.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param mthreshold maximum number of results per initial alias to show. Shows all by default.
#' @param filter_na logical indicating whether to filter out results without a
#' corresponding target.
#' @return The output is a data.frame which is a table closely corresponding to the
#' web interface output.
#' @author  Liis Kolberg <liis.kolberg@@ut.ee>, Uku Raudvere <uku.raudvere@@ut.ee>
#' @examples
#' gconvert(c("POU5F1", "SOX2", "NANOG"), organism = "hsapiens", target="AFFY_HG_U133_PLUS_2")
#' @export
gconvert = function(
  query,
  organism = "hsapiens",
  target = "ENSG",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE
) {
  url = file.path(gp_globals$base_url, "api", "convert", "convert")

  if (is.null(query)) {
    stop("Missing query")
  }

  if (is.list(query)) {
    stop("Parameter 'query' must be a vector of values")
  }

  body <- jsonlite::toJSON((
    list(
      organism = organism,
      query = query,
      target = target,
      numeric_ns = numeric_ns,
      output = "json"
    )
  ),
  auto_unbox = TRUE,
  null = "null")

  # Headers

  headers <-
    list("Accept" = "application/json",
         "Content-Type" = "application/json",
         "charset" = "UTF-8")

  oldw <- getOption("warn")
  options(warn = -1)
  h1 = RCurl::basicTextGatherer()
  h2 = RCurl::getCurlHandle() # Get info about the request

  # Request

  r = RCurl::curlPerform(
    url = url,
    postfields = body,
    httpheader = headers,
    customrequest = 'POST',
    verbose = FALSE,
    ssl.verifypeer = FALSE,
    writefunction = h1$update,
    curl = h2,
    .opts = gp_globals$rcurl_opts
  )
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2)[["response.code"]]
  txt <- h1$value()

  if (rescode != 200) {
    stop("Bad request, response code ", rescode)
  }

  res <- jsonlite::fromJSON(txt)
  df = res$result

  if (is.null(dim(df))) {
    message("No results to show\n", "Please make sure that the organism or namespace is correct")
    return(NULL)
  }

  col_names <- c("n_incoming", "incoming", "n_converted", "converted", "name", "description", "namespaces")
  df = df[,col_names]

  # rename columns
  colnames(df) <- c("input_number", "input", "target_number", "target", "name", "description", "namespace")
  df$target_number <- paste(df$input_number, df$target_number, sep=".")

  df = plyr::dlply(df, "input", function(x) {
    if (filter_na)
      x = x[x$target != "None",]
    if (nrow(x) == 0)
      return(NULL)
    if (nrow(x) > mthreshold)
      x = lapply(x, function(y) { as.character(y)[1:mthreshold] })
    return(x)
  })

  df <- plyr::ldply(df, function(x) data.frame(x, stringsAsFactors = F))

  return(df)
}


#' Orthology search.
#'
#' Interface to the g:Profiler tool g:Orth (\url{https://biit.cs.ut.ee/gprofiler/orth}) that, given a target organism, retrieves the genes of the target organism that are similar in sequence to the source organism genes in the input.
#'
#' @param query vector of gene IDs to be translated.
#' @param source_organism name of the source organism. Organism names are constructed by concatenating
#' the first letter of the name and the family name. Example: human - 'hsapiens',
#' mouse - 'mmusculus'.
#' @param target_organism name of the target organism. Organism names are constructed by concatenating
#' the first letter of the name and the family name. Example: human - 'hsapiens',
#' mouse - 'mmusculus'.
#' @param numeric_ns namespace to use for fully numeric IDs.
#' @param mthreshold maximum number of ortholog names per gene to show.
#' @param filter_na logical indicating whether to filter out results without a
#' corresponding target name.
#' @return The output is a data.frame which is a table closely corresponding to the
#' web interface output.
#' @author  Liis Kolberg <liis.kolberg@@ut.ee>, Uku Raudvere <uku.raudvere@@ut.ee>
#' @examples
#' gorth(c("Klf4","Pax5","Sox2","Nanog"), source_organism="mmusculus", target_organism="hsapiens")
#' @export

gorth <- function(
  query,
  source_organism = "hsapiens",
  target_organism = "mmusculus",
  numeric_ns = "",
  mthreshold = Inf,
  filter_na = TRUE
){
  url = file.path(gp_globals$base_url, "api", "orth", "orth")

  if (is.null(query)) {
    stop("Missing query")
  }

  if (is.list(query)) {
    stop("Parameter 'query' should be a vector of values")
  }

  body <- jsonlite::toJSON((
    list(
      organism = source_organism,
      target = target_organism,
      query = query,
      numeric_ns = numeric_ns,
      output = "json"
    )
  ),
  auto_unbox = TRUE,
  null = "null")

  # Headers

  headers <-
    list("Accept" = "application/json",
         "Content-Type" = "application/json",
         "charset" = "UTF-8")

  oldw <- getOption("warn")
  options(warn = -1)
  h1 = RCurl::basicTextGatherer()
  h2 = RCurl::getCurlHandle() # Get info about the request

  # Request

  r = RCurl::curlPerform(
    url = url,
    postfields = body,
    httpheader = headers,
    customrequest = 'POST',
    verbose = FALSE,
    ssl.verifypeer = FALSE,
    writefunction = h1$update,
    curl = h2,
    .opts = gp_globals$rcurl_opts
  )
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2)[["response.code"]]
  txt <- h1$value()

  if (rescode != 200) {
    stop("Bad request, response code ", rescode)
  }

  res <- jsonlite::fromJSON(txt)
  df = res$result

  if (is.null(dim(df))) {
    message("No results to show\n", "Please make sure that the organisms or the namespace are correct")
    return(NULL)
  }

  df$ensg_number = paste(df$n_incoming, df$n_converted, df$n_result, sep = ".")
  col_names = c("n_incoming", "incoming", "converted", "ensg_number", "name", "ortholog_ensg", "description")
  df <- df[,col_names]

  # rename columns
  colnames(df) <- c("input_number", "input", "input_ensg", "ensg_number", "ortholog_name", "ortholog_ensg", "description")

  df = plyr::dlply(df, "input", function(x) {
    if (filter_na)
      x = x[x$ortholog_name != "None",]
    if (nrow(x) == 0)
      return(NULL)
    if (nrow(x) > mthreshold)
      x = lapply(x, function(y) { as.character(y)[1:mthreshold] })
    return(x)
  })

  df <- plyr::ldply(df, function(x) data.frame(x, stringsAsFactors = F))

  return(df)
}


#' Convert SNP rs numbers to genes.
#'
#' Interface to the g:Profiler tool g:SNPense (\url{https://biit.cs.ut.ee/gprofiler/snpense}) that maps SNP rs identifiers to chromosome positions, genes and variant effects.
#' Available only for human SNPs.
#'
#' @param query vector of SNP IDs to be translated (should start with prefix 'rs').
#' @param filter_na logical indicating whether to filter out results without a
#' corresponding target name.
#' @return The output is a data.frame which is a table closely corresponding to the
#' web interface output. Columns 'ensgs' and 'gene_names' can contain list of multiple values.
#' @author  Liis Kolberg <liis.kolberg@@ut.ee>, Uku Raudvere <uku.raudvere@@ut.ee>
#' @examples
#' gsnpense(c("rs11734132", "rs7961894", "rs4305276", "rs17396340", "rs3184504"))
#' @export

gsnpense <- function(
  query, filter_na = TRUE
){
  url = file.path(gp_globals$base_url, "api", "snpense", "snpense")

  if (is.null(query)) {
    stop("Missing query")
  }

  if ( !any(startsWith(tolower(query), "rs")) ){
    stop("Query must contain SNP ids with 'rs' prefix.")
  }

  body <- jsonlite::toJSON((
    list(
      query = query,
      output = "json"
    )
  ),
  auto_unbox = TRUE,
  null = "null")

  # Headers

  headers <-
    list("Accept" = "application/json",
         "Content-Type" = "application/json",
         "charset" = "UTF-8")

  oldw <- getOption("warn")
  options(warn = -1)
  h1 = RCurl::basicTextGatherer()
  h2 = RCurl::getCurlHandle() # Get info about the request

  # Request

  r = RCurl::curlPerform(
    url = url,
    postfields = body,
    httpheader = headers,
    customrequest = 'POST',
    verbose = FALSE,
    ssl.verifypeer = FALSE,
    writefunction = h1$update,
    curl = h2,
    .opts = gp_globals$rcurl_opts
  )
  options(warn = 0)
  rescode = RCurl::getCurlInfo(h2)[["response.code"]]
  txt <- h1$value()

  if (rescode != 200) {
    stop("Bad request, response code ", rescode)
  }

  res <- jsonlite::fromJSON(txt)
  df = res$result

  if (is.null(dim(df))) {
    message("No results to show\n", "Please make sure that the input is of correct format")
    return(NULL)
  }

  # Add NA where relevant

  df[df$start == -1, -1] = NA

  if (filter_na) {
    df <- df[!is.na(df$chromosome),]
  }
  row.names(df) <- NULL
  return(df)
}


#' Get current user agent string.
#'
#' Get the HTTP User-Agent string.
#'
#' @export

get_user_agent = function() {
  gp_globals$rcurl_opts$useragent
}

#' Set custom user agent string.
#'
#' Set the HTTP User-Agent string. Useful for overriding the default user agent for
#' packages that depend on gprofiler2 functionality.
#'
#' @param ua the user agent string.
#' @param append logical indicating whether to append the passed string to the default user agent string.
#' @export

set_user_agent = function(ua, append = F) {
  rco = gp_globals$rcurl_opts
  rco$useragent = ifelse(
    append, paste(rco$useragent, ua, sep=""), ua)
  assign("rcurl_opts", rco, envir=gp_globals)
}

#' Get the TLS version for SSL
#'
#' @export

get_tls_version = function() {
  v = gp_globals$rcurl_opts$sslversion
  if (v == gp_globals$CURL_SSLVERSION_TLSv1_2) {
    "1.2"
  }
  else if (v == gp_globals$CURL_SSLVERSION_TLSv1_1) {
    "1.1"
  }
}

#' Set the TLS version to use for SSL
#'
#' Set the TLS version. Could be useful at environments where SSL was built without TLS 1.2 support
#'
#' @param v version: "1.2" (default), "1.1" (fallback)
#' @export

set_tls_version = function(v) {
  if (v == "1.1") {
    rco = gp_globals$rcurl_opts
    rco$sslversion = gp_globals$CURL_SSLVERSION_TLSv1_1
    assign("rcurl_opts", rco, envir=gp_globals)
  }
  else if (v == "1.2") {
    rco = gp_globals$rcurl_opts
    rco$sslversion = gp_globals$CURL_SSLVERSION_TLSv1_2
    assign("rcurl_opts", rco, envir=gp_globals)
  }
  else {
    print('Only "1.1" (fallback) or "1.2" (default) are allowed to be passed as a parameter.')
  }
}


#' Get the current base URL.
#'
#' @export

get_base_url = function() {
  gp_globals$base_url
}

#' Set the base URL.
#'
#' Function to change the g:Profiler base URL. Useful for overriding the default URL (http://biit.cs.ut.ee/gprofiler)
#' with the beta (http://biit.cs.ut.ee/gprofiler_beta) or an archived version (available starting from the version e94_eg41_p11, e.g. http://biit.cs.ut.ee/gprofiler_archive3/e94_eg41_p11).
#'
#' @param url the base URL.
#' @export

set_base_url = function(url) {
  url = as.character(url)
  schema = substr(url, 1, 4)

  if (schema != "http")
    stop("The URL must be absolute and use the HTTP(S) schema")

  assign("base_url", url, envir = gp_globals)
}

