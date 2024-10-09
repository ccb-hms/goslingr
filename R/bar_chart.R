#' Title
#'
#' @param gr 
#' @param chromosome_field 
#' @param genomic_fields 
#' @param width 
#' @param height 
#' @param y_field 
#' @param bin_size 
#' @param title 
#' @param subtitle 
#'
#' @return
#' @export
#'
#' @examples
#' 
#' # generate example GRanges object
#' library(GenomicRanges)
#' set.seed(42)
#'
#' # number of rows per chromosome
#' n <- 1000  
#' 
#' gr <- GRanges(
#'   seqnames = rep(paste0("chr", c(1:22, 'X', 'Y')), each = n),
#'   ranges = IRanges(
#'     start = rep(sample(1:250000000, n, replace = TRUE), each = 24),
#'     width = round(600 + rexp(n, rate = 1/1800)),
#'     peak = round(50 + rexp(n, rate = 1/500))
#'   )
#' )
#' 
#' # generate gosling spec and render visualization
#' spec <- bar_chart_spec(gr)
#' print(spec)
#' 
#' goslingr(spec)

bar_chart_spec <- function(
    gr, 
    chromosome_field = "seqnames", 
    genomic_fields = c("start", "end"),
    width = 800,
    height = 180,
    y_field = "peak",
    size = 5,
    title = "Basic Marks: bar",
    subtitle = "Tutorial Examples"
) {
  
  # Extract data from GRanges object
  values <- as.data.frame(gr)
  
  # construct gosling data spec
  gr_data <- list(
    type = "json",
    chromosomeField = chromosome_field,
    genomicFields = genomic_fields,
    values = values
  )
  
  spec <- GoslingSpec(list(
    title = title,
    subtitle = subtitle,
    tracks = list(
      list(
        layout = "linear",
        width = width,
        height = height,
        data = gr_data,
        mark = "bar",
        x = list(field = genomic_fields[1], type = "genomic", axis = "bottom"),
        xe = list(field = genomic_fields[2], type = "genomic"),
        y = list(field = y_field, type = "quantitative", axis = "right"),
        size = list(value = size)
      )
    )
  ))
  
  return(spec)
}
