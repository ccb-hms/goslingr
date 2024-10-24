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
#' # generate example GPos object
#' library(GenomicRanges)
#' set.seed(42)
#' 
#' # number of rows per chromosome
#' n <- 100
#'
#' gpos <- GPos(
#'   seqnames = Rle(rep(paste0("chr", 1:22), each = n)),
#'   pos = rep(sample(250e6, n), each = 22)
#' )
#'
#' npos <- length(gpos)
#' mcols(gpos)$counts <- sample(100, npos, replace=TRUE)
#' 
#' 
#' # generate gosling spec and render visualization
#' spec <- line_chart_spec(gpos)
#' print(spec)
#' 
#' goslingr(spec)

line_chart_spec <- function(
    gpos, 
    chromosome_field = "seqnames", 
    genomic_fields = "pos",
    width = 800,
    height = 180,
    y_field = "counts",
    size = 2,
    title = "Basic Marks: line",
    subtitle = "Tutorial Examples"
) {
  
  # Extract data from GPos object
  values <- as.data.frame(gpos)
  
  # construct gosling data spec
  gpos_data <- list(
    type = "json",
    chromosomeField = chromosome_field,
    genomicFields = list(genomic_fields),
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
        data = gpos_data,
        mark = "line",
        x = list(field = genomic_fields[1], type = "genomic", axis = "bottom"),
        y = list(field = y_field, type = "quantitative", axis = "right"),
        size = list(value = size)
      )
    )
  ))
  
  
  return(spec)
}
