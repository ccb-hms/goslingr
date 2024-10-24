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
#' gpos <- GPos(c("chr1:1-20000000", "chr1:5-10", "chr2:2-5", "chr3:500-605", "chr4:700-800"))
#' npos <- length(gpos)
#' mcols(gpos)$counts <- sample(10, npos, replace=TRUE)
#' 
#' temp_file <- file.path(tempdir(), 'gr_data.csv')
#' data.table::fwrite(as.data.frame(gpos), temp_file)
#' 
#' # generate gosling spec and render visualization
#' spec <- line_chart_spec('http://127.0.0.1:8000/gr_data.csv',  y_field = 'counts')
#' print(spec)
#' 
#' goslingr(spec)
# servr::daemon_stop(3)

line_chart_spec <- function(
    gpos, 
    chromosome_field = "seqnames", 
    genomic_fields = "pos",
    width = 800,
    height = 180,
    y_field = "peak",
    size = 2,
    title = "Basic Marks: line",
    subtitle = "Tutorial Examples"
) {
  
  # Extract data from GPos object
  # values <- as.data.frame(gpos)
  
  # construct gosling data spec
  gpos_data <- list(
    type = "csv",
    url = gpos,
    sampleLength = 200000,
    chromosomeField = chromosome_field,
    genomicFields = list(genomic_fields)
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
