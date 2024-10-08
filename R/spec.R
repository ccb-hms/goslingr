#' Constructor function for GoslingSpec
#' @export
GoslingSpec <- function(spec) {
  structure(spec, class = 'GoslingSpec')
}

#' Print method for GoslingSpec objects
#' @export
print.GoslingSpec <- function(x, ...) {
  cat("GoslingSpec(")
  print_list(x)
  cat(")")
}

#' Create gosling data spec from GRanges object
#'
#' @param gr GRanges object
#'
#' @return
#' @export
#'
#' @examples
create_gr_data <- function(gr) {
  list(
    type = "json",
    chromosomeField = "seqnames",
    genomicFields = c("start", "end"),
    values = as.data.frame(gr)
  )
}