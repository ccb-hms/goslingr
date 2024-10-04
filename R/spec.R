#' Add a gosling view
#'
#' @param list_obj 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
create_view <- function(...) {
  new_view <- list(...)
  return(new_view)
}


#' Add a gosling track
#'
#' @param list_obj 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
create_track <- function(...) {
  new_track <- list(...)
  return(new_track)
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