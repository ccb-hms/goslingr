#' Print a list as source code
#'
#' This function prints a nested list as formatted R source code. 
#' It's particularly useful when you have a nested R list and want to 
#' generate the R code to recreate or modify it.
#'
#' @param x A list.
#' @param indent Integer. The level of indentation (used for recursive calls). Default is 0.
#'
#' @return No return value, called for side effects (printing).
#' @export
#'
#' @examples
#' 
#' # example GRanges object
#' 
#' library(GenomicRanges)
#' 
#' gr <- GRanges(
#'   seqnames = paste0('chr', c(1:22, 'X', 'Y')),
#'   peak = rnorm(24, mean = 500, sd = 100),
#'   ranges = IRanges(
#'     start = rep(0, 24),
#'     end = rep(10, 24)
#'   )
#' )
#' 
#' # create a bar chart spec
#' spec <- bar_chart_spec(gr)
#' 
#' # print it for copying/editing in R
#' print(spec)
#' 
print_list <- function(x, indent = 0) {
  indent_str <- strrep(" ", indent * 2)
  
  # cat with no sep
  cat0 <- function(...) cat(..., sep = "")
  
  # coerce factors to character
  if (is.factor(x)) x <- as.character(x)
  
  # Handle lists
  if (is.list(x)) {
    list_type <- 'list'
    
    # handle data frames
    if (is.data.frame(x)) {
      list_type <- 'data.frame'
      x <- head(x, 5)
    }
    
    if (length(x) == 0) {
      cat0(list_type, "()")
    } else {
      cat0(list_type, "(\n")
      for (i in seq_along(x)) {
        # Print truncation comment if its a data.frame
        if (is.data.frame(x[[i]])) {
          cat0(indent_str, "  ", "# NOTE: truncated to first 5 rows of data.frame\n")
        }
        
        # Print name if it exists
        if (!is.null(names(x)[i])) {
          cat0(indent_str, "  ", names(x)[i], " = ")
        } else {
          cat0(indent_str, "  ")
        }
        # Recursively print the list item
        print_list(x[[i]], indent + 1)
        if (i < length(x)) cat(",")
        cat("\n")
      }
      cat0(indent_str, ")")
    }
  } 
  # Handle vectors (of length > 1)
  else if (is.vector(x) && length(x) > 1) {
    if (is.character(x)) {
      # For character vectors, wrap each element in quotes
      cat0('c(', paste(sprintf('"%s"', x), collapse = ", "), ')')
    } else {
      cat0('c(', paste(format(x), collapse = ", "), ')')
    }
  } 
  # Handle single character values
  else if (is.character(x)) {
    cat0('"', x, '"')
  } 
  # Handle single numeric values
  else if (is.numeric(x)) {
    cat(x)
  } 
  # Handle any other types
  else {
    cat(x)
  }
}