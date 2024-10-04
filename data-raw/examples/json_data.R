reference_data <- list(
  type = "json",
  chromosomeField = "c",
  genomicFields = list("p"),
  values = list(
    list(
      c = "chr2",
      p = 100000,
      v = 0.0001
    ),
    list(
      c = "chr5",
      p = 100000,
      v = 0.0004
    ),
    list(
      c = "chr10",
      p = 100000,
      v = 0.0009
    )
  )
)

spec <- list(
  title = "Rule Mark",
  subtitle = "Annotate visualization with horizontal and vertical lines",
  style = list(dashed = c(3, 3)),
  tracks = list(
    list(
      data = reference_data,
      mark = "rule",
      x = list(field = "p", type = "genomic"),
      y = list(field = "v", type = "quantitative", domain = c(0, 0.003)),
      strokeWidth = list(field = "v", type = "quantitative"),
      color = list(value = "red")
    )
  )
)

goslingr(spec)
