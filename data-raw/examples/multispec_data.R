spec <- list(
  title = "Basic Marks: bar",
  subtitle = "Tutorial Examples",
  tracks = list(
    list(
      layout = "linear",
      width = 800,
      height = 180,
      mark = "bar",
      x = list(field = "start", type = "genomic", axis = "bottom"),
      xe = list(field = "end", type = "genomic"),
      y = list(field = "peak", type = "quantitative", axis = "right"),
      size = list(value = 5),
      data = list(
        url = "https://resgen.io/api/v1/tileset_info/?d=UvVPeLHuRDiYA3qwFlm7xQ",
        type = "multivec",
        row = "sample",
        column = "position",
        value = "peak",
        categories = list("sample 1"),
        binSize = 5
      )
    )
  )
)

goslingr::goslingr(spec)