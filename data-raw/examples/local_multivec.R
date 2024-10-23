# example from gos: https://gosling-lang.github.io/gos/user_guide/local_data.html
library(goslingr)


# download local multivec file (4GB)
# https://s3.amazonaws.com/gosling-lang.org/data/cistrome.multires.mv5
local_multivec_fpath <- './path/to/cistrome.multires.mv5'


data <- list(
  type='multivec',
  url=local_multivec_fpath,
  row="sample",
  column="position",
  value="peak",
  categories=c("sample 1", "sample 2", "sample 3", "sample 4"),
  binSize=4
)

spec <- list(
  tracks = list(
    list(
      color = list(
        field = "peak",
        legend = TRUE,
        type = "quantitative"
      ),
      data = data,
      height = 100,
      mark = "rect",
      row = list(
        field = "sample",
        legend = TRUE,
        type = "nominal"
      ),
      width = 725,
      x = list(
        axis = "top",
        field = "start",
        type = "genomic"
      ),
      xe = list(
        field = "end",
        type = "genomic"
      )
    )
  )
)

goslingr(spec)
