
# goslingr

<!-- badges: start -->
<!-- badges: end -->

The goal of goslingr is to ...

## Installation

You can install the development version of goslingr like so:

``` r
# install.packages('remotes')
remotes::install_github('ccb-hms/goslingr')

# install gosling python package (one time)
library(goslingr)
install_gosling()

```

## Define a Gosling Specification

`gosling` visualizations are defined using a Gosling specification, which is an 
R list or JSON object. You can define your spec directly as an R list:


``` r
library(goslingr)
library(GenomicRanges)

# toy GRanges object
gr <- GRanges(
  seqnames = paste0('chr', c(1:22, 'X', 'Y')),
  peak = rnorm(24, mean = 500, sd = 100),
  ranges = IRanges(
    start = rep(0, 24),
    end = rep(10, 24)
  )
)

# data part of spec
gr_data_spec <- list(
  type = "json",
  chromosomeField = "seqnames",
  genomicFields = c("start", "end"),
  values = as.data.frame(gr)
)

# example bar chart spec
spec <- list(
  title = "Basic Marks: bar",
  subtitle = "Tutorial Examples",
  tracks = list(
    list(
      layout = "linear",
      width = 800,
      height = 180,
      data = gr_data_spec,
      mark = "bar",
      x = list(field = "start", type = "genomic", axis = "bottom"),
      xe = list(field = "end", type = "genomic"),
      y = list(field = "peak", type = "quantitative", axis = "right"),
      size = list(value = 5)
    )
  )
)

goslingr(spec)
```

You can alternatively define the specification using the python 
[gos](https://gosling-lang.github.io/gos/index.html) api and then extract
the JSON specification:

``` r
library(goslingr)

data <- gos$bigwig(
  url = "https://s3.amazonaws.com/gosling-lang.org/data/ExcitatoryNeurons-insertions_bin100_RIPnorm.bw",
  column = "position",
  value = "peak"
)

track <- gos$Track(data, height = 100)$mark_point()$encode(
  x = gos$X("position:G"),
  y = gos$Y("peak:Q")
)

view <- track$view()

# extract gosling specification
spec <- view$to_json()

# generate the plot
goslingr(spec)
```

Note that you need to substitute with R 
equivalents (`$` instead of `.`, `c()` instead of `[]`, `TRUE` instead of `True`, etc).

Finally, you can provide the specification as a JSON string (useful for getting
started from [gosling.js examples](http://gosling-lang.org/examples/)):

``` r
spec <- '{
  "title": "Basic Marks: bar",
  "subtitle": "Tutorial Examples",
  "tracks": [
    {
      "layout": "linear",
      "width": 800,
      "height": 180,
      "data": {
        "url": "https://resgen.io/api/v1/tileset_info/?d=UvVPeLHuRDiYA3qwFlm7xQ",
        "type": "multivec",
        "row": "sample",
        "column": "position",
        "value": "peak",
        "categories": ["sample 1"],
        "binSize": 5
      },
      "mark": "bar",
      "x": {"field": "start", "type": "genomic", "axis": "bottom"},
      "xe": {"field": "end", "type": "genomic"},
      "y": {"field": "peak", "type": "quantitative", "axis": "right"},
      "size": {"value": 5}
    }
  ]
}'

# generate the plot
goslingr(spec)
```

For full details on how to define a `gosling` specification check out the [docs](http://gosling-lang.org/docs/).

## Local Data

See [local data](https://gosling-lang.github.io/gos/user_guide/local_data.html) example.
`goslingr` will server data using `gos` if a local file path is detected. This
will work no matter how you define your specification.