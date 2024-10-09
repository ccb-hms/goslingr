library(goslingr)


# example 1
vis <- 
  gos$Track(
    gos$bigwig("https://s3.amazonaws.com/gosling-lang.org/data/HFFc6_H3K4me3.bigWig")
  )$encode(
    x='position:G',
    y='value:Q'
  )$view()


json_spec <- vis$to_json()
goslingr(json_spec)

# example 2
data <- gos$multivec(
  url="https://server.gosling-lang.org/api/v1/tileset_info/?d=cistrome-multivec",
  row="sample",
  column="position",
  value="peak",
  categories=c("sample 1", "sample 2", "sample 3", "sample 4"),
  binSize=5,
)

base_track <- gos$Track(data, width=800, height=100)

heatmap <- base_track$mark_rect()$encode(
  x=gos$X("start:G", axis="top"),
  xe="end:G",
  row=gos$Row("sample:N", legend=TRUE),
  color=gos$Color("peak:Q", legend=TRUE),
)

bars <- base_track$mark_bar()$encode(
  x=gos$X("position:G", axis="top"),
  y="peak:Q",
  row="sample:N",
  color=gos$Color("sample:N", legend=TRUE),
)

lines <- base_track$mark_line()$encode(
  x=gos$X("position:G", axis="top"),
  y="peak:Q",
  row="sample:N",
  color=gos$Color("sample:N", legend=TRUE),
)

view <- gos$vertical(heatmap, bars, lines)$properties(
  title="Visual Encoding",
  subtitle="Gosling provides diverse visual encoding methods",
  layout="linear",
  centerRadius=0.8,
  xDomain=gos$GenomicDomain(chromosome="1", interval=c(1, 3000500)),
)


json_spec <- view$to_json()
goslingr(json_spec)


# example 3 ----

# generate example GRanges object
library(GenomicRanges)
set.seed(42)

# number of rows per chromosome
n <- 1000

gr <- GRanges(
  seqnames = rep(paste0("chr", c(1:22, 'X', 'Y')), each = n),
  ranges = IRanges(
    start = rep(sample(1:250000000, n, replace = TRUE), each = 24),
    width = round(600 + rexp(n, rate = 1/1800)),
    peak = round(50 + rexp(n, rate = 1/500))
  )
)

# convert gr factor columns to character
gr <- as.data.frame(gr)
factor_cols <- sapply(gr, is.factor)
gr[factor_cols] <- lapply(gr[factor_cols], as.character)

# convert gr to json-like R structure
values <- purrr::transpose(gr)

data <- gos$JsonData(
  type = "json",
  chromosomeField =  "seqnames",
  genomicFields = c("start", "end"),
  values = values
)

bars <- gos$Track(
  data,
  width=800, 
  height=100
)$mark_bar()$encode(
  x=gos$X("start:G", axis="bottom"),
  xe="end:G",
  y=gos$Y("peak:Q", axis="right"),
  size=gos$SizeValue(5)
)

view <- bars$view(title="Basic Marks: Bar", subtitle="Tutorial Examples")

json_spec <- view$to_json()
goslingr(json_spec)
