library(GenomicRanges)
library(goslingr)

url <- "https://rb.gy/7y3fx"
temp_file <- file.path(tempdir(), "data.gz")
download.file(url, destfile = temp_file, method = "auto", mode = "wb")
df <- read.delim(
  temp_file,
  header = FALSE,
  comment.char = "#",
  sep = ""
)

gr <- GRanges(
  seqnames = df$V1,
  ranges = IRanges(df$V2, df$V3, peak = df$V5)
)

gr

# Extract data from GRanges object
values <- as.data.frame(gr)

# Create the desired JSON structure as a list
gr_data <- list(
  type = "json",
  chromosomeField = "seqnames",
  genomicFields = list("start", "end"),
  values = values
)

spec <- list(
  title = "Basic Marks: bar",
  subtitle = "Tutorial Examples",
  tracks = list(
    list(
      layout = "linear",
      width = 800,
      height = 180,
      data = gr_data,
      mark = "bar",
      x = list(field = "start", type = "genomic", axis = "bottom"),
      xe = list(field = "end", type = "genomic"),
      y = list(field = "peak", type = "quantitative", axis = "right"),
      size = list(value = 5)
    )
  )
)

print_list(spec)
goslingr(spec)


# simpler example GRanges object

library(GenomicRanges)

gr <- GRanges(
  seqnames = paste0('chr', c(1:22, 'X', 'Y')),
  peak = rnorm(24, mean = 500, sd = 100),
  ranges = IRanges(
    start = rep(0, 24),
    end = rep(10, 24)
  )
)

# create a bar chart spec
spec <- bar_chart_spec(gr)

# print it for copying/editing in R
print(spec)
goslingr(spec)
