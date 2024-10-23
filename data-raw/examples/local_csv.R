library(goslingr)

# get the file
fname <- 'UCSC.HG38.Human.CytoBandIdeogram.csv'
csv_url <- file.path('https://raw.githubusercontent.com/sehilyi/gemini-datasets/master/data', fname)

temp_file <- file.path(tempdir(), fname)
download.file(csv_url, temp_file)

spec <- list(
  tracks = list(
    list(
      width = 700,
      height = 70,
      data = list(
        # where the data is being served
        url = temp_file,
        chromosomeField = "Chromosome",
        type = "csv",
        genomicFields = list(
          "chromStart",
          "chromEnd"
        )
      ),
      mark = "rect",
      color = list(
        field = "Stain",
        type = "nominal",
        domain = list(
          "gneg",
          "gpos25",
          "gpos50",
          "gpos75",
          "gpos100",
          "gvar"
        ),
        range = list(
          "white",
          "#D9D9D9",
          "#979797",
          "#636363",
          "black",
          "#A0A0F2"
        )
      ),
      x = list(
        field = "chromStart",
        type = "genomic",
        domain = list(
          chromosome = "1"
        ),
        axis = "top"
      ),
      xe = list(
        field = "chromEnd",
        type = "genomic"
      ),
      size = list(
        value = 20
      ),
      stroke = list(
        value = "gray"
      ),
      strokeWidth = list(
        value = 0.5
      )
    )
  )
)

goslingr(spec)
