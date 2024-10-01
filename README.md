
# goslingr

<!-- badges: start -->
<!-- badges: end -->

The goal of goslingr is to ...

## Installation

You can install the development version of goslingr like so:

``` r
# install.packages('remotes')
remotes::install_github('ccb-hms/goslingr')
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(goslingr)

# import bar chart spec file
bar_chart_spec_fpath <- system.file('extdata', 'bar_chart_spec.json', package = 'goslingr')
bar_chart_spec <- jsonlite::read_json(bar_chart_spec_fpath)

# generate gosling bar plot from spec
goslingr(bar_chart_spec)

# import lollipops plot spec file
lollipop_spec_fpath <- system.file('extdata', 'lollipop_plots_spec.json', package = 'goslingr')
lollipop_spec <- jsonlite::read_json(lollipop_spec_fpath)

goslingr(lollipop_spec)
```

