
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IdentiFlyR

<!-- badges: start -->
<!-- badges: end -->

IdentiFlyR can be used to store identification data in an XML file.
Later, the data can be read from the XML file and used to identify an
unknown sample. An alternative to this package could be to provide all
the raw data and use it for identification. For large datasets this
could be less convenient, slower and undesirable. Identification is
based on linear discriminant analysis (LDA), also known as canonical
variate analysis (CVA). The focus of this package is on using LDA to
identify an unknown sample, not on exploratory analysis.  
Currently IdentiFLyR only supports geometric morphometric data in two
dimensions. \## Installation

You can install the development version of IdentiFlyR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("DrawWing/IdentiFlyR")
```

## Examples

Read identification data from XML file and plot means and confidence
ellipses in the first two linear discriminant functions. Identification
data consists of: reference, means, covariances and coefficients.

``` r
library(IdentiFlyR)
idData = xml2gmLdaData("https://zenodo.org/record/10512712/files/apis-mellifera-queens-workers-drones.dw.xml")
names(idData)
#> [1] "reference"    "means"        "covariances"  "coefficients"
# transfer means to LDA space
meansLda = idData$means %*% t(idData$coefficients)
covEllipses(meansLda, idData$covariances)
```

<img src="man/figures/README-read-xml-1.png" width="100%" />

Read raw coordinated of 19 landmarks.

``` r
wings <- read.csv("https://zenodo.org/record/8071014/files/IN-raw-coordinates.csv")
wings <- data.frame(wings, row.names = 1)  # move column 1 to row names
```

Classify the mean of all data. The mean of all rows has been classified
as “workers”.

``` r
id = gmLdaData2id(idData, wings, average = TRUE)
id$plot
#> Warning: Use of `LdTab$LD1` is discouraged.
#> ℹ Use `LD1` instead.
#> Warning: Use of `LdTab$LD2` is discouraged.
#> ℹ Use `LD2` instead.
#> Warning: Use of `LdTab$LD1` is discouraged.
#> ℹ Use `LD1` instead.
#> Warning: Use of `LdTab$LD2` is discouraged.
#> ℹ Use `LD2` instead.
```

<img src="man/figures/README-classify-mean-1.png" width="100%" />

``` r
id$id
#>               MD2            P
#> drone  73.8759245 8.318471e-18
#> queen  64.6842675 8.791401e-16
#> worker  0.5041501 4.776823e-01
```

Classify rows. All 350 rows were classified as “workers”.

``` r
id = gmLdaData2id(idData, wings, average = FALSE)
id$plot
```

<img src="man/figures/README-classify-rows-1.png" width="100%" />

``` r
head(id$id, 2)
#>                          group         P
#> IN-0001-000243-L.dw.png worker 0.5746041
#> IN-0001-000243-R.dw.png worker 0.3949792
table(id$id$group)
#> 
#> worker 
#>    350
```

Classify the first row. It was classified as “worker”.

``` r
id = gmLdaData2id(idData, wings[1,])
id$id
#>               MD2            P
#> drone  69.7285902 6.805252e-17
#> queen  67.5410674 2.063440e-16
#> worker  0.3150395 5.746041e-01
```

Create identification data and store it in an XML file.

``` r
wingsLin = read.csv("https://zenodo.org/record/7567336/files/Nawrocka_et_al2018-sample-aligned.csv")
grVec = wingsLin$lineage
wingsLin = wingsLin[, -c(1:4)] # remove unwanted columns
XML = gmLdaData2xml(wingsLin, grVec)
XML$addTag("prototype", close=TRUE, attrs = c(file="apis-worker-prototype.dw.png")) # add prototype for IdentiFly software
library(XML) 
XML::saveXML(XML$value(), file="apis-mellifera-lineage.dw.xml",
             prefix = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
#> [1] "apis-mellifera-lineage.dw.xml"
```
