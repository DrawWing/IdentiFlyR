lineages <- read.csv('https://zenodo.org/record/7567336/files/Nawrocka_et_al2018-sample-aligned.csv')
lineages = lineages[,-1] #remove first column with row numbers
lineages <- data.frame(lineages, row.names = 1)  # move column 1 to row names
lineages = lineages[,-1] #remove first column with row numbers

usethis::use_data(lineages, overwrite = TRUE)
