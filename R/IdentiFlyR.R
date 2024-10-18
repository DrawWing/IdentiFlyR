# read matrix from XML
xm2matrix = function(inXML){
  # extract header
  header = xml2::xml_find_first(inXML, ".//header")
  headerTxt = xml2::xml_text(header)
  headerTxt = trimws(headerTxt)
  headerLst = strsplit(headerTxt, "\t")

  # extract rows
  rowsXML = xml2::xml_find_all(inXML, ".//vector")
  rowLst = list()
  rowNames = vector()
  for (i in 1:length(rowsXML)) {
    rowTxt = xml2::xml_text(rowsXML[[i]])
    rowTxt = trimws(rowTxt)
    rowVal = strsplit(rowTxt, "\t")
    rowVal = as.numeric(rowVal[[1]])
    rowLst[[i]] = rowVal

    # extract row attribute
    rowAtr = xml2::xml_attr(rowsXML[[i]], "id")
    rowAtr = trimws(rowAtr)
    rowNames = c(rowNames, rowAtr)
  }
  mat = do.call(rbind, rowLst)
  colnames(mat) = headerLst[[1]]
  rownames(mat) = rowNames
  return(mat)
}

# convert vector to XML
vector2XML = function(inVec, inId, outXML){
  inVec = as.character(inVec)
  inVec = paste(inVec, collapse = "\t")
  inVec = c("\t", inVec, "\t") # pad with tabulators for spreadsheet
  outXML$addTag("vector", inVec, attrs = c(id=inId))
}

# convert matrix to XML
matrix2XML = function(inMat, inId, outXML){
  outXML$addTag("matrix", close = FALSE, attrs = c(id=inId))
  header = as.character(colnames(inMat))
  header = paste(header, collapse = "\t")
  header = c("\t", header, "\t")
  outXML$addTag("header", header)

  for (i in 1:nrow(inMat)) {
    vector2XML(inMat[i,], rownames(inMat)[i], outXML)
  }
  outXML$closeTag() # matrix
}

#' Create DW.XML file with classification data
#'
#' @param data a matrix-like R object, with at least two dimensions
#' @param grouping a vector with grouping variable
#'
#' @return XML document
#' @export
#'
#' @examples
#' outXML <- gmLdaData2xml(inData, grVec)
gmLdaData2xml = function(inData, grouping){

  # Error detection
  if((ncol(inData) %% 2) != 0)
    stop("number of columns should be even")
  p = ncol(inData)/2
  k = 2

  # create coordinates names used by IdentiFly
  xyNames = c("x1", "y1")
  for (i in 2:p) {
    xyNames = c(xyNames, paste0("x", i))
    xyNames = c(xyNames, paste0("y", i))
  }

  ## Alignment of landmark configurations

  # Convert 2D array into a 3D array
  inData.3D = geomorph::arrayspecs(inData, p, k)
  # Align the coordinates using Generalized Procrustes Analysis
  GPA = geomorph::gpagen(inData.3D, print.progress = FALSE)
  # Convert 3D array into a 2D array - opposite to arrayspecs
  data = geomorph::two.d.array(GPA$coords)
  colnames(data) = xyNames

  XML = XML::xmlOutputDOM("identifly", attrs=c(version="1.8"))
  XML$addTag("lda", close=FALSE)

  reference = colMeans(data)
  reference = as.matrix(t(reference))

  colnames(reference) = xyNames
  rownames(reference) = c("reference")
  matrix2XML(reference, "reference", XML)

  means = aggregate(data, by = list(grouping), FUN = mean)
  rownames(means) = means$Group.1
  means = subset (means, select = -c(Group.1))
  colnames(means) = xyNames
  matrix2XML(means, "means", XML)

  PCA <- prcomp(data)
  data = as.matrix(data) %*% PCA$rotation
  data <- data[,1:(2*p-4)]

  n = length(unique(grouping)) # number of groups
  LDA = MASS::lda(data, grouping, prior=rep(1/n, n))

  scores = as.matrix(data) %*% LDA$scaling
  scores = as.data.frame(scores)

  groups = rownames(means)
  covariances = lapply(groups, function(x)cov(scores[grouping==x,]))
  XML$addTag("covariances", close=FALSE)
  for (i in 1:length(groups)) {
    matrix2XML(covariances[[i]], groups[i], XML)
  }
  XML$closeTag() # covariances

  # coefficients = t(LDA$scaling)
  coefficients = t(PCA$rotation[,1:(2*p-4)] %*% LDA$scaling)
  colnames(coefficients) = xyNames
  matrix2XML(coefficients, "coefficients", XML)

  XML$closeTag() # lda

  return(XML)
}

#' Read classification data from XML file
#'
#' @param path file name with xml document
#'
#' @return list of classification data
#' @export
#'
#' @examples
#' idDataList <- xml2gmLdaData("Apis-mellifera-lineages.dw.xml")
xml2gmLdaData = function(path){
  # add error checking
  XML = xml2::read_xml(path)

  # extract lda
  lda.XML = xml2::xml_find_first(XML, ".//lda")
  matrices = xml2::xml_find_all(lda.XML, ".//matrix")

  # extract reference
  reference.XML = matrices[xml2::xml_attr(matrices, "id")=="reference"]
  reference.mat = xm2matrix(reference.XML)

  # extract means
  means.XML = matrices[xml2::xml_attr(matrices, "id")=="means"]
  means.mat = xm2matrix(means.XML)

  # extract covariances
  rowNames = rownames(means.mat)
  # create list of empty dataframes
  covariances = list(rep(matrix(), length(rowNames)))
  for (i in 1:length(rowNames)) {
    cov.XML = matrices[xml2::xml_attr(matrices, "id")==rowNames[i]]
    cov.mat = xm2matrix(cov.XML)
    covariances[[i]] = cov.mat
    names(covariances)[i] = rowNames[i]
  }

  # extract coefficients
  coefficients.XML = matrices[xml2::xml_attr(matrices, "id")=="coefficients"]
  coefficients.mat = xm2matrix(coefficients.XML)

  return(list("reference" = reference.mat,
              "means" = means.mat,
              "covariances" = covariances,
              "coefficients" = coefficients.mat))
}


#' Classify unknown samples using XML data
#'
#' @param path file name with xml document
#' @param data unknown data to classify
#' @param average average the data or classify each row
#'
#' @return list of id data
#' @export
#'
#' @examples
#' idData <- xml2id("Apis-mellifera-lineages.dw.xml", unknownDat)
xml2id = function(path, data, average = TRUE){
  id.data = xml2gmLdaData(path)

  # calculate means LD by matrix multiplication
  means.LD = id.data$means %*% t(id.data$coefficients)

  if(ncol(data) != ncol(id.data$reference))
    stop("number of columns differes between reference and data")
  p = ncol(data)/2
  k = 2

  # Convert from 2D to 3D array
  data.3D = geomorph::arrayspecs(data, p, k)

  # convert reference from 1D to 2D
  reference = matrix(id.data$reference, ncol=2, byrow = TRUE)

  if(average) {
    # Align the coordinates using Generalized Procrustes Analysis
    GPA = geomorph::gpagen(data.3D, print.progress = FALSE)
    # Consensus configuration of the unknown sample
    unknown.consensus = GPA$consensus

    # Align unknown consensus with reference
    unknown.OPA = shapes::procOPA(reference, unknown.consensus)
    unknown.aligned = unknown.OPA$Bhat

    # convert aligned data from 2D to 1D
    unknown.aligned = matrix(t(unknown.aligned), nrow=1)

    # calculate LD scores
    LD.tab = unknown.aligned %*% t(id.data$coefficients)

    id = classifyVecLD(LD.tab, means.LD, id.data$covariances)
    id = id$class

    LD.tab = as.data.frame(LD.tab)
    plot = covEllipses(means.LD, id.data$covariances) +
      ggplot2::geom_point(LD.tab, mapping=ggplot2::aes(x = LD1, y = LD2, colour = "zzz") ) +
      ggrepel::geom_label_repel(LD.tab,
                       mapping=ggplot2::aes(x = LD1, y = LD2, colour = "zzz", label = "unknown"),
                       nudge_x = 0.75, nudge_y = 0) +
      ggplot2::scale_color_manual(values = append(rainbow(nrow(id.data$means)), "black"))
  } else
  {
    # calculate lda scores for all wings
    LD.list = vector(mode = "list", length = nrow(data))
    for (r in 1:nrow(data)) {
      # Align unknown consensus with reference
      unknown.OPA = shapes::procOPA(reference, data.3D[,,r])
      unknown.aligned = matrix(t(unknown.OPA$Bhat), nrow=1) # convert from 2D to 1D
      LD.row = unknown.aligned %*% t(id.data$coefficients)
      LD.list[[r]] = LD.row
    }
    LD.tab = do.call(rbind, LD.list)
    rownames(LD.tab) = rownames(data)
    LD.tab = as.data.frame(LD.tab)

    id = classifyMatLD(LD.tab, means.LD, id.data$covariances)

    plot = covEllipses(means.LD, id.data$covariances) +
      ggplot2::geom_point(LD.tab, mapping=ggplot2::aes(x = LD1, y = LD2, colour = "zzz") ) +
      ggplot2::scale_color_manual(values = append(rainbow(nrow(id.data$means)), "black"))
  }

  return(list("id" = id,
              "plot" = plot,
              "LDx" = LD.tab))
}

# calculate Mahalnanobis distances for data in one vector
classifyVecLD = function(unknown.LD, means, covariances){
  groups = rownames(means)
  df = ncol(unknown.LD) - 1
  result.list = list() # empty list for results
  max.P = 0
  max.group = ""
  for (i in 1:length(groups)) {
    MD = mahalanobis(unknown.LD, means[i, ], as.matrix(covariances[[i]]))
    results = c(MD)
    P = pchisq(MD, df=df, lower.tail=FALSE)
    results = c(results, P)
    if(P > max.P){
      max.P = P
      max.group = groups[i]
    }
    result.list[[i]] = results
  }

  out.summary = paste("The sample was classified as", max.group, "with probability", max.P)
  out.max = data.frame("group" = max.group, "P" = max.P)

  out.class = do.call(rbind, result.list)
  colnames(out.class) = c("MD2", "P")
  rownames(out.class) = rownames(means)

  return(list("summary" = out.summary,
              "max.group" = out.max,
              "class" = out.class))
}

# calculate Mahalnanobis distances for each row in a matrix
classifyMatLD = function(unknown.LD, means, covariances){
  result.list = list() # empty list for results
  for (i in 1:nrow(unknown.LD)) {
    out.row = classifyVecLD(unknown.LD[i,], means, covariances)
    result.list[[i]] = out.row$max
  }
  id.results = do.call(rbind, result.list)
  rownames(id.results) = rownames(unknown.LD)

  return(id.results)
}

# plot ellipses for the two LD

#' Plot ellipses for the first two LD using data from covariance matrices
#'
#' @param means a matrix-like R object, with at least two dimensions
#' @param covariances a matrix-like R object, with at least two dimensions
#' @param x index of data for the x axis
#' @param y index of data for the y axis
#'
#' @return ggplot2 plot
#' @export
#'
#' @examples
#' idData = xml2gmLdaData("apis-mellifera-lineage.dw.xml")
#' idMeans = idData$means %*% t(idData$coefficients)
#' covEllipses(idMeans, idData$covariances)

covEllipses = function(means, covariances, x = 1, y = 2){
  chi = sqrt(qchisq(0.95, 2)) # Chi-Square value for probability 95%

  groups = rownames(means)
  ellipse.list = list()
  for (i in 1:length(groups)) {
    x0 = means[i, x]
    y0 = means[i, y]

    co = covariances[[i]]
    cov.mat = cbind(c(co[x,x], co[x,y]),
                    c(co[x,y], co[y,y]))
    # cov.mat = cov.mat[1:2, 1:2]
    cov.eigen = eigen(cov.mat)

    a = chi*sqrt(cov.eigen$values[1]);
    b = chi*sqrt(cov.eigen$values[2]);
    angle = atan2(cov.eigen$vectors[2,1], cov.eigen$vectors[1,1])

    ellipse = c(x0, y0, a, b, angle)
    ellipse.list[[i]] = ellipse
  }
  ellipses = do.call(rbind, ellipse.list)
  LDx = paste0("LD", x)
  LDy = paste0("LD", y)
  colnames(ellipses) = c(LDx, LDy, "a", "b", "angle")
  rownames(ellipses) = groups
  ellipses = as.data.frame(ellipses)

  ggplot2::ggplot(ellipses, ggplot2::aes(x = !! ggplot2::ensym(LDx),
                                         y = !! ggplot2::ensym(LDy),
                                         color = rownames(ellipses)) ) +
    ggplot2::geom_point() +
    ggplot2::coord_fixed() +
    ggforce::geom_ellipse(ellipses, mapping=ggplot2::aes(x0 = !! ggplot2::ensym(LDx),
                                                y0 = !! ggplot2::ensym(LDy),
                                                a = a, b = b, angle = angle)) +
    ggrepel::geom_label_repel( ggplot2::aes(label = rownames(ellipses), size = NULL), nudge_y = 0.75) +
    ggplot2::theme(legend.position="none")
}
