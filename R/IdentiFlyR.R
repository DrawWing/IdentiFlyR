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

#' Create XML document with classification data
#'
#' @param inData a matrix-like R object, with at least six columns
#' @param grVec a vector with grouping variable
#' @param file a character string with file name
#' @param prototype an optional character string with prototype file name. It is used by IdentiFly to describe landmarks.
#'
#' @export
#'
#' @examples
#' data(lineages)
#' grVec = lineages$lineage
#' coordinates = lineages[,-1] # remove the first column
#' gmLdaData2xml(coordinates, grVec, "apis-mellifera-lineage.dw.xml")
#' idData = xml2gmLdaData("apis-mellifera-lineage.dw.xml")
#' id = gmLdaData2id(idData, coordinates, average = FALSE)
#' id$plot
gmLdaData2xml = function(inData, grVec, file, prototype = ""){

  # Error detection
  if((ncol(inData) %% 2) != 0)
    stop("Number of columns should be even.")
  if(ncol(inData) < 6)
    stop("Number of columns should be 6 or greater.")
  # if(length(grVec != nrow(inData)))
  #   stop("Length of inVector should be equal to number of rows of inData.")

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

  means = stats::aggregate(data, by = list(grVec), FUN = mean)
  rownames(means) = means$Group.1
  means[,"Group.1"] = NULL
  colnames(means) = xyNames
  matrix2XML(means, "means", XML)

  PCA <- stats::prcomp(data)
  data = as.matrix(data) %*% PCA$rotation
  data <- data[,1:(2*p-4)]

  n = length(unique(grVec)) # number of groups
  LDA = MASS::lda(data, grVec, prior=rep(1/n, n))

  scores = as.matrix(data) %*% LDA$scaling
  scores = as.data.frame(scores)

  groups = rownames(means)
  covariances = lapply(groups, function(x)stats::cov(scores[grVec==x,]))
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

  if(prototype!="")
    XML$addTag("prototype", close=TRUE, attrs = c(file=prototype))

  XML::saveXML(XML$value(), file=file,
               prefix = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
}

#' Read classification data from XML file
#'
#' @param path a file name with xml document
#'
#' @return list of classification data
#' \itemize{
#'   \item reference - matrix with reference configuration.
#'   \item means - matrix with class means.
#'   \item covariances - list of covariance matrices.
#'   \item coefficients - matrix with LDA coefficients.
#' }
#'
#' @export
#'
#' @examples
#' xmlPath = system.file("extdata",
#'                       "apis-mellifera-queens-workers-drones.dw.xml",
#'                       package="IdentiFlyR")
#' idData = xml2gmLdaData(xmlPath)
#' data(lineages)
#' unknownData <- lineages[,-1]
#' id = gmLdaData2id(idData, unknownData, average = FALSE)
#' id$plot
xml2gmLdaData = function(path){
  XML = xml2::read_xml(path)

  # extract lda
  ldaXML = xml2::xml_find_first(XML, ".//lda")
  matrices = xml2::xml_find_all(ldaXML, ".//matrix")

  # extract reference
  referenceXML = matrices[xml2::xml_attr(matrices, "id")=="reference"]
  referenceMat = xm2matrix(referenceXML)

  # extract means
  meansXML = matrices[xml2::xml_attr(matrices, "id")=="means"]
  meansMat = xm2matrix(meansXML)

  # extract covariances
  rowNames = rownames(meansMat)
  # create list of empty dataframes
  covariances = list(rep(matrix(), length(rowNames)))
  for (i in 1:length(rowNames)) {
    covXML = matrices[xml2::xml_attr(matrices, "id")==rowNames[i]]
    covMat = xm2matrix(covXML)
    covariances[[i]] = covMat
    names(covariances)[i] = rowNames[i]
  }

  # extract coefficients
  coefficientsXML = matrices[xml2::xml_attr(matrices, "id")=="coefficients"]
  coefficientsMat = xm2matrix(coefficientsXML)

  return(list("reference" = referenceMat,
              "means" = meansMat,
              "covariances" = covariances,
              "coefficients" = coefficientsMat))
}


#' Classify unknown samples using identification data
#'
#' @param idData a list of identification data
#' @param data a matrix of unknown data to classify
#' @param average an optionally Boolean variable indicating if rows are averaged
#'
#' @return list of id data
#' @export
#'
#' @examples
#' xmlPath = system.file("extdata",
#'                       "apis-mellifera-queens-workers-drones.dw.xml",
#'                       package="IdentiFlyR")
#' idData = xml2gmLdaData(xmlPath)
#' data(lineages)
#' unknownData <- lineages[,-1]
#' id <- gmLdaData2id(idData, unknownData)
gmLdaData2id = function(idData, data, average = TRUE){
  # calculate means LD by matrix multiplication
  meansLD = idData$means %*% t(idData$coefficients)

  if(ncol(data) != ncol(idData$reference))
    stop("number of columns differes between reference and data")
  p = ncol(data)/2
  k = 2

  # Convert from 2D to 3D array
  data.3D = geomorph::arrayspecs(data, p, k)

  # convert reference from 1D to 2D
  reference = matrix(idData$reference, ncol=2, byrow = TRUE)

  if(average) {
    # Align the coordinates using Generalized Procrustes Analysis
    GPA = geomorph::gpagen(data.3D, print.progress = FALSE)
    # Consensus configuration of the unknown sample
    unknownConsensus = GPA$consensus

    # Align unknown consensus with reference
    unknownOPA = shapes::procOPA(reference, unknownConsensus)
    unknownAligned = unknownOPA$Bhat

    # convert aligned data from 2D to 1D
    unknownAligned = matrix(t(unknownAligned), nrow=1)

    # calculate LD scores
    LdTab = unknownAligned %*% t(idData$coefficients)

    id = classifyVecLD(LdTab, meansLD, idData$covariances)
    id = id$class

    LdTab = as.data.frame(LdTab)
    plot = covEllipses(meansLD, idData$covariances) +
      ggplot2::geom_point(LdTab, mapping=ggplot2::aes(x = LdTab$LD1, y = LdTab$LD2, colour = "zzz") ) +
      ggrepel::geom_label_repel(LdTab,
                                mapping=ggplot2::aes(x = LdTab$LD1, y = LdTab$LD2, colour = "zzz", label = "unknown"),
                                nudge_x = 0.75, nudge_y = 0) +
      ggplot2::scale_color_manual(values = append(grDevices::rainbow(nrow(idData$means)), "black"))
  } else
  {
    # calculate lda scores for all wings
    LD.list = vector(mode = "list", length = nrow(data))
    for (r in 1:nrow(data)) {
      # Align unknown consensus with reference
      unknownOPA = shapes::procOPA(reference, data.3D[,,r])
      unknownAligned = matrix(t(unknownOPA$Bhat), nrow=1) # convert from 2D to 1D
      LD.row = unknownAligned %*% t(idData$coefficients)
      LD.list[[r]] = LD.row
    }
    LdTab = do.call(rbind, LD.list)
    rownames(LdTab) = rownames(data)
    LdTab = as.data.frame(LdTab)

    id = classifyMatLD(LdTab, meansLD, idData$covariances)

    plot = covEllipses(meansLD, idData$covariances) +
      ggplot2::geom_point(LdTab, mapping=ggplot2::aes(x = LdTab$LD1, y = LdTab$LD2, colour = "zzz") ) +
      ggplot2::scale_color_manual(values = append(grDevices::rainbow(nrow(idData$means)), "black"))
  }

  return(list("id" = id,
              "plot" = plot,
              "LDx" = LdTab))
}

# calculate Mahalnanobis distances for data in one vector
classifyVecLD = function(unknown.LD, means, covariances){
  groups = rownames(means)
  df = ncol(unknown.LD) - 1
  resultList = list() # empty list for results
  maxP = 0
  maxGroup = ""
  for (i in 1:length(groups)) {
    MD = stats::mahalanobis(unknown.LD, means[i, ], as.matrix(covariances[[i]]))
    results = c(MD)
    P = stats::pchisq(MD, df=df, lower.tail=FALSE)
    results = c(results, P)
    if(P > maxP){
      maxP = P
      maxGroup = groups[i]
    }
    resultList[[i]] = results
  }

  outSummary = paste("The sample was classified as", maxGroup, "with probability", maxP)
  outMax = data.frame("group" = maxGroup, "P" = maxP)

  outClass = do.call(rbind, resultList)
  colnames(outClass) = c("MD2", "P")
  rownames(outClass) = rownames(means)

  return(list("summary" = outSummary,
              "max.group" = outMax,
              "class" = outClass))
}

# calculate Mahalnanobis distances for each row in a matrix
classifyMatLD = function(unknown.LD, means, covariances){
  resultList = list() # empty list for results
  for (i in 1:nrow(unknown.LD)) {
    outRow = classifyVecLD(unknown.LD[i,], means, covariances)
    resultList[[i]] = outRow$max
  }
  idResults = do.call(rbind, resultList)
  rownames(idResults) = rownames(unknown.LD)

  return(idResults)
}

#' Classify unknown data using data from XML file
#'
#' @param path a file name with xml document
#' @param data a matrix of unknown data to classify
#' @param average an optionally Boolean variable indicating if rows are averaged
#'
#' @return list of classification data
#' \itemize{
#'   \item reference - matrix with reference configuration.
#'   \item means - matrix with class means.
#'   \item covariances - list of covariance matrices.
#'   \item coefficients - matrix with LDA coefficients.
#' }
#'
#' @export
#'
#' @examples
#' xmlPath = system.file("extdata",
#'                       "apis-mellifera-queens-workers-drones.dw.xml",
#'                       package="IdentiFlyR")
#' data(lineages)
#' unknownData <- lineages[,-1]
#' id <- xml2id(xmlPath, unknownData)
xml2id = function(path, data, average = TRUE){
  idData = xml2gmLdaData(path)
  return(gmLdaData2id(idData, data, average))
}


#' Plot ellipses using covariance matrices without access to raw data
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
#' xmlPath = system.file("extdata",
#'                       "apis-mellifera-queens-workers-drones.dw.xml",
#'                       package="IdentiFlyR")
#' idData = xml2gmLdaData(xmlPath)
#' idMeans = idData$means %*% t(idData$coefficients)
#' covEllipses(idMeans, idData$covariances)
covEllipses = function(means, covariances, x = 1, y = 2){
  chi = sqrt(stats::qchisq(0.95, 2)) # Chi-Square value for probability 95%

  groups = rownames(means)
  ellipseList = list()
  for (i in 1:length(groups)) {
    x0 = means[i, x]
    y0 = means[i, y]

    co = covariances[[i]]
    covMat = cbind(c(co[x,x], co[x,y]),
                   c(co[x,y], co[y,y]))
    # covMat = covMat[1:2, 1:2]
    covEigen = eigen(covMat)

    a = chi*sqrt(covEigen$values[1]);
    b = chi*sqrt(covEigen$values[2]);
    angle = atan2(covEigen$vectors[2,1], covEigen$vectors[1,1])

    ellipse = c(x0, y0, a, b, angle)
    ellipseList[[i]] = ellipse
  }
  ellipses = do.call(rbind, ellipseList)
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
    ggforce::geom_ellipse(ellipses, mapping = ggplot2::aes(x0 = !! ggplot2::ensym(LDx),
                                                           y0 = !! ggplot2::ensym(LDy),
                                                           a = a, b = b, angle = angle)) +
    ggrepel::geom_label_repel( ggplot2::aes(label = rownames(ellipses), size = NULL), nudge_y = 0.75) +
    ggplot2::theme(legend.position="none")
}

#' Reorder 2D landmarks using column names
#'
#' @param inData a matrix-like R object, with two dimensions
#' @param newOrder a numerical vector with new order of landmarks
#'
#' @return a matrix with reordered columns
#' @export
#'
#' @examples
#' xmlPath = system.file("extdata",
#'                       "apis-mellifera-queens-workers-drones.dw.xml",
#'                       package="IdentiFlyR")
#' idData = xml2gmLdaData(xmlPath)
#' p <- 19 # number of landmarks
#' k <- 2 # number of dimensions
#' reference = idData$reference
#' reference3D <- geomorph::arrayspecs(reference, p, k)
#' geomorph::plotAllSpecimens(reference3D, label=TRUE)
#' Bustamente2IdentiFly =
#' c(18, 17, 15, 14, 16, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 19)
#' reference = reorder2D(reference, Bustamente2IdentiFly)
#' reference3D <- geomorph::arrayspecs(reference, p, k)
#' geomorph::plotAllSpecimens(reference3D, label=TRUE)
reorder2D = function(inData, newOrder){

  if(!is.matrix(inData))
    inData = as.matrix(inData)
  if(!is.matrix(inData))
    stop("inData is not a matrix.")

  if ((ncol(inData)%%2) != 0)
    stop("Number of columns is not even.")
  p <- ncol(inData)/2 # number of landmarks
  k <- 2 # number of dimensions

  # create coordinates names used by IdentiFly
  xyNames <- c("x1", "y1")
  for (i in 2:p) {
    xyNames <- c(xyNames, paste0("x", i))
    xyNames <- c(xyNames, paste0("y", i))
  }

  # is newOrder a correct sequence
  sorted = newOrder[order(newOrder)]
  if(!all(sorted == seq(from = 1, to = p, by = 1)))
    stop("Incorrect format of newOrder vector.\n",
         "After sorting it should be equal to seq(from = 1, to = p, by = 1)\n",
         "where p is number of landmarks: ncol(inData)/2.")

  # are all column names present?
  if(!all(xyNames %in% colnames(inData)))
    stop("Column names of inData should be x1, y1, x2, y2, ...")
  # are column names unique?
  if(any(duplicated(colnames(inData))))
    stop("Columnis of inData shoudl have uique names.")

  xy = paste0("x", newOrder[1])
  xy = c(xy, paste0("y", newOrder[1]))
  reordered = matrix(inData[,xy], ncol=2)
  for (i in 2:length(newOrder)) {
    xy = paste0("x", newOrder[i])
    xy = c(xy, paste0("y", newOrder[i]))
    reordered = cbind(reordered, matrix(inData[,xy], ncol=2))
  }
  colnames(reordered) = xyNames
  rownames(reordered) = rownames(inData)

  return(reordered)
}
