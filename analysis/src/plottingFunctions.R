plotCDF <- function(qps, add=FALSE, ...) {
  
  def.props <- qps$prop*as.numeric(gsub("[^0-9\\.]", "", names(qps$qps)))/100
  if(add) {
    lines(qps$qps, def.props, ...)
  } else {
    plot(qps$qps, def.props, ...)
  }
}