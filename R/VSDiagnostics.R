### CITE papers for the VSD and F/G statistics

### Performance Assessment of High-Dimensional Variable Identification
FFunc = function(trueActiveSet,selectedSet) {
  p1 = 2*length(intersect(trueActiveSet,selectedSet))
  p2 = length(trueActiveSet) + length(selectedSet)
  res = p1/p2
  return(res)
}

GFunc = function(trueActiveSet,selectedSet) {
  p1 = length(intersect(trueActiveSet,selectedSet))
  p2 = sqrt(length(trueActiveSet)*length(selectedSet))
  res = p1 / p2
  return(res)
}

### Variable Selection Diagnostics Measures for High-Dimensional Regression

# VSD+: active variables that were not selected
VSDPlusFunc = function(trueActiveSet,selectedSet) {
  res = length(setdiff(trueActiveSet,selectedSet))
  return(res)
}

# VSD-: how many inactive variables were selected
VSDMinusFunc = function(trueActiveSet,selectedSet) {
  res = length(setdiff(selectedSet,trueActiveSet))
  return(res)
}

# VSD: symmetric difference between active and selected sets
VSDFunc = function(trueActiveSet,selectedSet) {
  p1 = length(setdiff(trueActiveSet,selectedSet))
  p2 = length(setdiff(selectedSet,trueActiveSet))
  res = p1 + p2
  return(res)
}
