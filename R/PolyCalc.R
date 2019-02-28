# Output result of polymonial function. 
# For general use with any order of polynomial. 
# Written by Ivan Baldry.
# Translated to R by Leon Drygala
PolyCalc = function( coefficients, xValues){
  yValues <- rep(0, length(xValues)) + coefficients[1]
  for( i in 2:length(coefficients))
    yValues <- yValues + coefficients[i]*xValues^(i - 1)
  return = yValues
}