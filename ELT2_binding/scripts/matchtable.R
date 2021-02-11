matchtable = function(..., matchObj=NULL){
# Produces an object that can be used with upsetR
# ...: named list of factors or strings
# matchObj: previous call to matchtable to expand
# 
## Example: 
# matchtable(a=1:2, b=2:3)
#   a b
# 1 1 0
# 2 1 1
# 3 0 1
# 
## Add on to previous output:
# matchtable(c=3:4, matchObj=matchtable(a=1:2, b=2:3))
#   a b c
# 1 1 0 0
# 2 1 1 0
# 3 0 1 1
# 4 0 0 1
  
  things = list(...)
  
  allLevels = rownames(matchObj) # matchObj can be NULL
  
  for (thing in things) {
    allLevels = base::union(allLevels,thing)
  }
  newMatchObj = data.frame(row.names=allLevels)
  
  oldLevels = rownames(matchObj)
  for (thingName in names(matchObj)) {
    columnLevels = oldLevels[ matchObj[[thingName]] == 1]
    newMatchObj[ thingName ] = ifelse( allLevels %in% columnLevels, 1, 0)
  }
  for (thingName in names(things)) {
    newMatchObj[ thingName ] = ifelse( allLevels %in% things[[thingName]], 1, 0)
  }
  
  return(newMatchObj)
}