PriorityQueue <- function(decreasing = T) {
  keys <<- values <<- NULL
  insert <- function(key, value) {
    temp <- c(keys, key)
    ord <- order(temp, decreasing = decreasing)
    keys <<- temp[ord]
    values <<- c(values, c(value))[ord]
  }
  peak <- function() {
    head <- c(keys[[1]], values[[1]])
    return(head)
  }
  pop <- function() {
    head <- c(keys[[1]], values[[1]])
    values <<- values[-1]
    keys <<- keys[-1]
    return(head)
  }
  
  queueLength <- function() {
    return(length(values))
  }
  dump <- function(){
    return(list(keys,values))
  }
  empty <- function() length(keys) == 0
  list(insert = insert, pop = pop, empty = empty, dump = dump, 
       peak = peak, queueLength = queueLength)
}
