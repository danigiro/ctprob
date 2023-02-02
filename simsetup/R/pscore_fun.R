list2tibble_2Lm <- function(xlist, B){
  df_free <- NULL
  tmp <- lapply(names(xlist), function(comb){
    tmp1 <- lapply(names(xlist[[comb]]), function(k){
      list_tmp <- apply(xlist[[comb]][[k]], 2, function(x){ 
        t(matrix(x, ncol = B))
      }, simplify = FALSE)
      tibble(value = unname(list_tmp), serie = names(list_tmp), k = k, comb = comb)
    })
    do.call(rbind, tmp1)
  })
  do.call(rbind, tmp)
}

list2tibble_1Lm <- function(xlist, B){
  df_free <- NULL
  tmp1 <- lapply(names(xlist), function(k){
    list_tmp <- apply(xlist[[k]], 2, function(x){ 
      t(matrix(x, ncol = B))
    }, simplify = FALSE)
    tibble(value = unname(list_tmp), serie = names(list_tmp), k = k)
  })
  do.call(rbind, tmp1)
}

list2tibble_1Lv <- function(xlist){
  df_free <- NULL
  tmp1 <- lapply(names(xlist), function(k){
    list_tmp <- apply(rbind(xlist[[k]]), 2, function(x){ 
      as.vector(x)
    }, simplify = FALSE)
    tibble(value = unname(list_tmp), serie = names(list_tmp), k = k)
  })
  do.call(rbind, tmp1)
}

crps <- function(y, x){
  list(scoringRules::crps_sample(unlist(y), t(x[[1]])))
}

es <- function(y, x){
  y <- y[[1]]
  x <- x[[1]]
  if(is.vector(y)){
    #list(es_sample(y, t(x[,1,])))
    list(energyScore(y, x[,1,]))
  }else{
    #list(sapply(1:NROW(y), function(i) es_sample(y[i,], t(x[,i,]))))
    list(sapply(1:NROW(y), function(i) energyScore(y[i,], x[,i,])))
  }
}

energyScore <- function(target, sample){
  size = dim(sample)[1]
  randPermutation = sample(size)
  
  samples1 = sample
  samples2 = sample[randPermutation, ]
  
  mTarget = as.matrix(rep(1, size)) %*% as.vector(target)
  term1 = (samples1 - mTarget)
  normst1 = wordspace::rowNorms(term1)
  
  term2 = (samples1 - samples2)
  normst2 = wordspace::rowNorms(term2)
  
  score = mean(normst1) - 0.5 * mean(normst2)
  return(score)
}

arrange_value <- function(x, name){
  x <- setNames(x, name)
  list(simplify2array(x))
}