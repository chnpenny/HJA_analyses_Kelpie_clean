#### Function version for Baseiro irreplaceability ####

irrAB <- function(x, pc = NULL, type= c("total", "median"), tr = NULL, r = NULL, res = c("both", "alpha", "beta")){
  
  # x is a matrix of sites by species, with values of species pop/area, etc within sites
  # pc is % threshold of total area
  # tr is a vector of thresholds or single threshold.
  # r is an optional r template to write results to and return only beta raster
  # res value to return
  
  res <- match.arg(res)
  type <- match.arg(type)
  
  ## Total available area by species
  Rs <- colSums(x, na.rm = T) # range(Rs)
  
  # each species with % of range
  if(!is.null(pc)) {
    
    tr <- switch(type, 
                 total = pc * Rs,
                 median = pc * median(Rs))
  }
  Rs_t <- Rs - tr
  
  ## Alpha irreplaceability
  ## If site does not have species, already is 0
  ## If target is 0 (only if specific species targets)
  if(length(tr) == ncol(x)) x[, tr ==0] <- 0
  
  ## assign 1 where all cells are needed to complete target (where target is greater than total available )
  # x[, tr >= Rs] # this filter columns for which target is larger than total area
  # [x[, tr >= Rs] > 0] # this filters whole matrix for cell > 0
  x[, tr >= Rs][x[, tr >= Rs] > 0] <- 1
  
  # for remainder... 
  tmp <- t(t(x[, tr < Rs]) / Rs_t[tr < Rs]) # divide each row by Rs_t for those species with tr < Rs
  tmp[tmp > 1] <- 1 # adjust back to 1 if some sites have more than in total - target # range(tmp, na.rm =T)
  x[, tr < Rs] <- tmp
  
  ## Beta irreplaceability
  # rowwise
  beta <- apply(x, 1, function(z) 1 - prod(1- z)) # range(beta, na.rm = T)
  
  ## put back into raster
  if(!is.null(r)) {
    beta.r <- raster(r)
    beta.r[] <- beta
    return(beta.r)
  } else switch(res, 
                alpha = return(x),
                beta = return(beta),
                both =  return(list(alpha=x, beta=beta))
  )
  
}
