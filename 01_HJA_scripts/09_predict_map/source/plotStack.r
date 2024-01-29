### plot raster stack


plotStack <- function(x, by =6, nc = 2, nr = 3, ...) {
  
  require(terra)
  
  if(nc * nr != by) stop("by must equal nc x nr")
  
  n <- terra::nlyr(x)
  
  start <- (0:(n%/%by)*by)+1
  if(start[length(start)] > n) start <- start[-length(start)]
  end <- (start + by)-1
  if(end[length(end)] > n) end[length(end)] <- n
  #start
  #end
  
  for(i in seq_along(start)){
    terra::plot(x[[ start[i]:end[i] ]], nc = nc, nr = nr, ...)
  }
  
  
  
}



