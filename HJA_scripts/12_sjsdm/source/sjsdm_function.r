# from 'sjsdm_ukponds_20200320_max.r'

# functions to analyse sjsdm result

#loading function for Polygon drawing
deg2rad <- function(deg) {(deg * pi) / (180)}
rad2deg <- function(rad) {(rad * 180) / (pi)}
re_scale = function(x) return(cov2cor(x))
ff = function(x){(x-min(x))/(max(x)-min(x))}
add_curve = function(p1 = coords[1,], p2 = coords[3,], n = 10, spar = 0.7, col = "black", species = TRUE, lineSeq = 5.0, lwd = 1.0) {
 xxs1 = cos(deg2rad(p1[3]))* seq(0, lineSeq, length.out = n)
 xxs2 = cos(deg2rad(p2[3]))* seq(0, lineSeq, length.out = n)
 yys1 = sin(deg2rad(p1[3]))* seq(0, lineSeq, length.out = n)
 yys2 = sin(deg2rad(p2[3]))* seq(0, lineSeq, length.out = n)
 x = c(rev(xxs1), xxs2[-1])
 y = c(rev(yys1), yys2[-1])
 m = (p1[2] - p2[2])/(p1[1] - p2[1])
 a = rad2deg(atan(m))
 a = -(a+180)
 alpha = deg2rad(a)
 alpha2 = deg2rad(-a)
 rot = matrix(c(cos((alpha)), -sin((alpha)), sin((alpha)), cos((alpha))),2,2)
 rot2 = matrix(c(cos((alpha2)), -sin((alpha2)), sin((alpha2)), cos((alpha2))),2,2)
 tt = cbind(x,y) %*% rot
 sp = smooth.spline(tt[,1], tt[,2],spar = spar,df = 6, w = c(10.0, rep(0.1,nrow(tt)-2), 10.0))
 tt2 = cbind(sp$x, sp$y)
 b = tt2 %*% rot2
 lines(b[,1], b[,2], col = col, lwd = lwd)
 
 x1 = c(cos(deg2rad(p1[3]))*(lineSeq+0.1), cos(deg2rad(p1[3]))*(lineSeq+0.3))
 x2 = c(cos(deg2rad(p2[3]))*(lineSeq+0.1), cos(deg2rad(p2[3]))*(lineSeq+0.3))
 y1 = c(sin(deg2rad(p1[3]))* (lineSeq+0.1), sin(deg2rad(p1[3]))* (lineSeq+0.3))
 y2 = c(sin(deg2rad(p2[3]))* (lineSeq+0.1), sin(deg2rad(p2[3]))* (lineSeq+0.3))
 if(species){
 segments(x0 = x1[1], x1 = x1[2], y0 = y1[1], y1 = y1[2], col = "darkgrey")
 segments(x0 = x2[1], x1 = x2[2], y0 = y2[1], y1 = y2[2],  col = "darkgrey")
 }
}
curve_text = function(pos = 0, label = "", lineSeq = 5.0, reverse = FALSE,middle = FALSE,extend = 1.1,...){
  # inspired by plotrix package
  chars = (strsplit(label,split = "")[[1]])
  char_lens = strwidth(chars)*extend
  char_angles = char_lens / lineSeq
  changrang = range(char_angles)
  char_angles[char_angles < changrang[2]/2] = changrang[2]/2
  
  if(middle & reverse) pos = pos - rad2deg(sum(char_angles)/2)
  if(middle & !reverse) pos = pos + rad2deg(sum(char_angles)/2)
  
  if(reverse) {
    angles = c(deg2rad(pos), deg2rad(pos)+cumsum(char_angles)[-length(chars)])
    angles = angles + char_angles/2
  } else {
    angles = c(deg2rad(pos), deg2rad(pos)-cumsum(char_angles)[-length(chars)])
    angles = angles - char_angles/2
  }
  
  for(i in 1:length(chars)) text(label = chars[i], 
                                 x = cos((angles[i]))*(lineSeq), 
                                 srt = rad2deg(angles[i]) - 90+ 180*reverse,
                                 y = sin((angles[i]))*(lineSeq), 
                                 xpd = NA, adj = c(0.5, 0.5),...)
  return(max(angles))
  
}
