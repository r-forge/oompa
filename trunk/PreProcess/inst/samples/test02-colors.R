# colors-test.ssc

theta <- (0:360)*pi/180
x <- cos(theta)
y <- sin(theta)
xp <- x > 0
yp <- y > 0
colors <- list(ColorCoding(xp&yp, COLOR.BORING),
               ColorCoding(xp&!yp, COLOR.TOP.TEN),
               ColorCoding(!xp&yp, COLOR.BOTTOM.TEN),
               ColorCoding(!xp&!yp, COLOR.CONFIDENCE.CURVE))
plot(ColorCodedPair(x, y, colors))

plot(ColorCodedPair(theta, x, colors))

plot(ColorCodedPair(theta, y, colors),
     xlab='angle in radians', ylab='sine', main='colored sine')

fac <- factor(rep(c('left', 'right'), c(180, 181)))
colors <- colorCode(fac, c('blue', 'red'))
plot(ColorCodedPair(x, y, colors))

# cleanup
rm(x, y, xp, yp, theta, colors, fac)

colorList <- c(COLOR.BORING, COLOR.SIGNIFICANT,
               COLOR.EXPECTED, COLOR.OBSERVED,
               COLOR.PERMTEST, COLOR.FITTED,
               COLOR.CENTRAL.LINE, COLOR.CONFIDENCE.CURVE,
               COLOR.BAD.REPLICATE, COLOR.WORST.REPLICATE,
               COLOR.FOLD.DIFFERENCE, COLOR.BAD.REPLICATE.RATIO,
               COLOR.TOP.TEN, COLOR.BOTTOM.TEN,
               COLOR.TOP.TEN.SOLO, COLOR.BOTTOM.TEN.SOLO
               )
plot(c(1,4), c(1,4), type='n')
for (i in 1:4) {
  for (j in 1:4) {
    points(i,j, col=colorList[i + 4*(j-1)], pch=16, cex=4)
  }
}
rm(colorList, i, j)
