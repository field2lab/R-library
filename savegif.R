install.packages('animation', repos = 'http://yihui.name/xran')
library(animation)

saveGIF({
  for (i in 1:10) plot(runif(10), ylim = 0:1)
})

ani.options(convert ="C:/Program Files/ImageMagick-6.9.3-Q16/convert.exe")

movie3d(spin3d(axis = c(0, 0, 1)), duration = 3,
        dir = getwd())
writeWebGL(dir = "webGL", filename = file.path(dir, "index.html"))