'''Illustration showing how the scale parameter changes'''

x <- seq(-4,6,by = .1)
y <- dnorm(x,mean = 0, sd = 1)
b <- dnorm(x, mean = 0, sd = 0.75)
z <- dnorm(x, mean = 0, sd = 0.5)
a <- dnorm(x, mean = 0, sd = 0.25)

plot(x,y,type = 'l', xlab = 'SST (degrees Celsius)', ylab = 'Density', col = 'blue', ylim = c(0,1.6))
lines(x,y,type = 'l', xlab = 'SST (degrees Celsius)', col = 'blue', ylab = 'Density')
lines(x,b,type = 'l', xlab = 'SST (degrees Celsius)', col = 'red', ylab = 'Density')
lines(x,z, type = 'l', xlab = "SST (degrees Celsius", col = 'green', ylab = "Density")
lines(x,a, type = 'l', xlab = "SST (degrees Celsius", col = 'cyan', ylab = "Density")
legend(-4, 0.8, legend = c('1', '0.75', '0.5','0.25'), col = c('blue','red','green','cyan'),
       title = "Scale", lty = 1)