# Serial Interval Estimates

# Note from James
# "Note that it's the generation interval that you technically want (and
# what they [Science paper or Imperial College] estimate here) . the
# generation interval is
# the time between exposure in infector and infectee. What you can measure
# usually is the serial interval (time between symptoms in infector and
# infectee) . if you also have the distribution of time to symptoms, then
# you can deconvolve the serial interval into a generation interval. Anyway,
# the main reason to do this is that estimated serial intervals have the
# same mean but an exaggerated variance . this tends to lead to
# underestimates of R0."

# Distributions from Science paper
# https://science.sciencemag.org/content/early/2020/03/30/science.abb6936

# The distribution is lognormal with mean 5.5 days, median 5.2 days and
# standard deviation 2.1 days

mean <- 5.5
median <- 5.2
sd <- 2.1

# Test function
stats <- function(vector) {
 return(list(mean = mean(vector), median = median(vector), sd= sd(vector),
   iqrl = quantile(vector, probs = 0.25), 
   iqru = quantile(vector, probs = 0.75)))
}

# Fitting gamma to mean and variance of the serial interval
estGammaParams <- function(mu, var) {
  shape <- mu^2/var
  scale <- var/mu
  return(params = list(shape = shape, scale = scale))
}

gparams <- estGammaParams(mean, sd^2)

siGamma <- data.frame(x = 1:100,
  fit = dgamma(1:100, gparams$shape, scale = gparams$scale))

testGamma <- rgamma(1e5, gparams$shape, scale = gparams$scale)
stats(testGamma)
plot(siGamma$fit[1:20])

# In the paper they found the best fit for the generation interval
# was a Weibull distribution with:
wshape <- 2.826
wscale <- 5.665

siWeibull <- data.frame(x = 1:100,
  fit = dweibull(1:100, wshape, scale = wscale))

testWeibull <- rweibull(1e5, wshape, scale = wscale)
stats(testWeibull)
plot(siWeibull$fit[1:20])

# Original Imperial serial interval parameters --------------
# Mean = 6.5 and cv = 0.62 and the gammAlt functions for the EnvStats
# Package 

library(EnvStats)

ImperialShape <- 6.5
Imperialcv <- 0.62

siImperialGamma <- data.frame(x = 1:100,
  fit = dgammaAlt(1:100, ImperialShape, cv = Imperialcv))

testImperialGamma <- rgammaAlt(1e5, ImperialShape, cv = Imperialcv)
stats(testImperialGamma)
plot(siImperialGamma$fit[1:20])

# Original fit
# serial.interval <- read.csv("serial_interval.csv")
serial.interval.orig <- read.csv("serial_interval_Original.csv")
plot(serial.interval.orig$fit[1:20])

# Final Serial.Interval selected
# Originally went with the gamma fit but more recent serial interval data
# suggests a lower mean. Serial interval Distribution from WHO update 
# 15 April: 4.8 (IQR: 3.7 – 5.2) from 8 studies. Upper IQR seems a bit tight
# though. In version 1 of the Imperial model this produced
# a reaseonable R0 but version 2 and above seemed to result in an R0 > 3.
# Weibull distribution seems to be a better fit now. 

# Going with the Weibull(replace serial_interval.csv with the new file)
write.csv(siGamma, "serial_interval_Gamma.csv")
write.csv(siWeibull, "serial_interval_Weibull.csv")
