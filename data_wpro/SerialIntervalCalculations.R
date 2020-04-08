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
media <- 5.2
sd <- 2.1

# Test function
stats <- function(vector) {
 return(list(mean = mean(vector), median = median(vector), sd= sd(vector)))
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
# Not sure from report either shape = 6.5 and scale = 0.62 or
# Mean = 6.5 and variance = 0.62??

ImperialParams$shape <- 6.5
ImperialParams$scale <- 0.62

siImperialGamma <- data.frame(x = 1:100,
  fit = dgamma(1:100, ImperialParams$shape, scale = ImperialParams$scale))

testImperialGamma <- rgamma(1e5, ImperialParams$shape, scale = ImperialParams$scale)
stats(testImperialGamma)
plot(siImperialGamma$fit[1:20])

# Original fit
serial.interval <- read.csv("serial_interval.csv")
plot(serial.interval$fit[1:20])

# Final Serial.Interval selected
# Going with the gamma fit (replace serial_interval.csv with the new file)
write.csv(siGamma, "serial_interval_Update.csv")

