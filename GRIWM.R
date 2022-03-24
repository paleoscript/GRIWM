
## supplementary material for MAGYARI et al. 2022 (Sci. Rep.) ##
## https://doi.org/10.21203/rs.3.rs-778658/v1                 ##

## written by Attila VIRAG (1, 2) and Bence SZABO (1, 3)      ##
## 1: MTA-MTM-ELTE Research Group for Paleontology, Budapest  ##
## 2: Dept. of Mineralogy and Geology, University of Debrecen ##
## 3: Institute of Translational Medicine, Semmelweis Univ.   ##

#### REFERENCES ####
## ref(1): McINERNY et al. 2006 (Conserv. Biol.)              ##
## https://doi.org/10.1111/j.1523-1739.2006.00377.x           ##
## ref(2): BRADSHAW et al. 2012 (Quat. Sci. Rev.)             ##
## https://doi.org/10.1016/j.quascirev.2011.11.021            ##


#### FUNCTIONS ####


# calculates the probability of a sighting..
# after a given period since the last known occurrence
# McInerny et al. method, see ref(1)
prob <- function(td, tn, n) {
  
  return((1 - (n/tn))^td)
}
# n  = the number of sightings
# tn = the initial period of observation,..
# tn = oldest sighting - youngest sighting
# td = time since last observation
# returns a value between 0 and 1..
# if all input parameters are positive


# calculates the time since the last sighting..
# after which the probability of a new sighting..
# falls below a certain threshold value
# p < 0.05 imply extinction
# McInerny et al. method, see ref(1)
term <- function(tn, n, alpha = 0.05, step = 1) {
  
  td <- 0
  p  <- 1
  for (exp in 3:0) {
    while (p > alpha) {
      td <- td + step*10^exp
      p  <- prob(td, tn, n)
    }
    td <- td - step*10^exp
    p  <- prob(td, tn, n)
  }
  
  return(td)
}
# n     = the number of sightings
# tn    = the initial period of observation
# alpha = a pre-defined threshold value
# step  = the time value that is added to td..
# after each iterative step (time resolution)
# resolution refining steps are introduced..
# for performance optimization


# IWM (Inverse Weighted McInerny et al. method), see ref(3)
IWM <- function(sightings, alpha = 0.05, step = 1) {
  
  # sorting sightings from the most recent towards the oldest
  sightings <- sightings[order(sightings)]
  # number of sightings
  n  <- length(sightings)
  # most recent sighting
  t1 <- min(sightings)
  # initial period of observations
  tn <- max(sightings)-t1
  
  thetas <- numeric()
  omegas <- numeric()
  for(k in 2:n) {
    ti     <- sightings[k]
    tn     <- ti - t1
    thetas <- c(thetas, term(tn, k, alpha, step))
    omegas <- c(omegas, 1/tn)
  }
  
  term       <- sum(omegas*thetas)/sum(omegas)
  ans        <- list(term, thetas, omegas)
  names(ans) <- c("term", "thetas", "omegas")
  
  return(ans)
}
# sightings = a numeric array of radiocarbon ages
# alpha and step values are passed on to the term() function
# result$term = a terminal age calculated using the IWM method
# result$thetas = a list of ages using sequentially longer..
# time series in each iteration (from k = 2 to n;..
# where the youngest sighting is always included)
# result$omegas = a list of weights used for the ..
# weighted summation of the above mentioned ages


# GRIWM (Gaussian-Resampled IWM method), see ref(3)
GRIWM <- function(sightings, sdevs, iterations = 10000, 
                  conf = 0.95, alpha = 0.05, step = 1,
                  silent = TRUE) {
  
  # number of sightings
  n  <- length(sightings)
  
  if (!silent) {
    progress <- 0
    previous <- 0
  }
  
  terms <- numeric()
  m <- matrix(NA, nrow = iterations, ncol = n)
  for (i in 1:n) {
    # resample each radiometric date in the series..
    # from a Gaussian distribution
    m[, i] <- rnorm(iterations, sightings[i], sdevs[i])
  }
  
  terms <- numeric()
  for (i in 1:iterations) {
    terms[i] <- IWM(m[i, ], alpha, step)$term
    if (!silent) {
      progress <- progress + 100/iterations
      if (previous != round(progress)) {
        print(paste(round(progress), "%"))
      }
      previous <- round(progress)
    }
  }
  # terms <- numeric()
  # m <- matrix(NA, nrow = iterations, ncol = n)
  # for (i in 1:iterations) {
  #   for (j in 1:n) {
  #     # resample each radiometric date in the series..
  #     # from a Gaussian distribution
  #     m[i, j] <- rnorm(1, sightings[j], sdevs[j])
  #   }
  #   terms[i] <- IWM(m[i, ], alpha, step)$term
  #   if (!silent) {
  #     progress <- progress + 100 / iterations
  #     if (previous != round(progress, 1)) {
  #       print(paste(round(progress, 1), "%"))
  #     }
  #     previous <- round(progress, 1)
  #   }
  # }
  
  a          <- (1-conf)/2
  lower      <- quantile(terms, a)
  upper      <- quantile(terms, 1-a)
  ans        <- list(terms, lower, upper)
  names(ans) <- c("terms", "lower", "upper")
  
  return(ans)
}
# sightings = a numeric array of radiocarbon ages
# sdevs = a numeric array of standard deviations (sigmas).. 
# for the radiocarbon ages in the same order
# conf = a pre-defined confidence level for the resulting limits
# alpha and step values are passed on to the IWM() function
# if silent is TRUE, progress report is suppressed
# a high iteration number and a small step value..
# can cause performance issues
# result$terms = terminal ages for each iteration..
# calculated using the IWM method
# result$lower = lower limit for the confidence interval
# result$upper = upper limit for the confidence interval


#### EXAMPLE ####


#### STEP 1: LOAD INPUT SPREADSHEATS ####


# Set working directory (in RStudio)
setwd(paste(dirname(rstudioapi::getActiveDocumentContext()$path)))

# Loading radiocarbon data
mammoth  <- read.csv("mammoth.csv", sep = ";", header = FALSE)
reindeer <- read.csv("reindeer.csv", sep = ";", header = FALSE)

colnames(mammoth)  <- c("loc", "upp", "low", "age")
colnames(reindeer) <- c("loc", "upp", "low", "age")
# loc = locality
# upp = mean age + 2-sigma   (upper limit)
# low = mean age - 2-sigma   (lower limit)
# age = mean age

# Sorting data into ascending order by geological age
mammoth  <- mammoth[order(mammoth$age), ]
reindeer <- reindeer[order(reindeer$age), ]


#### STEP 2: VARIABLES ####


# number of sightings
mammoth_n   <- nrow(mammoth)
reindeer_n  <- nrow(reindeer)
# most recent sighting
mammoth_t1  <- min(mammoth$age)
reindeer_t1 <- min(reindeer$age)
# initial period of observations
mammoth_tn  <- max(mammoth$age)-mammoth_t1
reindeer_tn <- max(reindeer$age)-reindeer_t1
# standard deviations (1-sigma)
mammoth_sdevs  <- (mammoth$upp-mammoth$low)/4
reindeer_sdevs <- (reindeer$upp-reindeer$low)/4


#### STEP 3: McINERNY ET AL. METHOD ####


mammoth_term  <- mammoth_t1-term(mammoth_tn, mammoth_n)
reindeer_term <- reindeer_t1-term(reindeer_tn, reindeer_n)


#### STEP 4: IWM METHOD ####


mammoth_res  <- IWM(mammoth$age)
reindeer_res <- IWM(reindeer$age)

mammoth_t1 - mammoth_res$term
reindeer_t1 - reindeer_res$term

#### STEP 5: GRIWM METHOD ####


mammoth_conf  <- GRIWM(mammoth$age, mammoth_sdevs, silent = FALSE)
reindeer_conf <- GRIWM(reindeer$age, reindeer_sdevs, silent = FALSE)

mammoth_t1 - mean(mammoth_conf$terms)
mammoth_t1 - median(mammoth_conf$terms)
mammoth_t1 - mammoth_conf$lower
mammoth_t1 - mammoth_conf$upper

reindeer_t1 - mean(reindeer_conf$terms)
reindeer_t1 - median(reindeer_conf$terms)
reindeer_t1 - reindeer_conf$lower
reindeer_t1 - reindeer_conf$upper


#### FIG 1 ####


{
  xmin <- 5000
  xmax <- 50000
  plot(x = 0, y = 0, type = "n", 
       axes = FALSE, xlab = "", ylab = "", 
       xlim = c(xmax, xmin), ylim = c(-1, 5.5))
  axis(side = 1, at = seq(xmax, xmin, -1000), 
       labels = FALSE, tck = -0.02)
  axis(side = 1, at = seq(xmax, xmin, -5000), 
       labels = seq(xmax/1000, xmin/1000, -5))
  axis(side = 2, at = 0:1, las = 2)
  axis(side = 2, at = seq(0, 1, 0.1), 
       labels = FALSE, tck = -0.02)
  axis(side = 2, at = 3:4, las = 2, labels = 0:1)
  axis(side = 2, at = seq(3, 4, 0.1), 
       labels = FALSE, tck = -0.02)
  points(mammoth$age, rep(1.5, mammoth_n), pch = "×")
  points(reindeer$age, rep(4.5, reindeer_n), pch = "×")
  segments(x0 = max(mammoth$age), x1 = mammoth_t1, 
                    y0 = 1, lwd = 2)
  segments(x0 = max(reindeer$age), x1 = reindeer_t1, 
           y0 = 4, lwd = 2)
  title(main = "McInerny et al. (2006) method",
        xlab = "age (ka cal BP)", ylab ="probability")
  text(x = xmin+(xmax-xmin)/2, y = 2.5, "mammoths")
  text(x = xmin+(xmax-xmin)/2, y = 5.5, "reindeer")
  
  mammoth_probs <- numeric()
  for(i in 1:(mammoth_t1-xmin)) {
    mammoth_probs[i] <- prob(i, mammoth_tn, mammoth_n)
  }
  reindeer_probs <- numeric()
  for(i in 1:(reindeer_t1-xmin)) {
    reindeer_probs[i] <- prob(i, reindeer_tn, reindeer_n)
  }
  
  points(x = mammoth_t1 - 1:(mammoth_t1-xmin), 
         y = mammoth_probs, type = "l", lwd = 2)
  points(x = reindeer_t1 - 1:(reindeer_t1-xmin), 
         y = reindeer_probs + 3, type = "l", lwd = 2)
  segments(x0 = mammoth_term, y0 = 0, y1 = 1, lty = 2)
  segments(x0 = mammoth_term, y0 = -1, y1 = -0.6, lty = 2)
  segments(x0 = reindeer_term, y0 = 3, y1 = 4, lty = 2)
  segments(x0 = reindeer_term, y0 = -1, y1 = -0.6, lty = 2)
}


#### FIG 2 ####


{
  ymin <- 8000
  ymax <- 18000
  plot(x = 0, y = 0, type = "n", 
       axes = FALSE, xlab = "", ylab = "",
       ylim = c(ymax, ymin), xlim = c(0, mammoth_n))
  axis(side = 2, at = seq(ymax, ymin, -200), 
       labels = FALSE, tck = -0.02)
  axis(side = 2, at = seq(ymax, ymin, -1000), 
       labels = seq(ymax/1000, ymin/1000, -1), las = 2)
  axis(side = 1, at = seq(0, mammoth_n, 1), 
       labels = FALSE, tck = -0.02)
  axis(side = 1, at = seq(0, mammoth_n, 5), 
       labels = seq(0, mammoth_n, 5))
  title( xlab = "k (number of used observations)", 
         ylab = "age (ka cal BP)",
         main = "GRIWM method (Bradshaw et al. 2012)")
  
  polygon(x = c(-1, -1, mammoth_n + 1, mammoth_n + 1), 
          y = c(reindeer_t1 - reindeer_conf$lower,
                reindeer_t1 - reindeer_conf$upper,
                reindeer_t1 - reindeer_conf$upper,
                reindeer_t1 - reindeer_conf$lower),
          col = "grey90", border = NA)
  
  polygon(x = c(-1, -1, mammoth_n + 1, mammoth_n + 1), 
          y = c(mammoth_t1 - mammoth_conf$lower,
                mammoth_t1 - mammoth_conf$upper,
                mammoth_t1 - mammoth_conf$upper,
                mammoth_t1 - mammoth_conf$lower),
          col = "grey60", border = NA)
  
  points(2:(reindeer_n), reindeer_t1 - reindeer_res$thetas, 
         type = "o", pch = 16, lty = 2, cex = 0.7, col = "grey50")
  abline(h = reindeer_t1, lwd = 2, lty = 2, col = "grey50")
  abline(h = reindeer_t1 - reindeer_res$term, 
         lty = 2, col = "grey50")
  
  points(2:(mammoth_n), mammoth_t1 - mammoth_res$thetas, 
         type = "o", pch = 16, cex = 0.7)
  abline(h = mammoth_t1, lwd = 2)
  abline(h = mammoth_t1 - mammoth_res$term)

  legend("topright", cex = 0.7, inset = 0.01,
         fill   = c(rep(NA, 3), "grey60"),
         border = rep(NA, 4),
         pch = c(16, rep(NA, 3)),
         lty = c(rep(1, 3), NA),
         lwd = c(1, 2, 1, NA),
         col = c(rep("black", 3), NA),
         text.font = c(2, rep(1, 3)),
         legend = c("mammoth (IWM)", 
                    "last observation", 
                    "terminal age",
                    "95% conf. interval"))
  
  legend("topleft", cex = 0.7, inset = 0.01,
         fill   = c(rep(NA, 3), "grey90"),
         border = rep(NA, 4),
         pch    = c(16, rep(NA, 3)),
         lty    = c(rep(2, 3), NA),
         lwd    = c(1, 2, 1, NA),
         col    = c(rep("grey50", 3), NA),
         text.font = c(2, rep(1, 3)),
         legend    = c("reindeer (IWM)", 
                       "last observation", 
                       "terminal age",
                       "95% conf. interval"))
}


#### ADDITIONAL CALCULATIONS ####


mammoth1 <- mammoth[ 1:17, ]
mammoth2 <- mammoth[18:29, ]

mammoth1_n      <- nrow(mammoth1)
mammoth2_n      <- nrow(mammoth2)
mammoth1_t1     <- min(mammoth1$age)
mammoth2_t1     <- min(mammoth2$age)
mammoth1_tn     <- max(mammoth1$age)-mammoth1_t1
mammoth2_tn     <- max(mammoth2$age)-mammoth2_t1
mammoth1_sdevs  <- (mammoth1$upp-mammoth1$low)/4
mammoth2_sdevs  <- (mammoth2$upp-mammoth2$low)/4
mammoth1_term   <- mammoth1_t1-term(mammoth1_tn, mammoth1_n)
mammoth2_term   <- mammoth2_t1-term(mammoth2_tn, mammoth2_n)
mammoth1_res    <- IWM(mammoth1$age)
mammoth2_res    <- IWM(mammoth2$age)
mammoth1_conf   <- GRIWM(mammoth1$age, mammoth1_sdevs, silent = FALSE)
mammoth2_conf   <- GRIWM(mammoth2$age, mammoth2_sdevs, silent = FALSE)


#### FIG 3 ####


{
  xmin <- 5000
  xmax <- 50000
  plot(0, 0, type = "n", axes = FALSE,xlab = "", ylab = "", 
       xlim = c(xmax, xmin), ylim = c(-1, 5.5))
  axis(side = 1, at = seq(xmax, xmin, -1000), 
       labels = FALSE, tck = -0.02)
  axis(side = 1, at = seq(xmax, xmin, -5000), 
       labels = seq(xmax/1000, xmin/1000, -5))
  axis(side = 2, at = 0:1, las = 2)
  axis(side = 2, at = seq(0, 1, 0.1), 
       labels = FALSE, tck = -0.02)
  axis(side = 2, at = 3:4, las = 2, labels = 0:1)
  axis(side = 2, at = seq(3, 4, 0.1), 
       labels = FALSE, tck = -0.02)
  points(mammoth1$age, rep(1.5, mammoth1_n), pch = "×")
  points(mammoth2$age, rep(1.5, mammoth2_n), pch = "×")
  points(reindeer$age, rep(4.5, reindeer_n), pch = "×")
  segments(x0 = max(mammoth1$age), x1 = mammoth1_t1, 
           y0 = 1, lwd = 2)
  segments(x0 = max(mammoth2$age), x1 = mammoth2_t1, 
           y0 = 1, lwd = 2)
  segments(x0 = max(reindeer$age), x1 = reindeer_t1, 
           y0 = 4, lwd = 2)
  title(main = "McInerny et al. (2006) method",
        xlab = "age (ka cal BP)", ylab ="probability")
  text(x = xmin+(xmax-xmin)/2, y = 2.5, "mammoths")
  text(x = xmin+(xmax-xmin)/2, y = 5.5, "reindeer")
  
  mammoth1_probs <- numeric()
  for(i in 1:(mammoth1_t1-xmin)) {
    mammoth1_probs[i] <- prob(i, mammoth1_tn, mammoth1_n)
  }
  mammoth2_probs <- numeric()
  for(i in 1:(mammoth2_t1-20000)) {
    mammoth2_probs[i] <- prob(i, mammoth2_tn, mammoth2_n)
  }
  reindeer_probs <- numeric()
  for(i in 1:(reindeer_t1-xmin)) {
    reindeer_probs[i] <- prob(i, reindeer_tn, reindeer_n)
  }
  
  points(x = mammoth1_t1 - 1:(mammoth1_t1-xmin), 
         y = mammoth1_probs, type = "l", lwd = 2)
  points(x = mammoth2_t1 - 1:(mammoth2_t1-20000), 
         y = mammoth2_probs, type = "l", lwd = 2)
  points(x = reindeer_t1 - 1:(reindeer_t1-xmin), 
         y = reindeer_probs + 3, type = "l", lwd = 2)
  segments(x0 = mammoth1_term, y0 = 0, y1 = 1, lty = 2)
  segments(x0 = mammoth2_term, y0 = 0, y1 = 1, lty = 2)
  segments(x0 = mammoth1_term, y0 = -1, y1 = -0.6, lty = 2)
  segments(x0 = mammoth2_term, y0 = -1, y1 = -0.6, lty = 2)
  segments(x0 = reindeer_term, y0 = 3, y1 = 4, lty = 2)
  segments(x0 = reindeer_term, y0 = -1, y1 = -0.6, lty = 2)
}


#### FIG 4 ####


{
  ymin <- 12000
  ymax <- 35000
  plot(x = 0, y = 0, type = "n", 
       axes = FALSE, xlab = "", ylab = "",
       ylim = c(ymax, ymin), xlim = c(0, mammoth_n + 1))
  axis(side = 2, at = seq(ymax, ymin, -200), 
       labels = FALSE, tck = -0.02)
  axis(side = 2, at = seq(ymax, ymin, -1000), 
       labels = seq(ymax/1000, ymin/1000, -1), las = 2)
  axis(side = 1, at = seq(0, mammoth1_n, 1), 
       labels = FALSE, tck = -0.02)
  axis(side = 1, at = seq(0, mammoth1_n, 5), 
       labels = seq(0, mammoth1_n, 5))
  axis(side = 1, at = seq(mammoth1_n + 1, mammoth_n + 1, 1), 
       labels = FALSE, tck = -0.02)
  axis(side = 1, at = seq(mammoth1_n + 1, mammoth_n + 1, 5), 
       labels = seq(0, mammoth2_n, 5))
  title( xlab = "k (number of used observations)", 
         ylab = "age (ka cal BP)",
         main = "GRIWM method (Bradshaw et al. 2012)")

  polygon(x = c(1, 1, mammoth_n + 1, mammoth_n + 1), 
          y = c(mammoth2_t1 - mammoth2_conf$lower,
                mammoth2_t1 - mammoth2_conf$upper,
                mammoth2_t1 - mammoth2_conf$upper,
                mammoth2_t1 - mammoth2_conf$lower),
          col = "grey80", border = NA)
  
  polygon(x = c(1, 1, mammoth1_n, mammoth1_n), 
          y = c(mammoth1_t1 - mammoth1_conf$lower,
                mammoth1_t1 - mammoth1_conf$upper,
                mammoth1_t1 - mammoth1_conf$upper,
                mammoth1_t1 - mammoth1_conf$lower),
          col = "grey80", border = NA)
  
  points(2:(mammoth1_n), mammoth1_t1 - mammoth1_res$thetas, 
         type = "o", pch = 16, cex = 0.7)
  polygon(x = c(1, 1, mammoth1_n, mammoth1_n), 
          y = c(mammoth1_t1, max(mammoth1$age), 
                max(mammoth1$age), mammoth1_t1),
          col = "grey20", border = NA)
  segments(x0 = 1, x1 = mammoth1_n, 
           y0 = mammoth1_t1 - mammoth1_res$term)
  
  points((mammoth1_n + 3):(mammoth_n + 1), mammoth2_t1 - mammoth2_res$thetas, 
         type = "o", pch = 16, cex = 0.7)
  polygon(x = c(1, 1, mammoth_n + 1, mammoth_n + 1), 
          y = c(mammoth2_t1, ymax, ymax, mammoth2_t1),
          col = "grey20", border = NA)
  segments(x0 = 1, x1 = mammoth_n + 1, 
           y0 = mammoth2_t1 - mammoth2_res$term)
  
  points(rep(0, mammoth_n), mammoth$age, pch = "×", col = "black")
  
  legend("topright", cex = 0.7, inset = 0.01,
         fill   = c(rep(NA, 3), "grey60", "grey20"),
         border = rep(NA, 5),
         pch = c(4, 16, rep(NA, 2), NA),
         lty = c(NA, rep(1, 2), rep(NA, 2)),
         lwd = c(NA, rep(1, 2), rep(NA, 2)),
         col = c(rep("black", 3), NA, NA),
         legend = c("radiocarbon data",
                    "mammoth (IWM)",
                    "terminal age",
                    "95% conf. interval",
                    "known occupation"))
}

