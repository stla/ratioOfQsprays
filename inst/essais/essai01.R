library(ratioOfQsprays)

x <- qlone(1)

roq <- new("ratioOfQsprays", numerator = 1 - x^7, denominator = 1 - x)
ratioOfQsprays:::simplifyRatioOfQsprays(roq)

(1 - x^7) / (1 - x)

as.ratioOfQsprays(1-x^7) / as.ratioOfQsprays(1-x)
