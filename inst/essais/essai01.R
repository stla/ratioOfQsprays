library(ratioOfQsprays)

x <- qlone(1)

roq <- new("ratioOfQsprays", numerator = x^7 - 1, denominator = x - 1)
ratioOfQsprays:::simplifyRatioOfQsprays(roq)

(1 - x^7) / (1 - x)

as.ratioOfQsprays(1-x^7) / as.ratioOfQsprays(1-x)
