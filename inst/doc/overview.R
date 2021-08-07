## ----eval = FALSE--------------------------------------------------------
#  install.packages("biotools")
#  devtools::install_github("arsilva87/biotools")

## ------------------------------------------------------------------------
library(biotools)

## ------------------------------------------------------------------------
data(maize)
head(maize)   # firts 6 rows

## ------------------------------------------------------------------------
M <- lm(cbind(NKPR, ED, CD) ~ family, data = maize)
anova(M, test = "Wilks")

## ------------------------------------------------------------------------
mvpaircomp(M, "family", test = "Wilks", adjust = "bonferroni")

## ------------------------------------------------------------------------
res <- residuals(M)
boxM(data = res, grouping = maize$family)

## ------------------------------------------------------------------------
s <- singh(data = maize[, 1:3], cov = cov(res))
s
plot(s)

## ------------------------------------------------------------------------
d <- D2.dist(data = maize[, 1:3], cov = cov(res))
range(d)

## ------------------------------------------------------------------------
toc <- tocher(d, algorithm = "sequential")
toc

## ----digits=3------------------------------------------------------------
print(toc$distClust, digits = 2)

## ------------------------------------------------------------------------
cop <- cophenetic(toc)  # cophenetic matrix
cor(d, cop)             # cophenetic correlation coefficient

## ------------------------------------------------------------------------
mantelTest(d, cop, nperm = 900)

## ------------------------------------------------------------------------
D2.disc(data = maize[, 1:3], 
        grouping = maize$family, 
        pooled.cov = cov(res))

