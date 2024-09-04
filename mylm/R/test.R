#install.packages("car")
library(car)

data(SLID, package = "carData")
SLID <- SLID[complete.cases(SLID), ]

lm1 = lm(wages ~ education + age, SLID)
mylm1 = mylm(wages ~ education + age, SLID)
mylm(wages ~ education + age, SLID)
summary(lm1)
summary(mylm1)
min(mylm1$residuals)
