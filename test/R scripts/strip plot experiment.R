library(lmerTest)
df <- read.csv("data/Data battery cell Chapter 11.csv")
fm <- lmer(Y~ 1 + X1 + X2 + X3 + X4 + X5 + X6 +
              X1*X2 + X1*X3 + X1*X4 + X1*X5 + X1*X6 +
              X2*X3 + X2*X4 + X2*X5 + X2*X6 +
              X3*X4 + X3*X5 + X3*X6 +
              X4*X5 + X4*X6 +
              X5*X6 +
              (1 | Whole.Plots) + (1 | Subplots), data=df)
test = summary(fm, ddf = "Kenward-Roger")
write.csv(test[10], "Results battery cell lmertest.csv")