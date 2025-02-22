library(lmerTest)
df <- sleepstudy
fm <- lmer(Reaction ~ Days + (1 + Days || Subject), df)
test = summary(fm, ddf = "Kenward-Roger")
write.csv(test[10], "Results sleep study lmertest.csv")