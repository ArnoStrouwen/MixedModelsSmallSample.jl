library(lmerTest)
df <- read.csv("data/Data Pastry Dough Experiment Chapter 7.csv")
names(df)[names(df) == 'Flow.Rate'] <- 'FR'
names(df)[names(df) == 'Moisture.Content'] <- 'MC'
names(df)[names(df) == 'Screw.Speed'] <- 'SS'
names(df)[names(df) == 'Longitudinal.Expansion.Index'] <- 'LEI'
df[c(2)] <-  2*df[c(2)]/(max(df[c(2)])-min(df[c(2)])) - 2*min(df[c(2)])/(max(df[c(2)])-min(df[c(2)])) - 1
df[c(3)] <-  2*df[c(3)]/(max(df[c(3)])-min(df[c(3)])) - 2*min(df[c(3)])/(max(df[c(3)])-min(df[c(3)])) - 1
df[c(4)] <-  2*df[c(4)]/(max(df[c(4)])-min(df[c(4)])) - 2*min(df[c(4)])/(max(df[c(4)])-min(df[c(4)])) - 1

fm <- lmer(LEI ~ 1 + FR + MC + SS + FR*MC + FR*SS + MC*SS + I(FR^2) + I(MC^2) + I(SS^2) + (1 | Day), data=df)
test = summary(fm, ddf = "Kenward-Roger")
write.csv(test[10], "Results pastry dough lmertest.csv")