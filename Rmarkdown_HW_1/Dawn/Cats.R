cats <- read.csv("data/cats.csv")
cats
str(cats)
head(cats)
hist(cats$Bwt)
hist(cats$Bwt, las=1, main = "Distribution of Feline Body Weight", xlab = "Body Weight (kg)", col="blue")
plot(cats$Hwt ~ cats$Bwt, las=1, main = "Relationship between Feline body weight and heart rate", xlab = "Body weight (kg)", ylab = "Heart weight (kg)", pch = 16, cex = 1, col = "dodger blue")
abline(lm(cats$Hwt~cats$Bwt), col="gold", lwd=3)
summary(cats)

cat.scatter<-ggscatter(cats, x = "Hwt", y = "Bwt", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",mainlab = "Relationship between feline body weight and heart rate",
          xlab = "Body Weight (kg)", ylab = "Heart Weight (kg)")
cat.scatter
shapiro.test(cats$Hwt)
shapiro.test(cats$Bwt)
