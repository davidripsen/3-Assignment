#install.packages("ctsmr", repo = "http://ctsm.info/repo/dev")
#install.packages("pkgbuild")

# For git pushing
## git push https://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@github.com/davidripsen/3-Assignment.git

library(ctsmr)
library(splines)
source("CompEx3_E18/sdeTiTm.R")
# Load data
if (Sys.info()[7] == "davidipsen")
  {path <- "~/Documents/DTU/3. Semester (MSc)/Advanced Time Series/Assignments/3-Assignment/CompEx3_E18/"
} else {path <- "CompEx3_E18/"}
load(paste0(path,"Exercise3.RData"))
AllDat

fit1 <- sdeTiTm(AllDat,AllDat$yTi1,AllDat$Ph1)

summary(fit1,extended=TRUE)

Hour <- as.numeric(strftime(AllDat$date, format="%H"))

Pred <- predict(fit1)
plot(Pred[[1]]$state$pred$Ti - AllDat$yTi1 ~ Hour)
# What is going on 10 AM?
# Try to fir a varying effective window area


plot(AllDat$Gv ~ Hour)


idx <- (Hour>8 & Hour < 23) # It is impossible to fit a window area for the hours without any sun, so we limit the window area estimation to the hours with sun.
bs = bs(Hour[idx],df=5,intercept=TRUE) # Dvs. 4 knots / 5 basis splines

# What does the splines look like?
plot(bs[14:27,1],type='l')
lines(bs[ 14:27,2])
lines(bs[ 14:27,3])
lines(bs[ 14:27,4])
lines(bs[ 14:27,5])

bs1 <- bs2 <- bs3 <- bs4 <- bs5 <- bs6 <- numeric(dim(AllDat)[1])

bs1[idx] = bs[,1]
bs2[idx] = bs[,2]
bs3[idx] = bs[,3]
bs4[idx] = bs[,4]
bs5[idx] = bs[,5]

AllDat$bs1 = bs1
AllDat$bs2 = bs2
AllDat$bs3 = bs3
AllDat$bs4 = bs4
AllDat$bs5 = bs5


### You will have to implement sdeTITmAv ###
source(paste0(path,"sdeTiTmAv.R"))
fit2 <- sdeTiTmAv(AllDat,AllDat$yTi1,AllDat$Ph1)

plot(bs[14:27,1]*fit2$xm[3]+bs[14:27,2]*fit2$xm[4]+bs[14:27,3]*fit2$xm[5]+bs[14:27,4]*fit2$xm[6]+bs[14:27,5]*fit2$xm[7],type='l')
fit2
