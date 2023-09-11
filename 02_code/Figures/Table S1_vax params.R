# Table S1 - parameter values for vaccines
library(malariasimulation)

rtss_booster_boosted <- rtss_booster_profile
rtss_booster_boosted$cs <- c(6.37008, 0.35)

r21_profile <- rtss_profile
r21_profile$vmax <- 0.84#0.8441944
r21_profile$alpha <- 0.99#0.9879701
r21_profile$beta <- 455.69#455.6919
r21_profile$cs <- c(9.32, 0.839)
r21_profile$rho <- c(0.791, 0.607)
r21_profile$ds <- c(3.79, 0.165)
r21_profile$dl <- c(6.28, 0.433)

r21_booster_profile <- r21_profile
r21_booster_profile$rho <- c(0.0821, 0.547)
r21_booster_profile$cs <- c(9.24, 0.719)

r21_booster_profile2 <- r21_booster_profile
r21_booster_profile2$cs <- c(9.02, 0.845)


tbl <- data.frame(rtss = unlist(rtss_profile),
                  rtssbooster = unlist(rtss_booster_profile),
                  rtssboosterboosted = unlist(rtss_booster_boosted),
                  r21 = unlist(r21_profile),
                  r21booster1 = unlist(r21_booster_profile),
                  r21booster2 = unlist(r21_booster_profile2)) 
tbl$parameter <- rownames(tbl)
print(tbl, noSpaces = T) |> write.table("clipboard", sep = "\t", row.names = FALSE)
