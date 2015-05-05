### StarBeast summary results (marginal likelihood summaries)

get_sbeast_results <- function(file) {
    sbeastOrig <- read.csv(file="data/starbeastResults.csv")
    sbeast <- sbeastOrig[, c("groupings", "runs", "chainLength", "dataIncluded",
                             "PS_logLik", "SS_logLik")]
    sbeast <- sbeast[grep("[0-9]$", sbeast$runs), ] # runs with X at the ends are not meant to be included
    sbeast
}

get_sbeast_summary_all <- function(sbeast) {
    sbSummAll <- subset(sbeast, chainLength == "long" & dataIncluded == "all" & groupings != "random")
    meanESU1SS <- mean(subset(sbSummAll, groupings == "allESU1")$SS_logLik)
    meanESU1PS <- mean(subset(sbSummAll, groupings == "allESU1")$PS_logLik)
    sbSummAll$stdSS <- sbSummAll$SS_logLik - meanESU1SS
    sbSummAll$stdPS <- sbSummAll$PS_logLik - meanESU1PS
    sbSummAll$BFSS <- 2 * sbSummAll$stdSS
    sbSummAll$BFPS <- 2 * sbSummAll$stdPS
    sbSummAll$dataIncluded <- "All data"
    sbSummAll$plotGroupings <- sbSummAll$groupings
    levels(sbSummAll$plotGroupings)[levels(sbSummAll$plotGroupings) == "allESU1"] <- "M0"
    levels(sbSummAll$plotGroupings)[levels(sbSummAll$plotGroupings) == "allESU1split"] <- "oversplit"
    levels(sbSummAll$plotGroupings)[levels(sbSummAll$plotGroupings) == "noHawaii"] <- "M5"
    levels(sbSummAll$plotGroupings)[levels(sbSummAll$plotGroupings) == "noRedSea"] <- "M6"
    levels(sbSummAll$plotGroupings)[levels(sbSummAll$plotGroupings) == "noWpac"] <- "M4"
    levels(sbSummAll$plotGroupings)[levels(sbSummAll$plotGroupings) == "noWpacHawaii"] <- "M3"
    levels(sbSummAll$plotGroupings)[levels(sbSummAll$plotGroupings) == "noWpacRedSea"] <- "M2"
    levels(sbSummAll$plotGroupings)[levels(sbSummAll$plotGroupings) == "noWpacRedSeaHawaii"] <- "M1"
    sbSummAll
}


get_sbeast_summary_noCOI <- function(sbeast) {
    sbSummNoCOI <- subset(sbeast, chainLength == "long" & dataIncluded == "noCOI" & groupings != "random")
    meanESU1SSnoCOI <- mean(subset(sbSummNoCOI, groupings == "allESU1")$SS_logLik)
    meanESU1PSnoCOI <- mean(subset(sbSummNoCOI, groupings == "allESU1")$PS_logLik)
    sbSummNoCOI$stdSS <- sbSummNoCOI$SS_logLik - meanESU1SSnoCOI
    sbSummNoCOI$stdPS <- sbSummNoCOI$PS_logLik - meanESU1SSnoCOI
    sbSummNoCOI$dataIncluded <- "No COI"
    sbSummNoCOI
}

get_sbeast_summary_noMt <- function(sbeast) {
    sbSummNoMt <- subset(sbeast, dataIncluded == "noMt")
    meanESU1SSnoMt <- mean(subset(sbSummNoMt, groupings == "allESU1")$SS_logLik)
    sbSummNoMt$stdSS <- sbSummNoMt$SS_logLik - meanESU1SSnoMt
    sbSummNoMt
}
