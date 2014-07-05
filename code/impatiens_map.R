### ---- init-map ----
library(maps)
library(plyr)
library(ggplot2)


impMap <- impDB[nzchar(impDB$consensusESU), ]
center <- 200
impMap$Long.recenter <- ifelse(impMap$Long < center - 180, impMap$Long + 360, impMap$Long)

####
### By Scott Chamberlain http://r.789695.n4.nabble.com/Geographic-distance-between-lat-long-points-in-R-td3442338.html
## Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

## Calculates the geodesic distance between two points specified by
## radian latitude/longitude using the Haversine formula
gcd.hf <- function(long1, lat1, long2, lat2) {
    R <- 6371 # Earth mean radius [km]
    delta.long <- (long2 - long1)
    delta.lat <- (lat2 - lat1)
    a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
    c <- 2 * asin(min(1,sqrt(a)))
    d = R * c
    return(d) # Distance in km
}

## Fxn to calculate matrix of distances between each two sites
CalcDists <- function(latlongs) {
    name <- list(rownames(latlongs), rownames(latlongs))
    n <- nrow(latlongs)
    z <- matrix(0, n, n, dimnames = name)
    for (i in 1:n) {
        for (j in 1:n) z[i, j] <- gcd.hf(long1 = latlongs[i, 1],
                                         lat1 = latlongs[i, 2], long2 = latlongs[j, 1], lat2 = latlongs[j,2])
    }
    z <- as.dist(z)
    return(z)
}


thinCoords <- function(coords=impMap) {
    eachESU <- unique(coords$consensusESU)
    res <- vector("list", length(eachESU))  
    tmpDt <- coords 
    tmpDt <- tmpDt[!is.na(tmpDt$Lat) & !is.na(tmpDt$Long), ]
    tmpDt$Lat2 <- tmpDt$Lat
    tmpDt$Long2 <- tmpDt$Long
    tmpDt$Long2.recenter <- tmpDt$Long.recenter    
    d <- CalcDists(cbind(deg2rad(tmpDt$Long), deg2rad(tmpDt$Lat)))
    d <- as.matrix(d)
    d[upper.tri(d, diag = TRUE)] <- NA
    lbl <- cbind(dimnames(d)[[1]][row(d)], dimnames(d)[[1]][col(d)])
    d <- as.vector(d)
    whichDup <- lbl[which(d < 200), ]      
    if (length(whichDup) > 0) {
        if(length(whichDup) == 2) {
            whichDup <- matrix(whichDup, ncol = 2, byrow = TRUE) # to deal with single result
        }
        for (j in 1:nrow(whichDup)) {
            tmpDt$Lat2[as.numeric(whichDup[j,2])] <- tmpDt$Lat2[as.numeric(whichDup[j,1])]
            tmpDt$Long2[as.numeric(whichDup[j,2])] <- tmpDt$Long2[as.numeric(whichDup[j,1])]
            tmpDt$Long2.recenter[as.numeric(whichDup[j,2])] <- tmpDt$Long2.recenter[as.numeric(whichDup[j,1])]
        }
    }                                        
    tmpDt
}

impMapThin <- thinCoords(impMap)
impMapThin <- impMapThin[!duplicated(paste(impMapThin$consensusESU, impMapThin$Lat2, impMapThin$Long2)), ]

iwp <- map_data("world2")

### ---- impatiens-map ----
### Global impatiens map
## gMap <- ggplot(impMapThin) + annotation_map(iwp, fill = "gray50", colour = "gray50") +
##     geom_point(aes(x = Long2.recenter, y = Lat2, colour = consensusESU), data = impMapThin,
##                position = position_dodge(width = 5, height = 1)) +
##     xlim(c(15,310)) + ylim(c(-40, 40)) +
##     coord_map(projection = "mercator", orientation = c(90, 160, 0)) +
##     theme(panel.background = element_rect(fill = 'aliceblue'))

## gMap

### ---- impatiens-map-WA ----
## 1. WA + EP + Gala
c1Pal <- impPal[c("WA", "EP", "Gala")]
tmpMap <- subset(impMap, consensusESU %in% c("WA", "EP", "Gala"))
ggplot(impMap) + annotation_map(iwp, fill = "gray40", colour = "gray40") +
    geom_point(aes(x = Long.recenter, y = Lat, colour = consensusESU), data = tmpMap,
               position = position_dodge(width = 1.5), shape=16, size=3) +
    scale_colour_manual(values=c1Pal) + 
    xlim(c(250,310)) + ylim(c(-25, 25)) + ylab("Latitude") + xlab("Longitude") +
    coord_map(projection = "mercator", orientation = c(90, 160, 0)) +
    theme(legend.position="top", legend.title = element_blank(),
          panel.background = element_rect(fill = 'aliceblue'),
          panel.grid.major = element_line(colour = "white", size=0.1))

### ---- impatiens-map-group2 ----
## 2. tiger + ESU2 + redSeaTiger
c2Pal <- impPal[c("tiger", "ESU2", "tigerRedSea")]
tmpMap <- subset(impMapThin, consensusESU %in% c("tiger", "ESU2", "tigerRedSea"))

ggplot(tmpMap) + annotation_map(iwp, fill = "gray40", colour = "gray40") +
    geom_point(aes(x = Long2.recenter, y = Lat2, colour = consensusESU), data = tmpMap,
               position = position_dodge(width=1.5), shape=16, size=3) +
    scale_colour_manual(values=c2Pal) + 
    xlim(c(25,220)) + ylim(c(-25, 25)) +  ylab("Latitude") + xlab("Longitude") +
    coord_map(projection = "mercator", orientation = c(90, 160, 0)) +
    theme(legend.position = "top", legend.title = element_blank(),
          panel.background = element_rect(fill = 'aliceblue'),
          panel.grid.major = element_line(colour = "white", size = 0.1))

### ---- impatiens-map-group1 ----
## 3. ESU3 + RedSea + Gracilis + Hawaii + WPac + ESU1
c3Pal <- impPal[c("ESU1", "ESU3", "RedSea", "gracilis", "Hawaii", "Wpac")]
tmpMap <- subset(impMapThin, consensusESU %in% c("ESU1", "ESU3", "RedSea", "gracilis",
                                                 "Hawaii", "Wpac"))
ggplot(tmpMap) + annotation_map(iwp, fill = "gray40", colour = "gray40") +
    geom_point(aes(x = Long2.recenter, y = Lat2, colour = consensusESU), data = tmpMap,
               position = position_jitter(width = 1, height = 1), shape = 16, size=3) +
    scale_colour_manual(values=c3Pal) +
    xlim(c(30,220)) + ylim(c(-27, 29)) + ylab("Latitude") + xlab("Longitude") +
    coord_map(projection = "mercator", orientation = c(90, 160, 0)) +
    theme(legend.position = "top", legend.title = element_blank(),
          panel.grid.major = element_line(colour = "white", size = 0.1),
          panel.background = element_rect(fill = 'aliceblue'))

## pdf(paper = "USr", file = "impatiensMaps.pdf")
## print(gMap)
## print(c1Map)
## print(c2Map)
## print(c3Map)
## dev.off()

