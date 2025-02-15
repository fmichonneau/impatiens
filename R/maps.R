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
    if (n == 0) return(z) else {
        for (i in 1:n) {
            for (j in 1:n) z[i, j] <- gcd.hf(long1 = latlongs[i, 1],
                                             lat1 = latlongs[i, 2],
                                             long2 = latlongs[j, 1],
                                             lat2 = latlongs[j,2])
        }
        z <- as.dist(z)
    }
    return(z)
}

thinCoords <- function(coords=impMap, min_dist = 100) {
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
    whichDup <- lbl[which(d < min_dist), ]
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

get_impatiens_map_data <- function(impDB) {
    impMap <- impDB[nzchar(impDB$consensusESU), ]
    center <- 200
    impMap$Long.recenter <- ifelse(impMap$Long < center - 180, impMap$Long + 360, impMap$Long)
    impMapThin <- lapply(unique(impMap$consensusESU), function(x) {
                             impMap_tmp <- subset(impMap, consensusESU == x)
                             thinCoords(impMap_tmp)
                         })
    impMapThin <- do.call("rbind", impMapThin)
    impMapThin <- impMapThin[!duplicated(paste(impMapThin$consensusESU, impMapThin$Lat2, impMapThin$Long2)), ]
    impMapThin
}

draw_impatiens_map <- function(map_data, bg, palette, ESUs, x.lim, y.lim) {
    tmpMap <- subset(map_data, consensusESU %in% ESUs)
    pal <- palette[ESUs]

    g <- ggplot(tmpMap) + annotation_map(bg, fill = "gray40", colour = "gray40") +
           geom_point(aes(x = Long2.recenter, y = Lat2, colour = consensusESU),
                      position = position_dodge(width = 1.5), shape = 16, size = 3) +
           scale_colour_manual(values = pal) +
           xlim(x.lim) + ylim(y.lim) +  ylab("Latitude") + xlab("Longitude") +
           coord_map(projection = "mercator", orientation = c(90, 160, 0)) +
           theme(legend.position = "top", legend.title = element_blank(),
                 panel.background = element_rect(fill = 'aliceblue'),
                 panel.grid.major = element_line(colour = "white", size = 0.1))
    g
}
