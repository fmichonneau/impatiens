### ---- init-map ----
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
draw_impatiens_map(impatiens_map_data, bg = iwp, palette = impPal,
                   ESUs =  c("WA", "EP", "Gala"), x.lim = c(250, 310),
                   y.lim = c(-25, 25))

### ---- impatiens-map-group2 ----
## 2. tiger + ESU2 + redSeaTiger
draw_impatiens_map(impatiens_map_data, bg = iwp, palette = impPal,
                   ESUs =  c("tiger", "ESU2", "tigerRedSea"),
                   x.lim =c(25, 220), y.lim = c(-25, 25))

### ---- impatiens-map-group1 ----
## 3. ESU3 + RedSea + Gracilis + Hawaii + WPac + ESU1
draw_impatiens_map(impatiens_map_data, bg = iwp, palette = impPal,
                   ESUs = c("ESU1", "ESU3", "RedSea", "gracilis", "Hawaii", "Wpac"),
                   x.lim = c(30, 220), y.lim = c(-27, 29))
