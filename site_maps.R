######################################################################
#
#     Mapping for site locations
#
#     L. McKinley Nevins, laura.nevins@wsu.edu, 31 Aug., 2023
#
##      GNU General Public License, Version 3.0    ###################

library(maps)

#wind river location maps
maps::map(database = 'world', boundary = TRUE, xlim=c(-135,-107.5), ylim=c(20.0,60.0), col = "darkolivegreen4", interior = TRUE, fill=TRUE)
points(-121.957111, 45.819740, pch = 19, col = 'red2', cex = 2)
map.axes()


maps::map(database = 'state', region = c('washington', 'oregon', 'idaho', 'montana'), boundary = TRUE, xlim=c(-126,-115), ylim=c(45.0,49.0), col = "darkolivegreen4", interior = TRUE, fill=TRUE)
points(-121.957111, 45.819740, pch = 19, col = 'red2', cex = 2)
map.axes()

maps::map(database = 'state', region = c('washington', 'oregon', 'idaho', 'montana'), boundary = TRUE, xlim=c(-126,-115), ylim=c(45.0,49.0), col = "darkolivegreen4", interior = TRUE, fill=TRUE)
points(-121.957111, 45.819740, pch = 19, col = 'red2', cex = 2)
map.axes()


maps::map('county', 'washington', boundary = TRUE, xlim=c(-124,-120), ylim=c(45.3,47.0), col = "darkolivegreen4", interior = TRUE, fill=TRUE)
points(-121.957111, 45.819740, pch = 19, col = 'red2', cex = 2)
map.axes()






maps::map(database = 'state', region = c('washington', 'oregon', 'idaho', 'montana', 'california', 'nevada'), boundary = TRUE, xlim=c(-126,-112), ylim=c(40.0,49.0), col = "gray", interior = TRUE, fill=TRUE)
points(-117.166292, 46.729829, pch = 19, col = 'red2', cex = 2)
map.axes()


# Washington, Oregon, and California 
maps::map(database = 'state', region = c('washington', 'oregon', 'california'), boundary = TRUE, xlim=c(-126,-113), ylim=c(30.0,49.0), col = "green4", interior = TRUE, fill=TRUE, bg = "black")



