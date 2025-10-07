

# Generator of needle anatomy

# io_function
library(sf)
library(tidyverse)
`%!in%` <- compose(`!`, `%in%`) 

# Function to check if points are inside the polygon
point_outside <- function(points, polygon) {
  # Convert the points tibble to an sf object
  points_sf <- st_as_sf(points, coords = c("x", "y"), crs = st_crs(polygon))
  
  # Check which points are within the polygon
  points$inside <- st_within(points_sf, polygon, sparse = FALSE)
  
  # Return a logical vector indicating if each point is inside the polygon
  return(points%>%filter(!inside))
}


# Function to check if points are inside the polygon
Qinside <- function(points, polygon) {
  # Convert the points tibble to an sf object
  points_sf <- st_as_sf(points, coords = c("X", "Y"), crs = st_crs(polygon))
  
  # Return a logical vector indicating if each point is inside the polygon
  return(st_within(points_sf, polygon, sparse = FALSE))
}

alp_polygon <- function(dat, a = 100, str_polygon = T){
  # Get the alphahull shape and convert it to dataframe
  my.ashape = alphahull::ashape(x= dat$x, y = dat$y, alpha = a)
  # converted to sf polygons
  a <- data.frame(my.ashape$edges)[,c( 'x1', 'y1', 'x2', 'y2')]
  l <- st_linestring(matrix(as.numeric(a[1,]), ncol=2, byrow = T))
  for(i in 2:nrow(a)){
    l <- c(l, st_linestring(matrix(as.numeric(a[i,]), ncol=2, byrow = T)))
  }
  alphapoly <- st_sf(geom = st_sfc(l), crs = 2056) %>% st_polygonize() %>% st_collection_extract()
  pol = as.data.frame(matrix(unlist(alphapoly$geom[[1]]),ncol = 2))
  colnames(pol) = c("x", "y")
  
  if(str_polygon){
    return(alphapoly)
  }else{
    return(pol)
  }
}

add_euc <- function(data, pol){
  data$euc = 0
  for(i in data$X.1){
    pol$euc = sqrt((pol$x-data$X[i])^2+(pol$y-data$Y[i])^2)
    data$euc[i] = min(pol$euc)  
  }
  return(data)
}

count_ <- function(dat) {
  sapply(1:nrow(dat), function(i) {
    sum(sqrt((dat$X - dat$X[i])^2 + (dat$Y - dat$Y[i])^2))/nrow(dat)
  })
}


rem_dupl_points <- function(df){
  df$id_point <- paste0(round(df$X, 2),";",round(df$Y,2))
  df = df%>%
    mutate(dup = duplicated(id_point))%>%
    filter(!dup)
  return(df)
}



# Function to check if points lie within the ellipse
inside_ellipse <- function(data, level = 0.95) {
  # Calculate the covariance matrix
  cov_matrix <- cov(data[, c("X", "Y")])
  # Calculate the mean of X and Y
  center <- colMeans(data[, c("X", "Y")])
  # Calculate the Mahalanobis distance for each point
  distances <- mahalanobis(data[, c("X", "Y")], center, cov_matrix)
  # Chi-square critical value for the given level (e.g., 95% confidence level)
  chi_sq_threshold <- qchisq(level, df = 2)
  # Return a logical vector indicating whether each point is within the ellipse
  id = data$X.1[which(distances <= chi_sq_threshold)]
  return(id)
}



# Function to scale down a polygon towards its center of mass
scale_polygon <- function(polygon, scale_factor = 0.8) {
  # Calculate the centroid (center of mass) of the polygon
  centroid <- st_centroid(polygon)
  
  # Extract the coordinates of the polygon and the centroid
  polygon_coords <- st_coordinates(polygon)[, 1:2]
  centroid_coords <- st_coordinates(centroid)[1, 1:2]
  
  x_c = (polygon_coords[,1]-centroid_coords[1])*scale_factor +centroid_coords[1]
  y_c = (polygon_coords[,2]-centroid_coords[2])*scale_factor +centroid_coords[2]
  
  
  
  # Create a new scaled polygon
  scaled_polygon <- st_polygon(list(matrix( c(x_c, y = y_c), ncol = 2)))
  
  # Return the scaled polygon as an sf object with the same CRS as the original
  return(st_sfc(scaled_polygon, crs = st_crs(polygon)))
}

voro <- function(df2){
  
  # Get the voronio data
  vtess <- deldir::deldir(df2$x, df2$y, digits = 6)
  if(is.null(vtess)){return(NULL)}
  
  rs <- vtess$dirsgs[vtess$dirsgs$ind1 %in% ids |
                       vtess$dirsgs$ind2 %in% ids,]
  
  # Get the cooridnates for every nodes in the voronoi
  rs <- rs %>% arrange(ind1)
  rs2 <- data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind1)
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind1))
  rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind2))
  rs2 <- rbind(rs2, data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind2))
  
  # rs2%>%
  #   ggplot(aes(x,y))+
  #   geom_point(aes(colour = id_cell))+
  #   coord_fixed()+theme_classic()+
  #   viridis::scale_colour_viridis()
  # 
  # length(unique(rs2$id_cell))
  # length(unique(df2$id_cell))
  
  rs2 <- left_join(rs2, df2[,c("id_cell", "type")], by="id_cell")
  
  
  rs1 <- rs2 %>%
    dplyr::group_by(id_cell) %>%
    dplyr::mutate(my = mean(y),
                  mx = mean(x),
                  atan = atan2(y-my, x - mx)) %>%
    dplyr::arrange(id_cell, atan)%>%
    ungroup()
  
  rs1$id_point <- paste0(rs1$x,";",rs1$y)
  return(rs1)
}

calculate_polygon_area <- function(df) {
  # List unique polygon IDs
  unique_ids <- unique(df$id_cell)
  
  # Calculate the area for each polygon using sapply
  areas <- sapply(unique_ids, function(id) {
    # Subset the dataframe for the current polygon
    polygon_df <- df %>% filter(id_cell == id)
    # Create the polygon
    polygon <- st_polygon(list(matrix(c(polygon_df$x ,polygon_df$x[1], polygon_df$y,polygon_df$y[1]), ncol = 2)))
    
    # Calculate the area
    area <- st_area(st_sfc(polygon))
    
    # Calculate the perimeter
    perimeter <- st_length(st_cast(st_sfc(polygon), "LINESTRING"))
    
    # Calculate the circularity
    circularity <- (4 * pi * as.numeric(area)) / (as.numeric(perimeter) ^ 2)
    
    # Return a list with area and circularity
    return(c(area = as.numeric(area), perimeter = as.numeric(perimeter), circularity = circularity))
  })
  result <- data.frame(id_cell = unique_ids, area = areas["area",], perim = areas["perimeter",], circularity = areas["circularity",])
  
  
  return(result)
}


vertex <- function(rs1){
  nodes <- rs1 %>%
    mutate(id_point= paste0(rs1$x,";",rs1$y))%>%
    group_by(id_cell) %>%
    filter(!duplicated(id_point))%>%
    dplyr::mutate(xx = c(x[-1],x[1])) %>%
    dplyr::mutate(yy = c(y[-1],y[1]))
  
  nodes <- nodes %>%
    ungroup() %>%
    mutate(vertical = ifelse(x == xx, "true", "false")) %>%
    mutate(x1 = ifelse(x > xx, x, xx)) %>%
    mutate(x2 = ifelse(x > xx, xx, x)) %>%
    mutate(y1 = ifelse(x > xx, y, yy)) %>%
    mutate(y2 = ifelse(x > xx, yy, y)) %>%
    # If wall is perfectly vertical
    mutate(y1 = ifelse(x == xx,
                       ifelse(y > yy, yy, y), y1)) %>%
    mutate(y2 = ifelse(x == xx,
                       ifelse(y > yy, y, yy), y2)) %>%
    mutate(wall_length = sqrt((x2-x1)^2 + (y2-y1)^2)) %>%
    mutate(wall_length2 = sqrt((xx-x)^2 + (yy-y)^2),
           slope = (y2-y1)/(x2-x1),
           intercept = y1 - slope*x1)
  return(nodes)
}



write_anatomy_xml <- function(sim = NULL, path = NULL){
  
  if(is.null(sim)) warning("No simulation found. Please input a GRANAR simulation")
  if(is.null(path)) warning("No path found to save the XML file")
  
  
  if(length(sim$walls$x3) > 0){
    nodal <- sim$walls_nodes
  }else{
    nodal <- sim$nodes
  }

  cellgroups <- data.frame(id_group = c(1, 2, 3, 3, 4, 5, 13, 16, 12, 11, 4, 4, 11, 13, 11),
                           type = c("exodermis", "epidermis", "endodermis", "passage_cell",  "mesophyll",
                                    "stele", "Xylem", "pericycle", "companion_cell", "Phloem",
                                    "inter_cellular_space", "aerenchyma", "Cambium", "metaxylem", "wax"))
  
  xml <- '<?xml version="1.0" encoding="utf-8"?>\n'
  xml <- paste0(xml, '<granardata>\n')
  
  # Write the Metadata
  xml <- paste0(xml, '\t<metadata>\n')
  xml <- paste0(xml, '\t\t<parameters>\n')
  xml <- paste0(xml,paste0('\t\t\t<parameter io="',sim$output$io,'" ',
                           'name="',sim$output$name,'" ',
                           'type="',sim$output$type,'" ',
                           'value="',sim$output$value,'"/>\n', collapse = ""))
  xml <- paste0(xml, '\t\t</parameters>\n')
  xml <- paste0(xml, '\t</metadata>\n')
  
  # Write the cells information
  xml <- paste0(xml, '\t<cells count="',length(unique(nodal$id_cell)),'">\n')
  
  nodes_data <- merge(nodal, cellgroups, by="type")
  
  temp_wall <- plyr::ddply(nodes_data, plyr::.(id_cell, id_group), summarise, walls = paste0('\t\t\t\t<wall id="',
                                                                                             paste(id_wall-1, collapse='"/>\n\t\t\t\t<wall id="'),
                                                                                             '"/>\n'))
  xml <- paste0(xml, paste0('\t\t<cell id="',temp_wall$id_cell-1, '" group="', temp_wall$id_group, '" truncated="false" >\n',
                            '\t\t\t<walls>\n', temp_wall$walls, '\t\t\t</walls>\n',
                            '\t\t</cell>\n', collapse=""))
  xml <- paste0(xml, '\t</cells>\n')
  
  
  # Write the walls information
  xml <- paste0(xml, '\t<walls count="',length(unique(nodes_data$id_wall)),'">\n')
  
  walls <- sim$walls%>%
    dplyr::select(starts_with("x"), starts_with("y"))
  col_nam <- colnames(walls)
  
  substr1(col_nam[nchar(col_nam) == 2], 1.5) <- "0"
  col_nam <- paste0(substr1(col_nam, -1), "_", substr1(col_nam, 1))
  colnames(walls) <- col_nam
  
  sorted_name <- sort(col_nam)
  walls <- walls%>%dplyr::select(sorted_name)
  
  N <- max(readr::parse_number(col_nam))
  begin <- tibble(tag1 = '\t\t<wall id="',
                  id_wall = sim$walls$id_wall-1,
                  tag2 = '" group="0" edgewall="false" >\n\t\t\t<points>\n')
  middle <- tibble(tag_x1 = '\t\t\t\t<point x="',
                   x1 = round(walls[,sorted_name[1]],6),
                   tag_y1 = '" y="',
                   y1 = round(walls[,sorted_name[2]],6),
                   tag_end1 = '"/>\n')
  for(k in 2:N){
    h <- k*2-1 # odd number
    tmp_coord <- walls[,sorted_name[c(h,h+1)]]
    tmp_middle <- tibble(tag_x = '\t\t\t\t<point x="',
                         x = round(tmp_coord[,1],6),
                         tag_y = '" y="',
                         y = round(tmp_coord[,2],6),
                         tag_end = '"/>\n')
    tmp_col_name <- colnames(tmp_middle)
    colnames(tmp_middle) <- paste0(t(tmp_col_name), k)
    middle <- cbind(middle, tmp_middle)
  }
  taged_walls <- cbind(begin,middle)%>%
    mutate(tag_ending = '\t\t\t</points>\n\t\t</wall>\n')
  xml <- paste0(xml, paste0(t(taged_walls), collapse = ""))
  xml <- paste0(xml, '\t</walls>\n')
  xml <- stringr::str_remove_all(xml, '\t\t\t\t<point x=\"NA\" y=\"NA\"/>\n')
  
  # Write the cell group informations
  print(cellgroups)
  xml <- paste0(xml, '\t<groups>\n')
  xml <- paste0(xml, '\t\t<cellgroups>\n')
  for(i in c(1:nrow(cellgroups))){
    xml <- paste0(xml, '\t\t\t<group id="',cellgroups$id_group[i],'" name="',cellgroups$type[i],'" />\n')
  }
  xml <- paste0(xml, '\t\t</cellgroups>\n')
  xml <- paste0(xml, '\t\t<wallgroups>\n')
  xml <- paste0(xml, '\t\t\t<group id="0" name="unassigned" />\n')
  xml <- paste0(xml, '\t\t</wallgroups>\n')
  xml <- paste0(xml, '\t</groups>\n')
  
  xml <- paste0(xml, '</granardata>')
  
  if(!is.null(path)){
    cat(xml, file = path)
    return(TRUE)
  }else{
    return(xml)
  }
  
  
}


#' @title remove one str
#'
#' @param x string
#' @param y pattern
#' @keywords string
#' @export
#'
substr1 <- function(x,y) {
  z <- sapply(strsplit(as.character(x),''),function(w) paste(stats::na.omit(w[y]),collapse=''))
  dim(z) <- dim(x)
  return(z) }

#' @title remove one str
#'
#' @param x string
#' @param y pattern
#' @param value something
#' @keywords string
#' @export
#'

`substr1<-` <- function(x,y,value) {
  names(y) <- c(value,rep('',length(y)-length(value)))
  z <- sapply(strsplit(as.character(x),''),function(w) {
    v <- seq(w)
    names(v) <- w
    paste(names(sort(c(y,v[setdiff(v,y)]))),collapse='') })
  dim(z) <- dim(x)
  return(z) }