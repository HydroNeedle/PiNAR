

# Needle generation
# Needle Anatomy in R

# Adrien Heymans
# September 2024

source("~/GitHub/HydroNeedle/PiNAR/R/io_function.R")
setwd("~/GitHub/HydroNeedle/PiNAR/")
# Load cell location 
data <- read.csv("./data/needle_01.csv")

data %>%
  ggplot(aes(X,Y))+
  geom_point()+
  coord_fixed()+theme_classic()
data = rem_dupl_points(data)

pol = alp_polygon(dat = data%>%transmute(x = X, y= Y),a = 50, str_polygon = F)

data = add_euc(data, pol)

# Add primary tissue type
data$type = NA
data$type = ifelse(data$euc < 10, 
                   "epidermis",
                   ifelse(data$euc < 50, "exodermis",
                          ifelse(data$euc < 275, "mesophyll", 
                                 ifelse(data$euc < 305, "endodermis", "stele"))))
data$type = ifelse(data$type == "stele" & data$Area < 250, "vasc", data$type)

data$type = ifelse(data$type == "vasc" & data$X < 1150, "vasc_l", data$type)
data$type = ifelse(data$type == "vasc" & data$X >= 1150, "vasc_r", data$type)

for(v in c("vasc_l","vasc_r")){
  dat= data %>%filter(type == v)
  dat$value = count_(dat)
  filtered_data <- dat %>% filter(value < 220)
  points_inside = inside_ellipse(data = filtered_data)
  
  pole = data%>%filter(X.1 %in% points_inside)
  pole = rem_dupl_points(pole)
  
  alphapoly_pole = alp_polygon(dat = pole%>%transmute(x = X, y= Y), a = 200)
  data$type[Qinside(data,alphapoly_pole)] = 'vascu'
  
}

data = data %>%
  mutate(cambi_curv = -(X-1200)^2/1200+730)

data$type = ifelse(data$Y<data$cambi_curv-15 & data$type == "vascu", "Xylem", data$type)
data$type = ifelse(data$Y>=data$cambi_curv-15 & data$type == "vascu" & data$Y<=data$cambi_curv+15  , "Cambium", data$type)
data$type = ifelse(data$Y>data$cambi_curv+15 & data$type == "vascu", "Phloem", data$type)
data$type = ifelse(data$type == "vasc_l", "stele", data$type)
data$type = ifelse(data$type == "vasc_r", "stele", data$type)


# MAke endodermis

endo = data%>%filter(type == "endodermis")
endo = rem_dupl_points(endo)
alphapoly_endo = alp_polygon(dat= endo%>%transmute(x = X, y = Y), a = 500)

data$endo = "out"
data$endo[Qinside(data,scale_polygon(alphapoly_endo, scale_factor = 1.07))] = 'endo'
# Clear invasive meso cell in the stele
data = data%>%filter(!(type == "mesophyll" & endo == "endo"))

scaled_alphapoly <- scale_polygon(alphapoly_endo, scale_factor = 0.95)
data$endo[Qinside(data,scaled_alphapoly)] = 'true_stele'
# Clear endodermis
data = data%>%filter(!(type == "endodermis" | endo == "endo"))
# New regular
endo_perimeter = st_cast(alphapoly_endo, "LINESTRING")
total_length <- st_length(endo_perimeter)
n=round(total_length /40)
points <- st_line_sample(endo_perimeter, n = round(total_length /40), type = "regular")
pts = st_coordinates(points)[, 1:2]
data = rbind(data%>%dplyr::select(c("X.1", "X","Y", "type")), as_tibble(pts)%>%mutate(type = "endodermis", X.1 = (-nrow(pts)+1):0))

# Make epidermis
epi = data%>%filter(type == "epidermis")
epi = rem_dupl_points(epi)
alphapoly_epi = alp_polygon(dat= epi%>%transmute(x = X, y = Y), a = 800)

data$epi = "out"
data$epi[Qinside(data,alphapoly_epi)] = 'epi'
# Clear epidermis
data = data%>%filter(!(type == "epidermis"))
# New regular
epi_perimeter = st_cast(alphapoly_epi, "LINESTRING")
total_length <- st_length(epi_perimeter)
n=round(total_length /40)
points <- st_line_sample(epi_perimeter, n = round(total_length /20), type = "regular")
pts = st_coordinates(points)[, 1:2]
id_p = (max(data$X.1)+1):(nrow(pts)+max(data$X.1))
data = rbind(data%>%dplyr::select(c("X.1", "X","Y", "type")), as_tibble(pts)%>%mutate(type = "epidermis", X.1 = id_p))

scaled_alphapoly <- scale_polygon(alphapoly_epi, scale_factor = 1.02)
out_perimeter = st_cast(scaled_alphapoly, "LINESTRING")
total_length <- st_length(out_perimeter)
n=round(total_length /40)
points <- st_line_sample(out_perimeter, n = round(total_length /15), type = "regular")
pts = st_coordinates(points)[, 1:2]
id_p = (max(data$X.1)+1):(nrow(pts)+max(data$X.1))
data = rbind(data%>%dplyr::select(c("X.1", "X","Y", "type")), as_tibble(pts)%>%mutate(type = "out", X.1 = id_p))

scaled_alphapoly <- scale_polygon(alphapoly_epi, scale_factor = 1.05)
out_perimeter = st_cast(scaled_alphapoly, "LINESTRING")
total_length <- st_length(out_perimeter)
n=round(total_length /40)
points <- st_line_sample(out_perimeter, n = round(total_length /15), type = "regular")
pts = st_coordinates(points)[, 1:2]
id_p = (max(data$X.1)+1):(nrow(pts)+max(data$X.1))
data = rbind(data%>%dplyr::select(c("X.1", "X","Y", "type")), as_tibble(pts)%>%mutate(type = "out", X.1 = id_p))

# Example usage:
# Assuming `alphapoly_endo` is your original polygon


data%>%
  ggplot(aes(X,Y))+
  geom_point(aes(colour = type))+
  #geom_polygon(aes(x,y), fill = 'blue',alpha = 0.3, data = edge)+
  coord_fixed()+theme_classic()+
  viridis::scale_colour_viridis(discrete = T)

data = rem_dupl_points(data)
data$X.1 = 1:nrow(data)
data = add_euc(data, pol)

df2 = data%>%mutate(x = X/1000, y = Y/1000)

ids = 1:nrow(df2)
df2$id_cell = ids

rs1 = voro(df2)

rs1%>%
  ggplot(aes(x,y))+
  geom_polygon(aes(group = id_cell, fill = type), colour = 'white')+
  coord_fixed()+theme_classic()

rs1 <- rs1%>%filter(type != "out")

# Create a dataframe with the results

po_area = calculate_polygon_area(rs1)


rs0 = left_join(rs1, po_area, "id_cell")

rs0$euc = 0
for(i in unique(rs0$id_cell)){
  pol$euc = sqrt((pol$x-rs0$mx[rs0$id_cell == i][1])^2+(pol$y-rs0$my[rs0$id_cell == i][1])^2)
  rs0$euc[rs0$id_cell == i] = min(pol$euc)  
}

rs0%>%
  ggplot(aes(x,y))+
  geom_polygon(aes(group = id_cell, fill = euc), colour = 'white')+
  coord_fixed()+theme_classic()



rs0%>%
  ggplot(aes(x,-y))+
  geom_polygon(aes(group = id_cell, fill = type), colour = 'white')+
  coord_fixed()+theme_classic()




rs0$sorting <- c(1:nrow(rs0))
nodes <- vertex(rs0)
wally <- nodes[!duplicated(nodes[,c('x1', 'x2', 'y1', 'y2')]),] %>%
  dplyr::select(c(x1, x2, y1, y2))

wally$id_wall <- c(1:nrow(wally))
walls <- wally

nodes <- merge(nodes, walls, by=c("x1", "x2", "y1", "y2"))
nodes <- nodes %>%
  arrange(sorting)

all_cells = rs0%>%
  dplyr::group_by(id_cell)%>%
  dplyr::reframe(type = unique(type),
                   mx = mean(mx),
                   my = mean(my),
                   euc = mean(euc),
                   area = mean(area),
                   perim = mean(perim),
                   circularity = mean(circularity))

# adding the outputs by cell layers
out <- plyr::ddply(all_cells, plyr::.(type), summarise, n_cells=length(type),
                   layer_area = sum(area),
                   cell_area = mean(area)) %>%
  mutate(name = type) %>%
  dplyr::select(-type) %>%
  tidyr::gather(key = "type", value = "value", n_cells, layer_area, cell_area) %>%
  mutate(io = "output")%>%
  dplyr::select(io, everything())
output <-rbind(out, data.frame(io="output", name="all", type="layer_area", value = sum(all_cells$area)))

sim = list()

sim$nodes = nodes
sim$walls = walls
sim$output = output


nodes%>%
  ggplot(aes(x,-y))+
  geom_polygon(aes(group = id_cell, fill = type), colour = 'white')+
  coord_fixed()+theme_classic()


nodes%>%
  ggplot(aes(x,-y))+
  geom_polygon(aes(group = id_cell, fill = type), size = 1, colour = 'white', alpha = 0.7)+
  coord_fixed()+theme_classic()+
  viridis::scale_fill_viridis(discrete = T)


path = "./data/cellsetdata/needle_01.xml"
write_anatomy_xml(sim, path)
