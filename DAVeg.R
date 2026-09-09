
# David Augustine's data:
#cover for 48 polygons from a few days in June 2026:
plot_surveys <- read.csv('Augustine_Plot_Surveys.csv')

#pure shrub and bare soil points and polygons
#pure_ems <- read.csv('Augustine_Pure_EMs.csv')



#shape files for David Augustine's data (although I am unclear how they relate)
LL_pt <- sf::st_zm(sf::st_read("/Users/khuelsma/ARIDNEON/Augustine_VegPoint_2026.shp/VegPOINT_2026.shp"))
pt_vect <- terra::vect(LL_pt)

LL_poly <- sf::st_zm(sf::st_read("/Users/khuelsma/ARIDNEON/Augustine_VegPolygon_2026/VegPolygon2026.shp"))
poly_vect <- terra::vect(LL_poly)
