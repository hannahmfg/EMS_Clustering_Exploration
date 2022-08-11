##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                          EMS PROJECT EXPLORATIONS                        ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)
library(leaflet)
library(ClustGeo)
library(cluster)
library(factoextra)
library(osrm)
library(osrmr)
library(ggplot2)
library(GDAtools)
library(lightslateblue404)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                   Set Up                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat_nearest = readRDS("dat_nearest.rds")
data(wyoming_census)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        Travel Time Distance Matrices                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dist_all = readRDS("dist_mats.rds")
alb_travel = dist_all$Albany_County_distmat
big_h_travel = dist_all$Big_Horn_County_distmat
camp_travel = dist_all$Campbell_County_distmat
carb_travel = dist_all$Carbon_County_distmat
conv_travel = dist_all$Converse_County_distmat
crook_travel = dist_all$Converse_County_distmat
frem_travel = dist_all$Fremont_County_distmat
gosh_travel = dist_all$Goshen_County_distmat
hsp_travel = dist_all$Hot_Springs_County_distmat
john_travel = dist_all$Johnson_County_distmat
lar_travel = dist_all$Laramie_County_distmat
linc_travel = dist_all$Lincoln_County_distmat
natr_travel = dist_all$Natrona_County_distmat
niob_travel = dist_all$Niobrara_County_distmat
park_travel = dist_all$Park_County_distmat
platte_travel = dist_all$Platte_County_distmat
sher_travel = dist_all$Sheridan_County_distmat
subl_travel = dist_all$Sublette_County_distmat
sweet_travel = dist_all$Sweetwater_County_distmat
teton_travel = dist_all$Teton_County_distmat
uinta_travel = dist_all$Uinta_County_distmat
wash_travel = dist_all$Washakie_County_distmat
west_travel = dist_all$Weston_County_distmat

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        Maximal Number of Allocations                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#...............................................................................
#                                                                              .
#  (NASEMSO, 2020): With no additional information regarding the breakdown of  .
#  available EMS vehicles in the state of Wyoming and assuming the latter      .
#  statements hold true, Wyoming has about 142 ambulances to allocate across   .
#  the state. By proportionality, the break down by county is as follows:      .
# 
#   A_i = County_i/State*142   [state = 578579; amb = 142]
#
#   Albany County:     38880/state*amb = 9.542275  = 9
#   Big Horn County:   11882/state*amb = 2.916186  = 3
#   Campbell County:   46341/state*amb = 11.37342  = 11
#   Carbon County:     15247/state*amb = 3.742054  = 4
#   Converse County:   13921/state*amb = 3.416616  = 3
#   Crook County:       7584/state*amb = 1.861333  = 2
#   Fremont County:    39261/state*amb = 9.635784  = 10
#   Goshen County:     13342/state*amb = 3.274512  = 3
#   Hot Springs County: 4607/state*amb = 1.130691  = 2
#   Johnson County:     8487/state*amb = 2.082955  = 2
#   Laramie County:    99500/state*amb = 24.42017  = 24
#   Lincoln County:    19274/state*amb = 4.730396  = 5
#   Natrona County:    79858/state*amb = 19.59946  = 20
#   Niobrara County:    2422/state*amb = 0.5944288 = 1
#   Park County:       29148/state*amb = 7.153761  = 7
#   Platte County:      8582/state*amb = 2.106271  = 2
#   Sheridan County:   30485/state*amb = 7.4819    = 7
#   Sublette County:    9880/state*amb = 2.414837  = 2
#   Sweetwater County: 42343/state*amb = 10.3922   = 10
#   Teton County:      23464/state*amb =  5.758743 = 6
#   Uinta County:      20226/state*amb =  4.964045 = 5
#   Washakie County:    8027/state*amb =  1.970058 = 2
#   Weston County:      6927/state*amb =  1.70086  = 2
#
#...............................................................................

#...............................................................................
#                                                                              .
#  Some of these smaller counties can probably be included with other regions  .
#  of service.                                                                 .
#                                                                              .
#...............................................................................

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Albany County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 9 to allocate #
alb_dist = dat_nearest[county == "Albany County", c(3,47:49)]
pop_wt_alb = as.vector(unlist(wyoming_census[county == "Albany County", 8]))
albany = clusters_4_dayz(9,alb_dist, alb_travel, pop_wt_alb)
# I used alpha = 0.80
albany_map = maps_for_dayz(alb_dist,albany)
albany_map

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Big Horn County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3 to allocate #
bh_dist = dat_nearest[county == "Big Horn County",c(3,47:49)]
pop_wt_bh = as.vector(unlist(wyoming_census[county == "Big Horn County", 8]))
big_horn = clusters_4_dayz(3,bh_dist, big_h_travel, pop_wt_bh)
# I used 0.95.
big_horn_map = maps_for_dayz(bh_dist,big_horn)
big_horn_map

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Campbell County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 11 to allocate #
camp_dist = dat_nearest[county == "Campbell County",c(3,47:49)]
pop_wt_camp = as.vector(unlist(wyoming_census[county == "Campbell County", 8]))
campbell = clusters_4_dayz(11,camp_dist, camp_travel, pop_wt_camp)
# I used 0.75 #
campbell_map = maps_for_dayz(camp_dist,campbell)
campbell_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Carbon County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4 to allocate #
carb_dist = dat_nearest[county == "Carbon County", c(3,47:49)]
pop_wt_carb = as.vector(unlist(wyoming_census[county == "Carbon County", 8]))
carbon = clusters_4_dayz(4,carb_dist, carb_travel, pop_wt_carb)
# I used 0.95
carbon_map = maps_for_dayz(carb_dist,carbon)
carbon_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Converse County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3 to allocate #
conv_dist = dat_nearest[county == "Converse County", c(3,47:49)]
pop_wt_conv = as.vector(unlist(wyoming_census[county == "Converse County", 8]))
converse = clusters_4_dayz(3,conv_dist, conv_travel, pop_wt_conv)
# I used 0.90 #
converse_map = maps_for_dayz(conv_dist,converse)
converse_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Crook County                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 to allocate #
crook_dist = dat_nearest[county == "Crook County", c(3,47:49)]
options(osrm.server =  "http://0.0.0.0:5000/", osrm.profile = "car") 
crook_travel = osrmTable(loc = crook_dist[, .(ID,Near_Lon,Near_Lat)])$durations
pop_wt_crook = as.vector(unlist(wyoming_census[county == "Crook County", 8]))
crook = clusters_4_dayz(2,crook_dist, crook_travel, pop_wt_crook)
# This one is a little weird compared to the others. I chose 0.70 as the value
# For the mixing parameter.
crook_map = maps_for_dayz(crook_dist,crook)
crook_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               Fremont County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 10 to allocate #
frem_dist = dat_nearest[county == "Fremont County", c(3,47:49)]
pop_wt_frem = as.vector(unlist(wyoming_census[county == "Fremont County", 8]))
fremont = clusters_4_dayz(10,frem_dist, frem_travel, pop_wt_frem)
# I used 0.95
fremont_map = maps_for_dayz(frem_dist,fremont)
fremont_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Goshen County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3 to allocate #
gosh_dist = dat_nearest[county == "Goshen County", c(3,47:49)]
pop_wt_gosh = as.vector(unlist(wyoming_census[county == "Goshen County", 8]))
goshen = clusters_4_dayz(3,gosh_dist, gosh_travel, pop_wt_gosh)
# Like the other really small counties, this is a little weird.
# I chose 0.70 for this one.
goshen_map = maps_for_dayz(gosh_dist, goshen)
goshen_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                             Hot Springs County                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 to allocate #
hs_dist = dat_nearest[county == "Hot Springs County", c(3, 47:49)]
pop_wt_hs = as.vector(unlist(wyoming_census[county == "Hot Springs County", 8]))
hot_springs = clusters_4_dayz(2,hs_dist, hsp_travel, pop_wt_hs)
# I used alpha = 0.80
hot_springs_map = maps_for_dayz(hs_dist,hot_springs)
hot_springs_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               Johnson County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 to allocate #
john_dist = dat_nearest[county == "Johnson County", c(3, 47:49)]
pop_wt_j = as.vector(unlist(wyoming_census[county == "Johnson County", 8]))
johnson = clusters_4_dayz(2,john_dist, john_travel, pop_wt_j)
# I used alpha = 0.95
johnson_map = maps_for_dayz(john_dist,johnson)
johnson_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               Laramie County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 24 to allocate #
lar_dist = dat_nearest[county == "Laramie County", c(3,47:49)]
pop_wt_lar = as.vector(unlist(wyoming_census[county == "Laramie County", 8]))
laramie = clusters_4_dayz(24,lar_dist, lar_travel, pop_wt_lar)
# I used alpha = 0.90
laramie_c_map = maps_for_dayz(lar_dist,laramie)
laramie_c_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               Lincoln County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5 to allocate #
linc_dist = dat_nearest[county == "Lincoln County", c(3, 47:49)]
pop_wt_lc = as.vector(unlist(wyoming_census[county == "Lincoln County", 8]))
lincoln = clusters_4_dayz(5,linc_dist, linc_travel, pop_wt_lc)
# I used alpha = 0.75
lincoln_map = maps_for_dayz(linc_dist,lincoln)
lincoln_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               Natrona County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 20 to allocate #
nat_dist = dat_nearest[county == "Natrona County", c(3,47:49)]
pop_wt_nat = as.vector(unlist(wyoming_census[county == "Natrona County", 8]))
natrona = clusters_4_dayz(20, nat_dist, natr_travel, pop_wt_nat)
# I used alpha = 0.75
natrona_map = maps_for_dayz(nat_dist,natrona)
natrona_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Niobrara County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1 to allocate- this may be a little silly! #
nio_dist = dat_nearest[county == "Niobrara County", c(3,47:49)]
pop_wt_nio = as.vector(unlist(wyoming_census[county == "Niobrara County", 8]))
niobrara = clusters_4_dayz(1,nio_dist, niob_travel, pop_wt_nio)
# I used a value of 1. I think the plot doesn't show because it's not really 
# clustering, if we desire only one cluster.
niobrara_map = maps_for_dayz(nio_dist,niobrara)
niobrara_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Park County                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7 to allocate #
park_dist = dat_nearest[county == "Park County", c(3,47:49)]
pop_wt_park = as.vector(unlist(wyoming_census[county == "Park County", 8]))
park = clusters_4_dayz(7, park_dist, park_travel, pop_wt_park)
# I used alpha = 0.80
park_map = maps_for_dayz(park_dist, park)
park_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Platte County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 to allocate #
platte_dist = dat_nearest[county == "Platte County", c(3, 47:49)]
pop_wt_pl = as.vector(unlist(wyoming_census[county == "Platte County", 8]))
platte = clusters_4_dayz(2, platte_dist, platte_travel, pop_wt_pl)
# I used alpha = 0.70
platte_map = maps_for_dayz(platte_dist,platte)
platte_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Sheridan County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 7 to allocate #
sher_dist = dat_nearest[county == "Sheridan County", c(3,47:49)]
pop_wt_sher = as.vector(unlist(wyoming_census[county == "Sheridan County", 8]))
sheridan = clusters_4_dayz(7, sher_dist, sher_travel, pop_wt_sher)
# I used alpha = 0.70
sheridan_map = maps_for_dayz(sher_dist,sheridan)
sheridan_map

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Sublette County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 to allocate #
subl_dist = dat_nearest[county == "Sublette County", c(3,47:49)]
pop_wt_subl = as.vector(unlist(wyoming_census[county == "Sublette County", 8]))
sublette = clusters_4_dayz(2, subl_dist, subl_travel, pop_wt_subl)
# I used alpha = 0.80
sublette_map = maps_for_dayz(subl_dist, sublette)
sublette_map

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Sweetwater County                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 10 to allocate #
sweet_dist =  dat_nearest[county == "Sweetwater County", c(3,47:49)]
pop_wt_sw = as.vector(unlist(wyoming_census[county == "Sweetwater County", 8]))
sweetwater = clusters_4_dayz(10, sweet_dist, sweet_travel, pop_wt_sw)
# I used alpha = 0.80
sweetwater_map = maps_for_dayz(sweet_dist,sweetwater)
sweetwater_map
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Teton County                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 6 to allocate #
teton_dist = dat_nearest[county == "Teton County",c(3,47:49)]
pop_wt_tet = as.vector(unlist(wyoming_census[county == "Teton County", 8]))
teton = clusters_4_dayz(6, teton_dist, teton_travel, pop_wt_tet)
# Another weird one. I used alpha of 0.60.
teton_map = maps_for_dayz(teton_dist,teton)
teton_map

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Uinta County                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5 to allocate #
uinta_dist = dat_nearest[county == "Uinta County",c(3,47:49)]
pop_wt_uin = as.vector(unlist(wyoming_census[county == "Uinta County", 8]))
uinta = clusters_4_dayz(5, uinta_dist, uinta_travel, pop_wt_uin)
# I used alpha = 0.80
uinta_map = maps_for_dayz(uinta_dist, uinta)
uinta_map

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Washakie County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 to allocate #
wash_dist = dat_nearest[county == "Washakie County", c(3,47:49)]
pop_wt_wash= as.vector(unlist(wyoming_census[county == "Washakie County", 8]))
washakie = clusters_4_dayz(2, wash_dist, wash_travel, pop_wt_wash)
# I used alpha = 1
washakie_map = maps_for_dayz(wash_dist, washakie)
washakie_map

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Weston County                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2 to allocate #
west_dist = dat_nearest[county == "Weston County",c(3,47:49)]
pop_wt_west = as.vector(unlist(wyoming_census[county == "Weston County", 8]))
weston = clusters_4_dayz(2, west_dist, west_travel, pop_wt_west)
# Another weird one. I used 0.95.
weston_map = maps_for_dayz(west_dist,weston)
weston_map

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                      Concise Maps and Medoid Locations                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# all_maps = readRDS("All County Maps.rds")
all_maps = list(Albany_County = albany_map, Big_Horn_County = big_horn_map, 
                Campbell_County = campbell_map, Carbon_County = carbon_map,
                Converse_County = converse_map, Crook_County = crook_map, 
                Fremont_Map = fremont_map, Goshen_County = goshen_map,
                HotSprings_County = hot_springs_map, Johnson_County = johnson_map,
                Laramie_County = laramie_c_map, Lincoln_County = lincoln_map, 
                Natrona_County = natrona_map, Niobrara_County = niobrara_map, 
                Park_County = park_map, Platte_County = platte_map, 
                Sheridan_County = sheridan_map, Sublette_County = sublette_map, 
                Sweetwater_County = sweetwater_map,Teton_County = teton_map, 
                Washakie_County = washakie_map, Weston_County = weston_map)
#saveRDS(all_maps, file = "All County Maps.rds")

# all_medoids = readRDS("All County Medoid Coords.rds")
all_medoids_locs = list(Albany_County = albany$medoids, Big_Horn_County = big_horn$medoids, 
                        Campbell_County = campbell$medoids, Carbon_County = carbon$medoids,
                        Converse_County = converse$medoids, Crook_County = crook$medoids, 
                        Fremont_Map = fremont$medoids, Goshen_County = goshen$medoids,
                        HotSprings_County = hot_springs$medoids, Johnson_County = johnson$medoids,
                        Laramie_County = laramie$medoids, Lincoln_County = lincoln$medoids, 
                        Natrona_County = natrona$medoids, Niobrara_County = niobrara$medoids, 
                        Park_County = park$medoids, Platte_County = platte$medoids, 
                        Sheridan_County = sheridan$medoids, Sublette_County = sublette$medoids, 
                        Sweetwater_County = sweetwater$medoids,Teton_County = teton$medoids, 
                        Washakie_County = washakie$medoids, Weston_County = weston$medoids)

#saveRDS(all_medoids_locs, file = "All County Medoid Coords.rds")
