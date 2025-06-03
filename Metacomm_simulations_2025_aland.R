#CODE WRITTEN AND CONFIGURATED BY JULIA SIDESTÅL SPRING 2025



library(tidyverse)
library(mecobane)

#######clean data set- for baseline extinction #######
#getting patches with coordinates, FULL DATA
patch <- read_csv("patch_data.csv")

#getting data from patches
#cleaning up bcs of wrong numbers
patchdata <- read_csv("TotalPatchDat.csv") |> 
  mutate(MelC = ifelse(MelC == 0, 0, 1),
         PhoP = ifelse(PhoP == 0, 0, 1),
         CotP = ifelse(CotP == 0, 0, 1)
  )

#joining data with the right patches
completepatchdata <- inner_join(patch, patchdata, by = join_by(Patch)) |> 
  mutate(PLpresence = ifelse(PL == 0, 0, 1),
         VSpresence = ifelse(VS == 0, 0, 1))

###### getting baseline extinction rate for the species #######

#function for baseline extinction
baselineExtinct<- function(data, species) {
  data |> #first part is to tidy data
    arrange(Year) |> 
    select({{species}}, "Year", "Patch") |> 
    group_by(Patch) |> 
    #making cols with the numbers that come after
    mutate(Species_next = lead({{species}}),
           Year_next = lead(Year)) |> 
    ungroup() |> 
    #filtering for where the species is present
    filter(!is.na(Species_next), {{species}} == 1) |> 
    mutate(extinct = ifelse(Species_next == 0, 1, 0)) |> 
    group_by(Year) |> 
    summarise   (n = n(),
                 n_extinct = sum(extinct),
                 extinct_rate = n_extinct / n) |> 
    ungroup() |> 
    summarise( m = mean(extinct_rate)) |> 
    pull(m)
}
#for non basal species
######MELC#####
MelCPLVS <- completepatchdata |> 
  arrange(Year) |> 
  select("Year", "Patch", "MelC", "PLpresence", "VSpresence") |> 
  group_by(Patch) |> 
  filter(MelC == 1, PLpresence == 1, VSpresence == 1) |> 
  ungroup()

MelCbase <- semi_join(completepatchdata, MelCPLVS, by = join_by(Patch)) |> 
  select("Year", "Patch", "MelC", "PLpresence", "VSpresence") |> 
  group_by(Patch) |> 
  mutate(Yes = ifelse(MelC == 1 & PLpresence == 1 & VSpresence == 1, 1,
                      ifelse(MelC == 0 & PLpresence == 0 & VSpresence == 0, 0, NA)),
         Yes_next = lead(Yes),
         Year_next = lead(Year)) |> 
  ungroup() |> 
  filter(!is.na(Yes_next), Yes == 1) |> 
  mutate(extinct = ifelse(Yes_next == 0, 1, 0)) |> 
  group_by(Year) |> 
  summarise   (n = n(),
               n_extinct = sum(extinct),
               extinct_rate = n_extinct / n) |> 
  ungroup() |> 
  summarise( m = mean(extinct_rate)) |> 
  pull(m)

#######COTP######
CotPMelC <- completepatchdata |> 
  arrange(Year) |> 
  select("Year", "Patch", "MelC", "CotP") |> 
  group_by(Patch) |> 
  filter(MelC == 1, CotP == 1) |> 
  ungroup()

CotPbase <- semi_join(completepatchdata, CotPMelC, by = join_by(Patch)) |> 
  select("Year", "Patch", "MelC", "CotP") |> 
  group_by(Patch) |> 
  mutate(Yes = ifelse(CotP== 1 & MelC == 1, 1,
                      ifelse(CotP == 0 & MelC == 0, 0, NA)),
         Yes_next = lead(Yes),
         Year_next = lead(Year)) |> 
  ungroup() |> 
  filter(!is.na(Yes_next), Yes == 1) |> 
  mutate(extinct = ifelse(Yes_next == 0, 1, 0)) |> 
  group_by(Year) |> 
  summarise   (n = n(),
               n_extinct = sum(extinct),
               extinct_rate = n_extinct / n) |> 
  ungroup() |> 
  summarise( m = mean(extinct_rate)) |> 
  pull(m)
#######PHOP######
PhoPPL <- completepatchdata |> 
  arrange(Year) |> 
  select("Year", "Patch", "PhoP", "PLpresence") |> 
  group_by(Patch) |> 
  filter(PhoP== 1, PLpresence == 1) |> 
  ungroup()

PhoPbase <- semi_join(completepatchdata, PhoPPL, by = join_by(Patch)) |> 
  select("Year", "Patch", "PhoP", "PLpresence") |> 
  group_by(Patch) |> 
  mutate(Yes = ifelse(PhoP== 1 & PLpresence == 1, 1,
                      ifelse(PhoP == 0 & PLpresence == 0, 0, NA)),
         Yes_next = lead(Yes),
         Year_next = lead(Year)) |> 
  ungroup() |> 
  filter(!is.na(Yes_next), Yes == 1) |> 
  mutate(extinct = ifelse(Yes_next == 0, 1, 0)) |> 
  group_by(Year) |> 
  summarise   (n = n(),
               n_extinct = sum(extinct),
               extinct_rate = n_extinct / n) |> 
  ungroup() |> 
  summarise( m = mean(extinct_rate)) |> 
  pull(m)

#for every year and year +1, where butterfly plant1 and plant2 
#do same simulations above but with dependencies
#do same for fungal and wasp


MelCbase#0.0112
PhoPbase #0.0389
CotPbase #0
baselineExtinct(completepatchdatahdata, PLpresence) #0.0342
baselineExtinct(completepatchdata, VSpresence) #0.144

piNumbers <- tibble(PL = baselineExtinct(completepatchdata, PLpresence),
                    VS = baselineExtinct(completepatchdata, VSpresence),
                    MelC = MelCbase,
                    PhoP = PhoPbase,
                    CotP = CotPbase)

################################################################################
# Calculations for probability tables - for MelC 

#p(-c|-s,l)
MelCbase + (1- MelCbase)*0.5 #0.505


#p(-c|s,-l) same as above


################################################################################
#speciesinfo and coords

patchcoordinates <-
  read_csv("patch_data.csv") |>
  select(X, Y) |>
  drop_na() |> # removing NAs in case it affects the simulation
  rename("dim1" = X, "dim2" = Y)

# making data set with pi-numbers
speciesinfo <- tibble(
  species = c("PL", "VS", "Butterfly", "PLmildew", "Parasitoidwasp"),
  pi = unlist(piNumbers),
  xi = 1,
  kernel = "Exponential"
)
#interactions table
edgeList <- tibble(
  consumer = c("Butterfly", "Butterfly", "Parasitoidwasp", "PLmildew"),
  resource = c("PL", "VS", "Butterfly", "PL")) 

################################################################################

#simulation of random removal of clusters
#not full simulation with function
removalSim <- function(species, coord, edgeList = edgeList, id = 0:5, distance = 10) {
  commTab <- communityTable(
    speciesInputTable = species,
    landscapeTable = coord
  )
  
  tibble(
    id = id, commTabs = list(commTab)) |>
    mutate(commTabs = accumulate(commTabs, \(commTabs, lonelyPoint) { 
      #putting in the old function for the randomization
      lonely <- commTabs |> 
        slice_sample(n = 1)
      commTabs |> 
        mutate(
          dist1 = sqrt((dim1 - lonely$dim1[1])^2 + (dim2 - lonely$dim2[1])^2),
        ) |>
        filter(dist1 > distance) |>  # behåll bara punkter som är långt från båda
        select(!dist1)
      
    })) 
}

#test
removalSim(speciesinfo, patchcoordinates, id = 0:20, distance = 10) |> 
  print(n = 21)

################################################################################

#function for removal and creation of occupancyprob witch clusters
removalSimClustered <- function(species, coord, interactions = edgeList,  id = 0:5, distance = 5) {
  commTab <- communityTable(
    speciesInputTable = species,
    landscapeTable = coord
  )
  
  tibble(
    id = id, commTabs = list(commTab)) |>
    mutate(commTabs = accumulate(commTabs, \(commTabs, lonelyPoint) { 
      #putting in the old function for the randomization
      lonely <- commTabs |> 
        slice_sample(n = 1)
      
      commTabs |> 
        mutate(
          dist1 = sqrt((dim1 - lonely$dim1[1])^2 + (dim2 - lonely$dim2[1])^2),
        ) |>
        filter(dist1 > distance) |> 
        select(!dist1)
      
    })) |>
    mutate( commTabs = map(commTabs, ~ simMetacomm(.x, edgeList, sparse = TRUE, threshold = 1e-4))) |> 
    unnest(commTabs)
  
}

# tictoc::tic() # Start the clock
# clustertest_40_5 <- removalSimClustered(speciesinfo, interactions = edgeList,  patchcoordinates, id = 0:40, distance = 5)
# tictoc::toc() # Start the clock
# 
# clustertest_40_5 |>
#  write_csv("clustertest_40_5.csv")

################################################################################

#sparse removal simulation with occupancyprob
removalSimSparse <- function(species, interactions = edgeList, coord, id = 0:5) {
  commTab <- communityTable(
    speciesInputTable = species,
    landscapeTable = coord
  )
  
  
  tibble(
    id = id, commTabs = list(commTab)) |>
    mutate(commTabs = accumulate(commTabs, \(commTabs, sparsepoints ) { 
      want <- commTabs |> 
        distinct(patch) |> 
        slice_sample(prop = 0.5)
      
      inner_join(commTabs, want, by = join_by(patch))
      
      
    }))|> 
    mutate( commTabs = map(commTabs, ~ simMetacomm(.x, edgeList, sparse = TRUE, threshold = 1e-4))) |> 
    unnest(commTabs)
}



# tictoc::tic()
# testsparse_12 <- removalSimSparse(speciesinfo, edgeList, patchcoordinates, id = 0:12) |>
# print(n = 13)
# tictoc::toc()
# 
# 
# testsparse_12 |>
# write_csv("testsparse_12.csv")

################################################################################

#assign datasets

sparse_12 <- read_csv("testsparse_12.csv")

cluster_40_5 <- read_csv("clustertest_40_5.csv")


###############################Diagrams#########################################

specieslabels <- labeller(species = c("Butterfly" = "M cinxia",
                                      "Parasitoidwasp" = "C melitaearum",
                                      "PL" = "P lanceolata", 
                                      "VS" = "V spicata",
                                      "PLmildew" = "P plantaginis"))
# Occupancyprob - trend

sparse_12 |> 
  group_by(id, species) |> 
  summarise(occupancy = mean(occupancy),
            lambda = mean(lambda)) |> 
  ggplot(aes(x = id, y = occupancy)) +
  geom_point(size = 2) +
  facet_grid(.~species, labeller = specieslabels)+ 
  geom_smooth(se = FALSE, color = "grey") +
  labs( x = "Number of removal simulations",
        y = "Probability of occupancy ",
        title = "Sparse method") +
  scale_x_continuous(breaks = seq(0, 12, by = 3)) +
  theme_bw(base_size = 18)

cluster_40_5 |> 
  group_by(id, species) |> 
  summarise(occupancy = mean(occupancy),
            lambda = mean(lambda)) |> 
  ggplot(aes(x = id, y = occupancy)) +
  geom_point(size = 2) +
  facet_grid(.~species, labeller = specieslabels) + 
  geom_smooth(se = FALSE, color = "grey") +
  labs( x = "Number of removal simulations",
        y = "Probability of occupancy ",
        title = "Clustered method") +
  theme_bw(base_size = 18)

# Lambda - trend

sparse_12 |> 
  group_by(id, species) |> 
  summarise(occupancy = mean(occupancy),
            lambda = mean(lambda)) |> 
  ggplot(aes(x = id, y = lambda)) +
  geom_point(size = 2) +
  facet_grid(.~species, labeller = specieslabels) +
  scale_y_log10() +
  geom_hline(yintercept = 1, linetype = 2) + 
  geom_smooth(se = FALSE, color = "grey") +
  labs(x = "Number of removal simulations",
       y = "log(lambda) base 10",
       title = "Sparse method") +
  scale_x_continuous(breaks = seq(0, 12, by = 3)) +
  theme_bw(base_size = 18)


cluster_40_5 |> 
  group_by(id, species) |> 
  summarise(occupancy = mean(occupancy),
            lambda = mean(lambda)) |> 
  ggplot(aes(x = id, y = lambda)) +
  geom_point(size = 2) +
  facet_grid(.~species, labeller = specieslabels) +
  scale_y_log10() +
  geom_hline(yintercept = 1, linetype = 2) + 
  geom_smooth(se = FALSE, color = "grey") +
  labs(x = "Number of removal simulations",
       y = "log(lambda) base 10",
       title = "Clustered method") +
  theme_bw(base_size = 18)


# point diagrams

sparse_12 |> 
  filter(id %in% c(0, 3, 6, 9, 12)) |> 
  ggplot(aes(x = dim1, y = dim2, alpha = occupancy)) +
  geom_point(color = "darkgreen", fill = "forestgreen", shape = 21) +
  facet_grid(species~id, labeller = specieslabels) +
  labs(title = "Sparse method") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.spacing = unit(0.2, "lines")) +
  labs( x = element_blank(),
        y = element_blank())


cluster_40_5 |> 
  filter(id %in% c(0, 10, 20, 30, 40)) |> 
  ggplot(aes(x = dim1, y = dim2, alpha = occupancy)) +
  geom_point(color = "darkgreen", fill = "forestgreen", shape = 21) +
  facet_grid(species~id, labeller = specieslabels) +
  labs(title = "Clustered method") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.spacing = unit(0.2, "lines")) +
  labs( x = element_blank(),
        y = element_blank())




## just points

completepatchdata |>
  ggplot(aes(x = X, y = Y)) +
  geom_point(color = "darkgreen", fill = "forestgreen", shape = 21, size = 3) +
  labs(title = " Patches on Åland") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.spacing = unit(0.2, "lines")) +
  labs( x = element_blank(),
        y = element_blank())


# One species
cluster_40_5 |> 
  filter(id %in% c(0, 10, 20, 30, 40)) |> 
  ggplot(aes(x = dim1, y = dim2, alpha = occupancy)) +
  geom_point(color = "darkgreen", fill = "forestgreen", shape = 21) +
  facet_grid(species~id, labeller = specieslabels) +
  labs(title = "Clustered method") +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.spacing = unit(0.2, "lines")) +
  labs( x = element_blank(),
        y = element_blank())
  
  