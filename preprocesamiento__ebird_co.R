# ------------------------------------------------------------------------------
# Introducción y Configuración

# Este esquema de preprocesamiento se basa en el GitHub de 
# [Best Practices for Using eBird Data](https://ebird.github.io/ebird-best-practices/intro.html).

# A pesar de los puntos fuertes de los datos de eBird, las observaciones de 
# especies recogidas a través de proyectos de ciencia ciudadana presentan una 
# serie de retos que no se encuentran en los datos científicos convencionales. 
# A continuación se exponen algunos de los principales retos asociados a estos 
# datos; retos que se abordarán a lo largo de este libro:

# Sesgo Taxonómico: los participantes suelen tener preferencias por 
# determinadas especies, lo que puede llevar a registrar preferentemente unas 
# especies sobre otras (Greenwood 2007; Tulloch y Szabo 2012). Restringir los 
# análisis a listas de control completas mitiga en gran medida este problema.

# Sesgo Espacial: la mayoría de los participantes en encuestas de ciencia 
# ciudadana toman muestras cerca de sus casas (Luck et al. 2004), en zonas de 
# fácil acceso como los bordes de las carreteras (Kadmon, Farber y Danin 2004) o
# en zonas y hábitats de alta biodiversidad conocida (Prendergast et al. 1993). 
# Un método sencillo para reducir el sesgo espacial que describimos consiste en 
# crear una cuadrícula de igual superficie sobre la región de interés y muestrear 
# un número determinado de listas de control dentro de cada celda de la cuadrícula.

# Sesgo Temporal: los participantes toman muestras preferentemente cuando están 
# disponibles, como los fines de semana (Courter et al. 2013), y en épocas del 
#año en las que esperan observar más aves, especialmente durante la migración 
# primaveral (Sullivan et al. 2014). Para abordar el sesgo del fin de semana, 
# recomendamos utilizar una escala temporal de una semana o varias semanas para 
# la mayoría de los análisis.

# Precisión Espacial: la ubicación espacial de una lista de control de eBird 
# se da como un único punto de latitud-longitud; sin embargo, esto puede no ser 
# preciso por dos razones principales. En primer lugar, para las listas de control 
# viajeras, esta ubicación representa sólo un punto del viaje. En segundo lugar, 
# las listas de control de eBird a menudo se asignan a un hotspot (una ubicación 
# común para todos los observadores de aves que visitan un lugar popular de 
# observación de aves) en lugar de a su verdadera ubicación. Por estas razones, 
# no es apropiado alinear las ubicaciones de eBird con variables de hábitat muy 
# precisas, y recomendamos resumir las variables dentro de un vecindario alrededor 
# de la ubicación de la lista de control.

# Desbalance de Clases: las especies de aves raras o difíciles de detectar pueden 
# tener datos con un elevado desequilibrio de clases, con muchas más listas de 
# control con no detecciones que con detecciones. Para estas especies, un modelo
# de distribución que prediga que la especie está ausente en todas partes tendrá
# una gran precisión, pero ningún valor ecológico. Seguiremos los métodos para 
# abordar el desequilibrio de clases propuestos por Robinson et al. (2018).

# Variación de la Detectabilidad: la detectabilidad describe la probabilidad de 
# que una especie presente en una zona sea detectada e identificada. La 
# detectabilidad varía según la estación, el hábitat y la especie 
# (Johnston et al. 2014, 2018). Además, los datos de eBird se recopilan con una 
# alta variación en el esfuerzo, la hora del día, el número de observadores y 
# las condiciones externas, como el clima, todo lo cual puede afectar la 
# detectabilidad de las especies (Ellis y Taylor 2018; Oliveira et al. 2018). 
# Por lo tanto, es particularmente importante tener en cuenta la detectabilidad 
# al comparar entre estaciones, hábitats o especies.Dado que eBird utiliza un 
# protocolo semiestructurado, que recoge variables asociadas con la variación en
# la detectabilidad, podremos tener en cuenta una mayor proporción de esta 
# variación en nuestros análisis.

# ------------------------------------------------------------------------------
# PREPROCESSING
# ------------------------------------------------------------------------------
# Required libraries
library(auk)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(hrbrthemes)
library(lubridate)
library(readr)
library(rnaturalearth)
library(sf)

# ------------------------------------------------------------------------------
# GIS data: country of interest
# File storage
gpkg_file <- "L:/Usuarios/pc.laboratorio.dz/Desktop/project_mp/ebird-best-practices-data/data/gis-data.gpkg"
dir.create(dirname(gpkg_file),showWarnings = FALSE, recursive = TRUE)

# Countries borderlines
ne_countries <- ne_download(scale = 50, category = "cultural",
                            type = "admin_0_countries_lakes",
                            returnclass = "sf") |>
  select(country = ADMIN, country_code = ISO_A2)

# Country selection
ne_colombia <- ne_countries |> filter(country_code == "CO")

# Se guarda el shape como geopackage
unlink(gpkg_file)
write_sf(ne_colombia, gpkg_file, "CO")

# ------------------------------------------------------------------------------
# Loading data: Checklists
f_sed_co <- "L:/Usuarios/pc.laboratorio.dz/Desktop/project_mp/ebd_CO_smp_relJan-2024_sampling.txt"
checklists_co <- read_sampling(f_sed_co)
n_entries_sed <- nrow(checklists_co)
glimpse(checklists_co)


# Distribution of distance traveling for traveling protocol checklists
checklists_traveling_co <- filter(checklists_co, protocol_type == "Traveling")

histogram <- ggplot(checklists_traveling_co) +
  aes(x = effort_distance_km) +
  geom_histogram(binwidth = 1, 
                 aes(y = after_stat(count / sum(count))),
                 fill="#0f0cc7", color="#e9ecef", alpha=0.9) +
  scale_y_continuous(limits = c(0, NA), labels = scales::label_percent()) +
  theme_ipsum() +
  labs(x = "Distancia recorrida [km]",
       y = "% de listas de eBird",
       title = "Distribución de la distancia recorrida en las listas de eBird en CO")+
  theme(plot.title = element_text(size=15))

print(histogram)

# Getting the percentages from each column
percentages <- ggplot_build(histogram)$data[[1]]$y*100

# Cummulative sum from percentages
cum_sum_per <- cumsum(percentages)

# Cummulative sum plot
df_lollipop <- data.frame(
  x = seq_along(cum_sum_per),
  y = cum_sum_per
)

# Threshold value
thresh <- 95

# Column to assign the categories regarding the threshold value
df_lollipop$group <- ifelse(df_lollipop$y < thresh, paste("Below Threshold (", thresh, "%)", sep = ""), 
                            paste("Above Threshold (", thresh, "%)", sep = ""))

# Dynamic legends
names_legend <- c(paste("Below Threshold (", thresh, "%)", sep = ""), 
                  paste("Above Threshold (", thresh, "%)", sep = ""))

# Lollipop plot
ggplot(df_lollipop, aes(x = x, y = y, color = group)) +
  geom_segment(aes(xend = x, yend = 0)) +
  geom_point(size = 2.5, fill = alpha("orange", 0.3), alpha = 0.7, shape = 21, stroke = 2) +
  geom_point(data = df_lollipop, aes(color = group), size = 4) +
  scale_color_manual(values = c("#0f0cc7", "red"), labels = names_legend) +
  theme_ipsum() +
  labs(x = "Distancia recorida [km]", 
       y = "Suma acumulada del % de listas de eBird", 
       title = "Suma acumulada del % de listas de eBird vs Distancia recorrida",
       color = "Grupo") + 
  theme(plot.title = element_text(size = 15))

# ANALYSIS: Around 95% of the Traveled checklists are about 9 km

# ------------------------------------------------------------------------------
# Loading data: Observations
f_ebd_co <- "L:/Usuarios/pc.laboratorio.dz/Desktop/project_mp/ebd_CO_smp_relJan-2024.txt"
observations_co <- read_ebd(f_ebd_co)
n_entries_ebd <- nrow(observations_co)
glimpse(observations_co)

# ------------------------------------------------------------------------------
# Loading endangered status: Colombia
lista_estado_co <-  read.csv('L:/Usuarios/pc.laboratorio.dz/Desktop/project_mp/ status_endemic_birds_co.csv',
                             header=TRUE)

# ------------------------------------------------------------------------------
# Species Filtering: 

observations_co <- observations_co |> 
  filter(scientific_name %in% c("Anisognathus melanogenys", "Anthocephala berlepschi", "Anthocephala floriceps", 
                                "Arremon basilicus", "Atlapetes blancae", "Atlapetes flaviceps", 
                                "Atlapetes fuscoolivaceus", "Atlapetes melanocephalus", "Bangsia aureocincta", 
                                "Bangsia melanochlamys", "Bolborhynchus ferrugineifrons", "Bucco noanamae",
                                "Campylopterus phainopeplus", "Capito hypoleucus", "Cercomacroides parkeri",
                                "Chaetocercus astreans", "Chlorochrysa nitidissima", "Chlorostilbon olivaresi",
                                "Chrysuronia lilliae", "Cistothorus apolinari", "Clibanornis rufipectus",
                                "Coeligena orina", "Coeligena phalerata", "Coeligena prunellei", "Cranioleuca hellmayri",
                                "Crax alberti", "Dacnis hartlaubi", "Diglossa gloriosissima", "Drymophila caudata",
                                "Drymophila hellmayri", "Eriocnemis isabellae", "Eriocnemis mirabilis", "Euphonia concinna",
                                "Grallaria bangsi", "Grallaria kaestneri", "Grallaria milleri", "Grallaria urraoensis",
                                "Habia cristata", "Habia gutturalis", "Hapalopsittaca fuertesi", "Henicorhina anachoreta",
                                "Henicorhina negreti", "Hypopyrrhus pyrohypogaster", "Leptotila conoveri", "Lipaugus weberi",
                                "Macroagelaius subalaris", "Megascops gilesi", "Melanerpes pulcher", "Myiarchus apicalis",
                                "Myioborus flavivertex", "Myiotheretes pernix", "Myiothlypis basilica", "Myiothlypis conspicillata",
                                "Odontophorus hyperythrus", "Odontophorus strophium", "Ortalis columbiana", "Ortalis garrula",
                                "Oxypogon cyanolaemus", "Oxypogon guerinii", "Oxypogon stuebelii", "Penelope perspicax",
                                "Phylloscartes lanyoni", "Picumnus granadensis", "Psarocolius cassini", "Pyrrhura calliptera",
                                "Pyrrhura viridicata", "Rallus semiplumbeus", "Ramphomicron dorsale", "Saucerottia castaneiventris",
                                "Saucerottia cyanifrons", "Scytalopus alvarezlopezi", "Scytalopus canus", "Scytalopus latebricola",
                                "Scytalopus rodriguezi", "Scytalopus sanctaemartae", "Scytalopus stilesi", "Synallaxis fuscorufa",
                                "Synallaxis subpudica", "Thryophilus nicefori", "Thryophilus sernai", "Troglodytes monticola",
                                "Vireo approximans", "Vireo caribaeus"))

# ------------------------------------------------------------------------------

# Check default values from SED and EBD: Shared Checklists and Taxonomic Roll-up
formals(read_sampling)
formals(read_ebd)

# ------------------------------------------------------------------------------
# Data time filtering: checklists and observations
checklists_co <- checklists_co |> 
  filter(all_species_reported,
         between(year(observation_date), 2003, 2023))
n_entries_sed_fil <- nrow(checklists_co)

observations_co <- observations_co |> 
  filter(all_species_reported,
         between(year(observation_date), 2003, 2023))
n_entries_ebd_fil <- nrow(observations_co)

# Quantification of discarded data
dis_entries_sed_per <- (100-(n_entries_sed_fil/n_entries_sed)*100)
dis_entries_ebd_per <- (100-(n_entries_ebd_fil/n_entries_ebd)*100)
print(paste("SED: % of discarded entries after filtering =", dis_entries_sed_per))
print(paste("EBD: % of discarded entries after filtering =", dis_entries_ebd_per))

# ------------------------------------------------------------------------------
# Data spatial filtering: checklists and observations

# Convert checklist locations to points geometries
checklists_co_sf <- checklists_co |> 
  select(checklist_id, latitude, longitude) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = "WGS84")

# Boundary of study region, buffered by 1 km
study_region_buffered <- read_sf("L:/Usuarios/pc.laboratorio.dz/Desktop/project_mp/ebird-best-practices-data/data/gis-data.gpkg",
                                 layer = "CO") |>
  st_transform(crs = st_crs(checklists_co_sf))

# Spatially subset the checklists to those in the study region
in_region <- checklists_co_sf[study_region_buffered, ]
rm(checklists_co_sf)

# Join to checklists and observations to remove checklists outside region
checklists_co <- semi_join(checklists_co, in_region, by = "checklist_id")
observations_co <- semi_join(observations_co, in_region, by = "checklist_id")

# Quantification of discarded data
dis_entries_sed_per_sf <- (100-(nrow(checklists_co)/n_entries_sed)*100)
dis_entries_ebd_per_sf <- (100-(nrow(observations_co)/n_entries_ebd)*100)
print(paste("SED: % of discarded entries after spatial filtering =", dis_entries_sed_per_sf-dis_entries_sed_per))
print(paste("EBD: % of discarded entries after spatial filtering =", dis_entries_ebd_per_sf-dis_entries_ebd_per))

# ------------------------------------------------------------------------------
# Zero filling
zf <- auk_zerofill(observations_co, checklists_co, collapse = TRUE)
rm(observations_co,checklists_co)

# Joining to endangered status in Colombia
zf <- left_join(zf, lista_estado_co, by = "scientific_name")


# Function to convert time observation to hours since midnight
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# Clean up variables
zf <- zf |> 
  mutate(
    # Convert count to integer and X to NA
    # ignore the warning "NAs introduced by coercion"
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for stationary counts
    effort_distance_km = if_else(protocol_type == "Stationary", 
                                 0, effort_distance_km),
    # convert duration to hours
    effort_hours = duration_minutes / 60,
    # speed km/h
    effort_speed_kmph = effort_distance_km / effort_hours,
    # convert time to decimal hours since midnight
    hours_of_day = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# Quantification of discarded data once is zero-filled
n_entries_zf_data <- nrow(zf)

# ------------------------------------------------------------------------------
# Accounting for variation in effort

# Additional filtering
zf <- zf |> 
  filter(protocol_type %in% c("Stationary", "Traveling"),
         effort_hours <= 6,
         effort_distance_km <= 10,
         effort_speed_kmph <= 100,
         number_observers <= 10)

# Quantification of discarded data once is zero-filled
dis_entries_zf_per <- (100-(nrow(zf)/n_entries_zf_data)*100)
print(paste("% of discarded entries after zero-filling and filtering =", dis_entries_zf_per))

# ------------------------------------------------------------------------------
# Endemic Species with higher counts

counts_per_specie <- aggregate(observation_count ~ scientific_name, data=zf, FUN=sum)
lista_estado_co <- left_join(lista_estado_co, counts_per_specie, by = "scientific_name")
ordered_counts_per_specie <- counts_per_specie[order(-counts_per_specie$observation_count), ]
print(ordered_counts_per_specie)

# ------------------------------------------------------------------------------
# Filtrar las especies de acuerdo a su estado de amenaza

lc_data_co <- zf |> filter(status_co %in% "Preocupacion menor")
nt_data_co <- zf |> filter(status_co %in% "Casi amenazado")
vu_data_co <- zf |> filter(status_co %in% "Vulnerable")
en_data_co <- zf |> filter(status_co %in% "En peligro")
cr_data_co <- zf |> filter(status_co %in% "En peligro critico")

# ------------------------------------------------------------------------------
# Data storage
path <- "L:/Usuarios/pc.laboratorio.dz/Desktop/project_mp/"
write.csv(lc_data_co, paste0(path,"2003_2023_aves_endemicas_co_lc.csv"), row.names=TRUE)
write.csv(nt_data_co, paste0(path,"2003_2023_aves_endemicas_co_nt.csv"), row.names=TRUE)
write.csv(vu_data_co, paste0(path,"2003_2023_aves_endemicas_co_vu.csv"), row.names=TRUE)
write.csv(en_data_co, paste0(path,"2003_2023_aves_endemicas_co_en.csv"), row.names=TRUE)
write.csv(cr_data_co, paste0(path,"2003_2023_aves_endemicas_co_cr.csv"), row.names=TRUE)
write.csv(lista_estado_co, paste0(path,"observaciones_endemicas_ordenadas.csv"), row.names=TRUE)