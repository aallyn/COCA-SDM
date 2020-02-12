### Season, differences by region
df<- data.frame(rbind(res.plot.all[[2]], res.plot.all[[3]]))
df<- df %>%
  filter(., as.character(Climate_Only) == "Mean")
df$Region_Only<- factor(df$Region_Only, levels = c("South", "GOM"), labels = c("Southern New England/Mid Atlantic Bight", "Gulf of Maine"))
df.null<- cbind(expand.grid(COMNAME = levels(df$COMNAME), SEASON = unique(df$SEASON), Climate_Only = unique(df$Climate_Only), Region_Only = levels(df$Region_Only), Change = NA))
df.null<- df.null %>%
  left_join(., func.groups, by = "COMNAME")
df<- rbind(df[,], df.null)
df$duplicated<- paste(df$COMNAME, df$SEASON, df$Region_Only, df$Functional.Group)
df<- df[!duplicated(df$duplicated),] %>%
  arrange(., COMNAME, SEASON) 
spp.keep<- df %>%
  group_by(COMNAME) %>%
  summarize_at(vars(Change), n_distinct, na.rm = T) %>%
  filter(., Change == 4) %>%
  dplyr::select(., COMNAME)
df<- df %>%
  filter(., as.character(COMNAME) %in% as.character(spp.keep$COMNAME))

spp.keep2<- c("ATLANTIC COD", "RED HAKE", "WHITE HAKE", "SILVER HAKE", "ATLANTIC HERRING", "NORTHERN SAND LANCE", "ALEWIFE", "ATLANTIC MENHADEN", "AMERICAN LOBSTER")
cod.dat<- df %>%
  filter(., as.character(COMNAME) %in% as.character(spp.keep2) & Region_Only == "Gulf of Maine")

spp.keep3<- c("AMERICAN LOBSTER", "ATLANTIC COD", "SMOOTH DOGFISH", "SPINY DOGFISH", "LITTLE SKATE", "RED HAKE", "HADDOCK", "SCUP", "TAUTOG", "BLACK SEA BASS")
lob.dat<- df %>% 
  filter(., as.character(COMNAME) %in% as.character(spp.keep3) & Region_Only == "Gulf of Maine")

# Plot them
dat.use<- cod.dat
dodge <- position_dodge(width = 1)

dat.use$COMNAME<- str_to_title(dat.use$COMNAME)
dat.use$COMNAME<- factor(dat.use$COMNAME, levels = rev(unique(dat.use$COMNAME)))

dat.use.df<- dat.use %>%
  drop_na(Change)

dat.use.df$COMNAME.Plot<- factor(dat.use.df$COMNAME, levels = rev(unique(dat.use.df$COMNAME)), labels = rev(unique(to_sentence_case(as.character(dat.use.df$COMNAME)))))
dat.use.df$COMNAME.Plot<- factor(dat.use.df$COMNAME.Plot, levels = rev(c("Atlantic cod", "Alewife", "Atlantic herring", "Atlantic menhaden", "Northern sand lance", "Red hake", "Silver hake", "White hake", "American lobster")))
#dat.use.df$COMNAME.Plot<- factor(dat.use.df$COMNAME.Plot, levels = rev(c("American lobster", "Atlantic cod", "Haddock", "Red hake", "Black sea bass", "Scup", "Tautog", "Little skate", "Smooth dogfish", "Spiny dogfish"))) 

plot.means<- ggplot(data = dat.use.df, aes(x = COMNAME.Plot, y = Change, fill = SEASON)) + 
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.6)) +
  scale_fill_manual(name = "Season", values  = c('#d95f02', '#1b9e77')) +
  ylab("Percent change in relative biomass") + 
  xlab("Species") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(text = element_text(size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black")) +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) +
  ggtitle("2055 projected distribution and abundance changes for\nAtlantic cod and some of its key prey species")
  #ggtitle("2055 projected distribution and abundance changes for\nAmerican lobster and some of its key predator species")

ggplot2::ggsave(filename = "~/Desktop/AtlanticCodandPrey.jpg", plot = plot.means, width = 11, height = 8, units = "in")



library(sf)
temp<- st_read("~/Box/RES Data/Shapefiles/NELME_Regions/NELME_sf.shp")
temp2<- raster("~/Box/RES Data/Shapefiles/NEShelf_Etopo1_bathy.tiff")

temp3<- mask(temp2, temp)
plot(temp3)
mean(temp3, na.rm = TRUE)
