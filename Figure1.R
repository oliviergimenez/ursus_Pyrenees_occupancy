# packages
library(rgdal)
 
# read in shapefiles
sousmassif.rg <- readOGR(".", "sous_massif_08_12_MaJ2016")
print(proj4string(sousmassif.rg))
deptmt <- readOGR(".", "ensemble_departement")
print(proj4string(deptmt))
x = proj4string(sousmassif.rg) 
# on reprojette pour avoir les lat/long
sousmassif.rg <- spTransform(sousmassif.rg, CRS ("+init=epsg:4326") )
deptmt <- spTransform(deptmt, CRS ("+init=epsg:4326") )

ppi <- 300
tiff("figure1.tiff", width=6*ppi, height=6*ppi, res=ppi)

# subsections and counties
plot(sousmassif.rg, axes=T, border="gray")#,xlab='longitude',ylab='latitude')
plot(deptmt, border="black",add=T)
box()

# scale
lines(c(-0.9924132,-0.9924132+0.622964),c(41.6,41.6))
lines(c(-0.9924132,-0.9924132),c(41.6,41.6+.05))
lines(rep(-0.9924132+0.622964,2),c(41.6,41.6+.05),type='l')
text(-0.9924132,41.6+.1,'0',cex=1)
text(-0.9924132+0.622964,41.6+.1,'50',cex=1)
text(-0.9924132+0.622964/2,41.6-0.07,'km',cex=1)

# add ref for counties
# use x <- drawExtent() to locate best spot
text(-0.6797100,43.33997,'64',cex=1.1)
text(0.1919815,43.26513,'65',cex=1.1)
text(1.4411108,43.53238,'31',cex=1.1)
text(2.3428967,43.21026,'11',cex=1.1)
text(1.7371359,43.08394,'09',cex=1.1) 
text(2.8331056,42.67632,'66',cex=1.1)

text(0.9830889,42.03091, '64: Pyrénées Atlantiques',adj = c(0,0))
text(0.9851593,41.94119, '65: Hautes Pyrénées',adj = c(0,0))
text(0.9872234,41.85148, '31: Haute Garonne',adj = c(0,0))
text(0.9892812,41.76178, '09: Ariège',adj = c(0,0))
text(0.9913327,41.67208, '11: Aude',adj = c(0,0))
text(0.9933780,41.58240, '66: Pyrénées Orientales',adj = c(0,0))

# add countries
text(0.2187799,42.46478,'SPAIN',adj=c(0,0))
text(2.1683432,43.82339,'FRANCE',adj=c(0,0))

dev.off()
