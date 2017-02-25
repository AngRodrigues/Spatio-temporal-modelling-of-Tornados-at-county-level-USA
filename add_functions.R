shp2raster <- function(
  shpname="",    # single file name like "coolstuff.shp". Ignored if shp is given.
  shp=NULL,      # Shapefile Object Spatial*DataFrame. If NULL, it reads shpname with rgdal::readOGR.
  ncells=99,     # Approximate number of cells in either direction to determine cellsize.
  cellsize=NA,   # Cell size in coordinate units (usually degrees or m). Computed from ncells if NA.
  ncellwarn=1000,# Warn if there will be more cells than this. To prevent e.g. accidental degrees instead of km.
  column="",     # Name of column to use for z dimension in raster. Empty string for interactive selection.
  ascname=NA,    # Output file name. If NA, inferred from shpname or shp.
  verbose=FALSE, # Report readOGR progress?
  ...)           # More arguments passed to raster::rasterize, like overwrite=TRUE
{
  # if shp is missing/default, read shpname:
  if(is.null(shp)) 
  {
    shp <- rgdal::readOGR(dsn=shpname, 
                          layer=basename(tools::file_path_sans_ext(shpname)),
                          verbose=verbose)
    if(is.na(ascname)) ascname <- sub(".shp", ".asc", shpname)
  } else
    if(is.na(ascname)) ascname <- paste0(deparse(substitute(shp)),".asc")
    # target raster extend and resolution:
    e <- extent(shp) 
    if(is.na(cellsize)) cellsize <- mean(c((e@xmax-e@xmin), (e@ymax-e@ymin))/ncells)
    nx <- (e@xmax-e@xmin)/cellsize # this seems revertive from the previous line, but
    ny <- (e@ymax-e@ymin)/cellsize # is needed because ncells differ in both directions
    cont <- TRUE # continue by default
    if(max(nx,ny)>ncellwarn) cont <- readline(paste0("Raster will be large: nx=",
                                                     round(nx,1), ", ny=",round(ny,1)," (with cellsize=", round(cellsize,4),", xmin=",
                                                     round(e@xmin,2), ", xmax=",round(e@xmax,2),"). Continue? y/n: "))
    cont <- tolower(cont) %in% c("y", "yes", "t", "true", "")
    if(!cont) return(list(nx=nx, ny=ny, cellsize=cellsize, extend_shp=e))
    r <- raster(ncol=nx, nrow=ny)
    extent(r) <- extent(shp)
    resdif <- abs((yres(r) - xres(r)) / yres(r) )
    if(resdif > 0.01) stop("Horizontal (",round(xres(r),3),") and vertical (", round(yres(r),3),
                           ") resolutions are too different (diff=",round(resdif,3), ", but must be <0.010).\n",
                           "  Use a smaller cell size to achieve this (currently ",round(cellsize,1),").")
    # column selection
    n <- names(shp)
    if(!column %in% n) message("Column '",column, "' is not in Shapefile. Select one of\n", 
                               paste(strwrap(toString(n)), collapse="\n"))
    while(!column %in% n) column <- readline(paste0("Nonexistent column '",column, 
                                                    "'. Type desired name, then hit ENTER: "))
    # actually convert and write to file:
    ras <- raster::rasterize(shp, r, column, filename=ascname, proj=shp@proj4string, ...)
    # return output
    ras
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

brier.score <- function(x, m){
  with(m, {mean(x^2) - 2 * mean(x * mean) + mean(mean^2 + sd^2)})
}
