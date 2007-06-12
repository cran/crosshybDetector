#Copyright 2004, W. Wolski, all rights reserved.
.First.lib <- function(lib, pkg) library.dynam("crosshybDetector",pkg,lib)
.Last.lib <- function(libpath) library.dynam.unload("crosshybDetector", libpath)

