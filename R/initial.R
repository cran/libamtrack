.onLoad <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 packageStartupMessage("This is libamtrack 0.5.3 'Green Wombat' (2012-04-27).\nType '?libamtrack' for help.\n")
}
.Last.lib <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}
