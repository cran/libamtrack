.onLoad <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 packageStartupMessage("This is libamtrack 0.5.2 'Red Wombat' (2011-12-22).\nType '?libamtrack' for help.\n")
}
.Last.lib <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}
