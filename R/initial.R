.onLoad <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
}
.onAttach <- function(lib, pkg){
 packageStartupMessage("This is libamtrack 0.5.4 'Yellow Wombat' (2013-01-22).\nType '?libamtrack' for help.\n")
}
.onUnload <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}
