.onLoad <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
}
.onAttach <- function(lib, pkg){
 packageStartupMessage("This is libamtrack 0.6.0 'Black Armadillo' (2015-05-21).\nType '?libamtrack' for help.\n")
}
.onUnload <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}
