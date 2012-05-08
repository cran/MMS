## .onLoad <- function(lib, pkg)
##   {
##     require(methods)
##     require(emulator)
##   }

.onAttach <- function(libname, pkgname){ packageStartupMessage("Loaded MMS ",as.character(packageDescription("MMS")[["Version"]]))}