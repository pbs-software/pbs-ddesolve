#openVig--------------------------------2010-05-04
# Open package vignette 'pkgIntro.pdf' if it exists.
#-----------------------------------------------RH
openVig = function(pkg="ddesolve"){
	pkgnam = as.character(substitute(pkg))
	vignam = system.file(paste("doc/",pkgnam,"Intro.pdf",sep=""), package=pkgnam )
	if (vignam!="" && file.exists(vignam)) {
		if (Sys.info()["sysname"]=="Windows") 
			shell.exec(vignam)
		else
			try(shell(vignam),silent=TRUE)
		cat("Vignette located at:\n",vignam,"\n")
	}
	else stop("There is no 'Intro.pdf' file for this package")
}
