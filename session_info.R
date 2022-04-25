fileConn<-file("session_info.txt")
pkgs <- c()
for (package_name in sort(loadedNamespaces())) {
  pkgs <- c(pkgs,paste(package_name, packageVersion(package_name)))
}
writeLines(pkgs, fileConn)
close(fileConn)


