# bluster and org.Hs.eg.db have no URL in their DESCRIPTION, so pkgdown
# cannot link library() to a package page. Send those two to Bioconductor.
setHook(packageEvent("downlit", "onLoad"), function(...) {
  ns <- asNamespace("downlit")
  orig <- ns$href_package
  extra <- c(
    bluster = "https://bioconductor.org/packages/bluster",
    "org.Hs.eg.db" = "https://bioconductor.org/packages/org.Hs.eg.db"
  )
  unlockBinding("href_package", ns)
  assign("href_package", function(package) {
    if (package %in% names(extra)) extra[[package]] else orig(package)
  }, envir = ns)
}, action = "append")
