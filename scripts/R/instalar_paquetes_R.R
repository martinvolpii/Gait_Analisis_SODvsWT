# Paquetes de R para el análisis estadístico final. Se ejecuta una sola vez.
pkgs <- c("lme4", "lmerTest", "emmeans", "ggplot2", "pbkrtest")
falta <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
repos <- getOption("repos")
if (is.null(repos) || identical(unname(repos[["CRAN"]]), "@CRAN@")) repos <- c(CRAN = "https://cloud.r-project.org")
if (length(falta)) install.packages(falta, repos = repos) else message("Ya estan todos instalados.")
