# Verificación: LRT del modelo mixto en R (lme4, bobyqa) para las 18 variables, comparado con Python.
# Carpeta con los resultados de Python (Analisis_marcha_v2). Se busca en este orden:
#   1) argumento:  Rscript analisis_rodilla.R "C:/ruta/Analisis_marcha_v2"
#   2) variable de entorno MARCHA_BASE
#   3) <carpeta de trabajo>/resultados/Analisis_marcha_v2   (en RStudio: setwd() a la raíz del repositorio)
# En Linux/macOS con localización "C" R no lee UTF-8; se fuerza una localización UTF-8 si hace falta.
if (!l10n_info()[["UTF-8"]] && .Platform$OS.type != "windows") {
  for (loc in c("C.UTF-8", "en_US.UTF-8", "es_ES.UTF-8")) {
    if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break
  }
}
args <- commandArgs(trailingOnly = TRUE)
BASE <- if (length(args) >= 1) args[1] else Sys.getenv("MARCHA_BASE", "")
if (BASE == "") {
  cand <- c(file.path(getwd(), "resultados", "Analisis_marcha_v2"))
  BASE <- cand[dir.exists(cand)][1]
}
if (is.na(BASE) || !dir.exists(BASE)) stop("No encuentro la carpeta Analisis_marcha_v2. Pase la ruta como argumento.")
suppressPackageStartupMessages(library(lme4))
DIAS <- c("P30","P37","P44","P51","P65","P85")
arch <- list.files(BASE, pattern = "_por_animal\\.csv$", recursive = TRUE, full.names = TRUE)
arch <- arch[!grepl("/R/", arch)]
ctrl <- lmerControl(optimizer = "bobyqa")
out <- do.call(rbind, lapply(arch, function(f) {
  w <- read.csv2(f, fileEncoding = "UTF-8-BOM", check.names = FALSE, colClasses = c(animal = "character"))
  d <- do.call(rbind, lapply(DIAS, function(k) data.frame(genotipo = factor(w$genotipo, c("WT","SOD1")),
                                                          animal = factor(w$animal), dia = factor(k, DIAS), y = as.numeric(w[[k]]))))
  d <- d[!is.na(d$y), ]
  fit <- function(fo) suppressMessages(lmer(fo, d, REML = FALSE, control = ctrl))
  full <- fit(y ~ genotipo * dia + (1 | animal)); add <- fit(y ~ genotipo + dia + (1 | animal)); nul <- fit(y ~ dia + (1 | animal))
  a1 <- anova(add, full); a2 <- anova(nul, add)
  data.frame(codigo = sub("_.*", "", basename(f)), p_int_R = a1$`Pr(>Chisq)`[2], p_gen_R = a2$`Pr(>Chisq)`[2])
}))
py <- read.csv2(file.path(BASE, "RESUMEN_todas_las_variables.csv"), fileEncoding = "UTF-8-BOM", check.names = FALSE)
out <- merge(out, py[, c("codigo", "variable", "p_interaccion_GxD", "p_genotipo_principal")], by = "codigo")
names(out)[5:6] <- c("p_int_Python", "p_gen_Python")
out$q_int_R <- p.adjust(out$p_int_R, "BH"); out$q_gen_R <- p.adjust(out$p_gen_R, "BH")
out$difiere <- abs(log10(out$p_int_R / out$p_int_Python)) > 0.1 | abs(log10(out$p_gen_R / out$p_gen_Python)) > 0.1
out <- out[order(out$p_int_R), c("codigo","variable","p_int_Python","p_int_R","q_int_R","p_gen_Python","p_gen_R","q_gen_R","difiere")]
print(format(out, digits = 3), row.names = FALSE)
