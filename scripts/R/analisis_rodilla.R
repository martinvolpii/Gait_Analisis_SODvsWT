# =============================================================================
# Ángulo de rodilla longitudinal, WT vs SOD1 (P30–P85). Análisis en R.
# Réplica del análisis de Python (estadistica_v2.py) con lme4 / lmerTest / emmeans.
#
# Unidad experimental: el animal. Dato de entrada: media por animal y día de todos
# los ciclos válidos (archivos *_por_animal.csv de la carpeta 04_...).
#
# Pruebas:
#   1) Modelo lineal mixto: y ~ genotipo * día + (1 | animal)
#      - Prueba de razón de verosimilitudes (LRT, ajuste ML) para la interacción
#        genotipo × día, el efecto principal de genotipo y el efecto global.
#        Es la misma prueba que en Python (statsmodels), sirve para comparar.
#      - Prueba F tipo III con gl de Kenward-Roger (ajuste REML). Más conservadora
#        con n pequeño; es la que conviene reportar en el manuscrito.
#   2) Contrastes SOD1 − WT por día desde el modelo mixto (emmeans, gl de
#      Kenward-Roger), con Holm sobre los 6 días.
#   3) Welch por día + Holm sobre 6 días, Mann-Whitney exacto y g de Hedges
#      (igual que en Python).
#   4) Supuestos: Shapiro-Wilk de los residuos y Brown-Forsythe (Levene con
#      mediana) entre las 12 celdas genotipo × día.
#
# Uso:  Rscript scripts/R/analisis_rodilla.R resultados/Analisis_marcha_v2
# =============================================================================

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
DIR_IN  <- file.path(BASE, "04_Goniometria_rodilla_tobillo_pie")
DIR_OUT <- file.path(BASE, "R")
dir.create(DIR_OUT, showWarnings = FALSE)

# ---- paquetes (se instalan solo si faltan) -----------------------------------
pkgs  <- c("lme4", "lmerTest", "emmeans", "ggplot2", "pbkrtest")
falta <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(falta)) install.packages(falta, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(emmeans); library(ggplot2)
})

DIAS <- c("P30", "P37", "P44", "P51", "P65", "P85")
VARS <- data.frame(
  codigo  = c("V4a", "V4c", "V4b"),
  archivo = c("V4a_rodilla_max", "V4c_rodilla_med", "V4b_rodilla_min"),
  titulo  = c("Extensión máxima", "Ángulo medio", "Flexión máxima (ángulo mínimo)"),
  stringsAsFactors = FALSE)

# ---- lectura: formato ancho (un animal por fila) a formato largo -------------
leer <- function(arch) {
  w <- read.csv2(file.path(DIR_IN, paste0(arch, "_por_animal.csv")),
                 fileEncoding = "UTF-8-BOM", check.names = FALSE,
                 colClasses = c(animal = "character"))
  l <- do.call(rbind, lapply(DIAS, function(d)
    data.frame(genotipo = w$genotipo, animal = w$animal, dia = d, y = as.numeric(w[[d]]))))
  l <- l[!is.na(l$y), ]
  l$genotipo <- factor(l$genotipo, c("WT", "SOD1"))
  l$dia      <- factor(l$dia, DIAS)
  l$animal   <- factor(l$animal)
  rownames(l) <- NULL
  l
}

hedges_g <- function(w, s) {             # (SOD1 − WT) / DE combinada, corrección de Hedges
  nw <- length(w); ns <- length(s)
  sp <- sqrt(((nw - 1) * var(w) + (ns - 1) * var(s)) / (nw + ns - 2))
  (mean(s) - mean(w)) / sp * (1 - 3 / (4 * (nw + ns) - 9))
}

por_dia <- function(d) {
  r <- do.call(rbind, lapply(DIAS, function(k) {
    w <- d$y[d$dia == k & d$genotipo == "WT"]; s <- d$y[d$dia == k & d$genotipo == "SOD1"]
    data.frame(dia = k, n_WT = length(w), n_SOD1 = length(s),
               media_WT = mean(w), de_WT = sd(w), media_SOD1 = mean(s), de_SOD1 = sd(s),
               dif_SOD1_menos_WT = mean(s) - mean(w),
               p_welch = t.test(s, w, var.equal = FALSE)$p.value,
               p_mann_whitney = wilcox.test(s, w, exact = TRUE)$p.value,
               hedges_g = hedges_g(w, s))
  }))
  r$p_welch_holm_6dias <- p.adjust(r$p_welch, "holm")
  r
}

ctrl <- lmerControl(optimizer = "bobyqa")
lrt  <- function(m_red, m_full) { a <- anova(m_red, m_full); c(chi2 = a$Chisq[2], gl = a$Df[2], p = a$`Pr(>Chisq)`[2]) }

analizar <- function(d) {
  # LRT con ML (mismo criterio que Python)
  full_ml <- lmer(y ~ genotipo * dia + (1 | animal), d, REML = FALSE, control = ctrl)
  add_ml  <- lmer(y ~ genotipo + dia + (1 | animal), d, REML = FALSE, control = ctrl)
  nul_ml  <- lmer(y ~ dia + (1 | animal),            d, REML = FALSE, control = ctrl)
  L_int <- lrt(add_ml, full_ml); L_gen <- lrt(nul_ml, add_ml); L_glo <- lrt(nul_ml, full_ml)

  # F tipo III, Kenward-Roger (REML, contrastes de suma)
  full_reml <- lmer(y ~ genotipo * dia + (1 | animal), d, REML = TRUE, control = ctrl,
                    contrasts = list(genotipo = "contr.sum", dia = "contr.sum"))
  aov3 <- tryCatch(anova(full_reml, type = 3, ddf = "Kenward-Roger"),
                   error = function(e) anova(full_reml, type = 3, ddf = "Satterthwaite"))

  # efecto principal de genotipo (modelo aditivo, REML) e ICC del animal
  add_reml <- lmer(y ~ genotipo + dia + (1 | animal), d, REML = TRUE, control = ctrl)
  b  <- fixef(add_reml)["genotipoSOD1"]
  ci <- confint(add_reml, parm = "genotipoSOD1", method = "Wald")
  vc <- as.data.frame(VarCorr(add_reml))
  icc <- vc$vcov[vc$grp == "animal"] / sum(vc$vcov)

  # contrastes SOD1 − WT por día desde el modelo mixto, Holm sobre 6 días
  emm <- emmeans(full_reml, ~ genotipo | dia, lmer.df = "kenward-roger")
  ct  <- as.data.frame(test(pairs(emm, reverse = TRUE), by = NULL, adjust = "holm"))
  ct_raw <- as.data.frame(test(pairs(emm, reverse = TRUE), by = NULL, adjust = "none"))
  ct$p_sin_ajuste <- ct_raw$p.value

  # supuestos
  res <- residuals(full_reml)
  celda <- interaction(d$genotipo, d$dia)
  bf <- anova(lm(abs(d$y - ave(d$y, celda, FUN = median)) ~ celda))

  list(LRT = rbind(interaccion_GxD = L_int, genotipo_principal = L_gen, genotipo_global = L_glo),
       aov3 = aov3, efecto = c(SOD1_menos_WT = unname(b), IC95_inf = ci[1], IC95_sup = ci[2], ICC_animal = icc),
       contrastes = ct, shapiro_p = shapiro.test(res)$p.value, brown_forsythe_p = bf$`Pr(>F)`[1],
       n_obs = nrow(d), n_animales = nlevels(droplevels(d$animal)), modelo = full_reml)
}

# ---- ejecución ---------------------------------------------------------------
sink(file.path(DIR_OUT, "salida_rodilla_R.txt"), split = TRUE)
cat("Análisis de rodilla en R —", format(Sys.time()), "\n", R.version.string, "\n")
cat("lme4", as.character(packageVersion("lme4")), "| lmerTest", as.character(packageVersion("lmerTest")),
    "| emmeans", as.character(packageVersion("emmeans")), "\n\n")

datos <- list(); resultados <- list(); tabla_dia <- list(); resumen <- list()
for (i in seq_len(nrow(VARS))) {
  v <- VARS[i, ]; d <- leer(v$archivo); datos[[v$codigo]] <- d
  A <- analizar(d); resultados[[v$codigo]] <- A
  P <- por_dia(d)
  P <- merge(P, data.frame(dia = as.character(A$contrastes$dia),
                           mm_dif = A$contrastes$estimate, mm_gl_KR = A$contrastes$df,
                           mm_p = A$contrastes$p_sin_ajuste, mm_p_holm_6dias = A$contrastes$p.value),
             by = "dia", sort = FALSE)
  P <- P[match(DIAS, P$dia), ]
  tabla_dia[[v$codigo]] <- cbind(codigo = v$codigo, P)

  cat("=============================================================\n")
  cat(v$codigo, "—", v$titulo, sprintf("(n obs = %d, animales = %d)\n", A$n_obs, A$n_animales))
  cat("\nLRT (ML):\n"); print(round(A$LRT, 5))
  cat("\nF tipo III (REML, gl Kenward-Roger):\n"); print(A$aov3)
  cat("\nEfecto principal SOD1 − WT (modelo aditivo, REML):\n"); print(round(A$efecto, 4))
  cat("\nPor día — Welch, Mann-Whitney exacto, Hedges g, Holm 6 días; contraste del modelo mixto (KR, Holm):\n")
  print(format(P, digits = 4), row.names = FALSE)
  cat(sprintf("\nSupuestos: Shapiro-Wilk residuos p = %.3f | Brown-Forsythe 12 celdas p = %.3f\n\n",
              A$shapiro_p, A$brown_forsythe_p))
  f3 <- as.data.frame(A$aov3)
  resumen[[v$codigo]] <- data.frame(
    codigo = v$codigo, variable = v$titulo,
    p_interaccion_LRT = A$LRT["interaccion_GxD", "p"], p_genotipo_LRT = A$LRT["genotipo_principal", "p"],
    F_interaccion_KR = f3["genotipo:dia", "F value"], gl_den_KR = f3["genotipo:dia", "DenDF"],
    p_interaccion_KR = f3["genotipo:dia", "Pr(>F)"], p_genotipo_KR = f3["genotipo", "Pr(>F)"],
    efecto_SOD1_menos_WT = A$efecto[1], IC95_inf = A$efecto[2], IC95_sup = A$efecto[3], ICC_animal = A$efecto[4],
    shapiro_p = A$shapiro_p, brown_forsythe_p = A$brown_forsythe_p)
}
R <- do.call(rbind, resumen); rownames(R) <- NULL
R$q_FDR_interaccion_LRT_3 <- p.adjust(R$p_interaccion_LRT, "BH")   # FDR solo entre las 3 de rodilla

# ---- comparación con Python ----------------------------------------------------
py_f <- file.path(BASE, "RESUMEN_todas_las_variables.csv")
if (file.exists(py_f)) {
  py <- read.csv2(py_f, fileEncoding = "UTF-8-BOM", check.names = FALSE)
  py <- py[match(R$codigo, py$codigo), ]
  cmp <- data.frame(codigo = R$codigo,
                    p_int_Python = py$p_interaccion_GxD, p_int_R_LRT = signif(R$p_interaccion_LRT, 3),
                    p_int_R_KR = signif(R$p_interaccion_KR, 3),
                    p_gen_Python = py$p_genotipo_principal, p_gen_R_LRT = signif(R$p_genotipo_LRT, 3),
                    efecto_Python = round(py$efecto_SOD1_vs_WT, 3), efecto_R = round(R$efecto_SOD1_menos_WT, 3))
  cat("=============================================================\nComparación Python vs R:\n")
  print(cmp, row.names = FALSE)
}
sink()

write.csv2(R, file.path(DIR_OUT, "rodilla_R_modelo_mixto.csv"), row.names = FALSE)
write.csv2(do.call(rbind, tabla_dia), file.path(DIR_OUT, "rodilla_R_por_dia.csv"), row.names = FALSE)

# ---- figura (ggplot2) -----------------------------------------------------------
COL <- c(WT = "#2a78d6", SOD1 = "#eb6834")
etiqueta <- setNames(sprintf("%s\ninteracción genotipo × día: LRT p = %.4f; F (KR) p = %.4f",
                             VARS$titulo, R$p_interaccion_LRT, R$p_interaccion_KR), VARS$codigo)
D_all <- do.call(rbind, lapply(VARS$codigo, function(k) cbind(codigo = k, datos[[k]])))
D_all$panel <- factor(etiqueta[D_all$codigo], etiqueta)
M <- aggregate(y ~ codigo + panel + genotipo + dia, D_all,
               function(x) c(m = mean(x), ee = sd(x) / sqrt(length(x))))
M <- do.call(data.frame, M); names(M)[5:6] <- c("m", "ee")
T <- do.call(rbind, tabla_dia)
T <- T[T$p_welch < 0.05, ]
lim <- aggregate(y ~ codigo, D_all, function(x) c(lo = min(x), hi = max(x))); lim <- do.call(data.frame, lim)
T$ypos <- lim$y.hi[match(T$codigo, lim$codigo)] + 0.04 * (lim$y.hi - lim$y.lo)[match(T$codigo, lim$codigo)]
T$lab  <- ifelse(T$p_welch_holm_6dias < 0.05,
                 sprintf("p = %.3f\nHolm = %.3f", T$p_welch, T$p_welch_holm_6dias), sprintf("p = %.3f", T$p_welch))
T$panel <- factor(etiqueta[T$codigo], etiqueta); T$dia <- factor(T$dia, DIAS)
# margen superior para las etiquetas
techo <- data.frame(panel = factor(etiqueta, etiqueta), dia = factor("P30", DIAS),
                    y = lim$y.hi[match(VARS$codigo, lim$codigo)] + 0.22 * (lim$y.hi - lim$y.lo)[match(VARS$codigo, lim$codigo)])

pd <- position_dodge(width = 0.35)
g <- ggplot(M, aes(dia, m, colour = genotipo, group = genotipo)) +
  geom_point(data = D_all, aes(dia, y, colour = genotipo, shape = genotipo),
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.35, seed = 1),
             size = 1.6, alpha = 0.6, stroke = 0.7, inherit.aes = FALSE) +
  geom_line(position = pd, linewidth = 0.9) +
  geom_errorbar(aes(ymin = m - ee, ymax = m + ee), width = 0.15, position = pd, linewidth = 0.7) +
  geom_point(aes(shape = genotipo), fill = "white", size = 2.6, position = pd, stroke = 1.1) +
  geom_text(data = T, aes(dia, ypos, label = lab), inherit.aes = FALSE, size = 2.6, vjust = 0, lineheight = 0.9) +
  geom_blank(data = techo, aes(dia, y), inherit.aes = FALSE) +
  facet_wrap(~ panel, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = COL, labels = c(WT = "WT (n = 6; en P30, 5)", SOD1 = "SOD1 (n = 6)"), name = NULL) +
  scale_shape_manual(values = c(WT = 21, SOD1 = 22), labels = c(WT = "WT (n = 6; en P30, 5)", SOD1 = "SOD1 (n = 6)"), name = NULL) +
  labs(x = "Edad", y = "Ángulo de rodilla, cadera–rodilla–tobillo (°)\nmenor = más flexionada",
       title = "Ángulo de rodilla longitudinal, WT vs SOD1 (media ± EE; puntos = animales)",
       caption = paste("Modelo mixto y ~ genotipo × día + (1 | animal) en R (lme4/lmerTest). Por día: Welch; Holm sobre 6 días.",
                       "Media de todos los ciclos válidos, sin suavizado.")) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 9),
        legend.position = c(0.01, 0.02), legend.justification = c(0, 0), legend.background = element_blank(),
        panel.grid.major.y = element_line(colour = "#e6e6e3", linewidth = 0.4),
        plot.caption = element_text(colour = "#6b6b68", size = 7, hjust = 0), plot.title.position = "plot")

ggsave(file.path(DIR_OUT, "Fig_rodilla_longitudinal_R.png"), g, width = 15, height = 5.2, dpi = 300, bg = "white")
tryCatch(ggsave(file.path(DIR_OUT, "Fig_rodilla_longitudinal_R.pdf"), g, width = 15, height = 5.2, device = cairo_pdf),
         error = function(e) ggsave(file.path(DIR_OUT, "Fig_rodilla_longitudinal_R.pdf"), g, width = 15, height = 5.2))
cat("\nListo. Resultados en:", DIR_OUT, "\n")
