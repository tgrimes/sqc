# Example datasets for the sqc package. Run from the package root:  Rscript data-raw/make-data.R
# Homework datasets are deliberately not included here.

# ---- filler: four-head filler, mixed-head sampling plan --------------------
# 40 sampled cycles; one bottle per head each. Target 500 mL, within-head SD 1 mL,
# head 3 runs 3 mL (= 3 SD) high.
set.seed(6)
head_mean <- c(500, 500, 503, 500)
filler <- expand.grid(cycle = 1:40, head = 1:4)
filler$fill_ml <- round(rnorm(nrow(filler), head_mean[filler$head], 1), 1)
filler <- filler[order(filler$cycle, filler$head), ]
rownames(filler) <- NULL
attr(filler, "out.attrs") <- NULL   # left over from expand.grid()

# ---- filler_rotating: same process, rational subgroup plan ------------------
# 5 consecutive cycles from head 1, then 5 from head 2, head 3, head 4, then head 1 again
# (100 sampled bottles, 20 subgroups of 5).
set.seed(6)
filler_rotating <- data.frame(cycle = 1:100)
filler_rotating$head <- ((filler_rotating$cycle - 1) %/% 5) %% 4 + 1
filler_rotating$fill_ml <- round(rnorm(100, head_mean[filler_rotating$head], 1), 1)

save(filler, file = "data/filler.rda", compress = "xz", version = 2)
save(filler_rotating, file = "data/filler_rotating.rda", compress = "xz", version = 2)
