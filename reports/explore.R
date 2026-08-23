## Use this script to run exploratory code maybe before to put it into
## the pipeline


# setup -----------------------------------------------------------
library(targets)
library(here)

# load all your custom functions
tar_source()


# Code here below -------------------------------------------------
# use `tar_read(target_name)` to load a target anywhere (note that
# `target_name` is NOT quoted!)

# =====================================================================
# Mini test script for preprocess_cmr()
# Self-contained: builds the two datasets inline, then calls the function.
# =====================================================================

library(tidyverse)
library(janitor)


# ---------------------------------------------------------------------
# 1. mort_count: one row per death, 9 towns with deaths (town 15016 has none)
#    'Ischaemic heart diseases' carries weight 0.5 (a 50/50 split cause).
# ---------------------------------------------------------------------
mort_count <- tibble::tribble(
  ~comune_residenza_2023, ~sesso, ~eta, ~cause,                      ~group,                                 ~mechanism,                      ~weight,
  "15002", 2L, 35L, "Lung cancer",              "Cancer",                             "Lifestyle and NCDs",           1.0,
  "15002", 1L, 67L, "Ischaemic heart diseases", "Diseases of the circulatory system", "Lifestyle and NCDs",           0.5,
  "15002", 1L, 57L, "Lung cancer",              "Cancer",                             "Lifestyle and NCDs",           1.0,
  "15002", 1L, 21L, "Transport Accidents",      "Injuries",                           "Environment and Safety",       1.0,
  "15005", 1L, 52L, "Ischaemic heart diseases", "Diseases of the circulatory system", "Lifestyle and NCDs",           0.5,
  "15005", 1L, 65L, "Lung cancer",              "Cancer",                             "Lifestyle and NCDs",           1.0,
  "15005", 1L, 48L, "Transport Accidents",      "Injuries",                           "Environment and Safety",       1.0,
  "15007", 2L, 70L, "Measles",                  "Infectious diseases",                "Immunisation and Prophylaxis", 1.0,
  "15007", 1L, 44L, "Ischaemic heart diseases", "Diseases of the circulatory system", "Lifestyle and NCDs",           0.5,
  "15007", 2L, 61L, "Lung cancer",              "Cancer",                             "Lifestyle and NCDs",           1.0,
  "15009", 1L, 33L, "Transport Accidents",      "Injuries",                           "Environment and Safety",       1.0,
  "15009", 2L, 58L, "Measles",                  "Infectious diseases",                "Immunisation and Prophylaxis", 1.0,
  "15009", 1L, 72L, "Ischaemic heart diseases", "Diseases of the circulatory system", "Lifestyle and NCDs",           0.5,
  "15010", 2L, 41L, "Lung cancer",              "Cancer",                             "Lifestyle and NCDs",           1.0,
  "15010", 1L, 69L, "Ischaemic heart diseases", "Diseases of the circulatory system", "Lifestyle and NCDs",           0.5,
  "15010", 1L, 27L, "Transport Accidents",      "Injuries",                           "Environment and Safety",       1.0,
  "15011", 2L, 63L, "Measles",                  "Infectious diseases",                "Immunisation and Prophylaxis", 1.0,
  "15011", 1L, 55L, "Lung cancer",              "Cancer",                             "Lifestyle and NCDs",           1.0,
  "15012", 1L, 38L, "Ischaemic heart diseases", "Diseases of the circulatory system", "Lifestyle and NCDs",           0.5,
  "15012", 2L, 74L, "Lung cancer",              "Cancer",                             "Lifestyle and NCDs",           1.0,
  "15012", 1L, 29L, "Transport Accidents",      "Injuries",                           "Environment and Safety",       1.0,
  "15014", 2L, 66L, "Measles",                  "Infectious diseases",                "Immunisation and Prophylaxis", 1.0,
  "15014", 1L, 50L, "Ischaemic heart diseases", "Diseases of the circulatory system", "Lifestyle and NCDs",           0.5,
  "15015", 1L, 45L, "Lung cancer",              "Cancer",                             "Lifestyle and NCDs",           1.0,
  "15015", 2L, 71L, "Transport Accidents",      "Injuries",                           "Environment and Safety",       1.0,
  "15015", 1L, 60L, "Ischaemic heart diseases", "Diseases of the circulatory system", "Lifestyle and NCDs",           0.5
)

# ---------------------------------------------------------------------
# 2. population: like the real ';' file but already read in.
#    Codice comune is NUMERIC (tests the as.character() coercion),
#    two years (tests the pop_year filter), several age bands (summed),
#    and all 10 towns - including 15016, which has no deaths.
# ---------------------------------------------------------------------
population <- tidyr::expand_grid(
  `Codice comune` = c(15002L, 15005L, 15007L, 15009L, 15010L,
                      15011L, 15012L, 15014L, 15015L, 15016L),
  Eta             = c("000", "040", "070"),
  anno            = c(2022L, 2023L),
  sesso           = c(1L, 2L)
) |>
  dplyr::mutate(
    Comune = paste0("town_", `Codice comune`),
    # give 2023 a slightly different size than 2022 so the year filter matters
    numero = 500L + (`Codice comune` %% 50L) * 10L + (anno - 2022L) * 25L
  )

# ---------------------------------------------------------------------
# 3. Run it
# ---------------------------------------------------------------------
res <- preprocess_cmr(mort_count, population)   # population as a data frame
# res <- preprocess_cmr(mort_count, "test_population.csv")  # or from a CSV path

print(res, width = Inf)

# ---------------------------------------------------------------------
# 4. Quick sanity checks
# ---------------------------------------------------------------------
stopifnot(
  nrow(res) == 10,                              # all towns kept
  all(c("comune", "population", "total") %in% names(res)),
  # town 15016 has no deaths -> every rate column is 0
  res |> dplyr::filter(comune == "15016") |>
    dplyr::select(-comune, -population) |> unlist() |> sum() == 0
)

# the C_, G_ and M_ blocks should each reproduce `total` per comune
res |>
  dplyr::mutate(
    sum_C = rowSums(dplyr::across(dplyr::starts_with("C_"))),
    sum_G = rowSums(dplyr::across(dplyr::starts_with("G_"))),
    sum_M = rowSums(dplyr::across(dplyr::starts_with("M_")))
  ) |>
  dplyr::select(comune, population, total, sum_C, sum_G, sum_M)



# =====================================================================
# Mini test script for preprocess_smr()  (indirect age-sex standardisation)
# Self-contained: builds the two datasets inline, then calls the function.
#
# Strata are SINGLE YEARS OF AGE x SEX, so both the deaths (eta, sesso) and the
# population (Eta single-year, sesso) must be stratified that way.
# =====================================================================

library(tidyverse)
library(janitor)
set.seed(1)

# 5 towns, 6-digit ISTAT codes in the DEATHS (leading zero), 5-digit in the POP
towns <- tibble::tibble(
  code6 = c("015002", "015005", "015007", "015009", "015016"),  # deaths form
  code5 = c( 15002L,   15005L,   15007L,   15009L,   15016L),    # population form (numeric)
  name  = c("Abbiategrasso", "Albairate", "Arconate", "Assago", "Besate")
)

# 3 causes spanning the axes; IHD is a 50/50 split cause (weight 0.5)
causes <- tibble::tribble(
  ~cause,                      ~group,                                 ~mechanism,                       ~weight,
  "Lung cancer",               "Cancer",                              "Lifestyle and NCDs",            1.0,
  "Ischaemic heart diseases",  "Diseases of the circulatory system",  "Lifestyle and NCDs",            0.5,
  "Measles",                   "Infectious diseases",                 "Immunisation and Prophylaxis",  1.0
)

# ---------------------------------------------------------------------
# 1. mort_count: one row per death. Towns 1-4 have deaths; town 015016 has none.
#    eta = single year of age, sesso = 1/2.
# ---------------------------------------------------------------------
mort_count <- purrr::map_dfr(towns$code6[1:4], function(cd) {
  n <- sample(6:14, 1)
  idx <- sample(seq_len(nrow(causes)), n, replace = TRUE)
  tibble::tibble(
    comune_residenza_2023 = cd,
    sesso = sample(1:2, n, replace = TRUE),
    eta   = sample(40:74, n, replace = TRUE),
    cause = causes$cause[idx],
    group = causes$group[idx],
    mechanism = causes$mechanism[idx],
    weight    = causes$weight[idx]
  )
})

# ---------------------------------------------------------------------
# 2. population: ';' layout, codes as NUMBERS, single-year Eta ("040".."074"),
#    two years (tests pop_year filter), both sexes, all 5 towns incl. 015016.
# ---------------------------------------------------------------------
population <- tidyr::expand_grid(
  `Codice comune` = towns$code5,
  age             = 40:74,
  anno            = c(2022L, 2023L),
  sesso           = c(1L, 2L)
) |>
  dplyr::left_join(towns, by = c(`Codice comune` = "code5")) |>
  dplyr::mutate(
    Eta    = sprintf("%03d", age),                       # single year, zero-padded
    Comune = name,
    # a smooth age profile so strata are populated; 2023 slightly bigger
    numero = as.integer(round(
      (250 - abs(age - 55) * 2) * (1 + (`Codice comune` %% 7) / 10) +
        (anno - 2022L) * 5
    ))
  ) |>
  dplyr::select(`Codice comune`, Eta, anno, sesso, Comune, numero)

# ---------------------------------------------------------------------
# 3. Run it
# ---------------------------------------------------------------------
res <- preprocess_smr(mort_count, population)   # population as a data frame
# res <- preprocess_smr(mort_count, "test_population.csv")  # or a CSV path

print(res, width = Inf)

# ---------------------------------------------------------------------
# 4. Sanity checks
# ---------------------------------------------------------------------
stopifnot(
  nrow(res) == 5,                                       # all towns kept
  all(c("comune", "population", "total_smr", "total_isr") %in% names(res)),
  # death-free town 015016: every standardised value is 0
  res |> dplyr::filter(comune == "015016") |>
    dplyr::select(-comune, -population) |> unlist() |> sum() == 0
)

# ---- KEY INVARIANT of indirect standardisation -----------------------
# Across all areas, sum(expected) == sum(observed). We must compute expected
# DIRECTLY from the standard schedule, NOT reconstruct it as observed / SMR:
# a death-free area has SMR = 0 yet a POSITIVE expected, and 0 / 0 would wrongly
# drop that area's expected and break the check.
#
# Recompute the standard schedule and expected deaths here, independently of the
# function, then compare totals.
pop_strata <- population |>
  dplyr::filter(anno == 2023) |>
  dplyr::transmute(
    .area = sprintf("%06d", as.integer(`Codice comune`)),
    .age  = as.integer(Eta),
    .sex  = as.integer(sesso),
    .pop  = as.numeric(numero)
  ) |>
  dplyr::group_by(.area, .age, .sex) |>
  dplyr::summarise(.pop = sum(.pop), .groups = "drop")

m_chk <- mort_count |>
  dplyr::transmute(
    .area = as.character(comune_residenza_2023),
    .age  = as.integer(eta),
    .sex  = as.integer(sesso),
    .w    = weight
  )

std_denom <- pop_strata |>
  dplyr::group_by(.age, .sex) |>
  dplyr::summarise(std_pop = sum(.pop), .groups = "drop")

std_rate <- m_chk |>
  dplyr::group_by(.age, .sex) |>
  dplyr::summarise(std_deaths = sum(.w), .groups = "drop") |>
  dplyr::left_join(std_denom, by = c(".age", ".sex")) |>
  dplyr::mutate(std_rate = std_deaths / std_pop)

total_expected <- pop_strata |>
  dplyr::left_join(std_rate, by = c(".age", ".sex")) |>
  dplyr::mutate(exp = .pop * std_rate) |>
  dplyr::summarise(expected = sum(exp, na.rm = TRUE)) |>
  dplyr::pull(expected)

total_observed <- sum(m_chk$.w)

cat("\nIndirect-standardisation invariant (all-cause):\n")
cat("  total observed deaths :", total_observed, "\n")
cat("  total expected deaths :", round(total_expected, 6),
    " (should equal observed)\n")

stopifnot(abs(total_expected - total_observed) < 1e-6)


# =====================================================================
# Mini test script for Moran I's calculation before and after BYM
# =====================================================================


library(geostan)
library(spdep)
set.seed(1)


sg <- sf::st_drop_geometry(tar_read(smr_geo))
smr_geo_bym2 <- tar_read(smr_geo_bym2)
C_matrix <- tar_read(C_matrix)

# 1. on the total SMR (quick, but noisy in small comuni)
geostan::moran_plot(sg$total_smr, C_matrix)          # visual + the I in the title

# 2. better: a variance-stabilised log relative risk, NA-guarded
lrr <- log((sg$total_obs + 0.5) / (sg$total_exp + 0.5))
geostan::moran_plot(lrr, C_matrix)

# Monte Carlo permutation test for a p-value (999 permutations)
lw <- spdep::mat2listw(as(C_matrix, "matrix") * 1, style = "W")  # *1 turns logical -> numeric
spdep::moran.mc(lrr, lw, nsim = 9999, zero.policy = TRUE)

# After fit
spdep::moran.mc(smr_geo_bym2$bym2_rr, lw, nsim = 999, zero.policy = TRUE)
geostan::moran_plot( sf::st_drop_geometry(smr_geo_bym2)$bym2_rr, C_matrix)          # visual + the I in the title
