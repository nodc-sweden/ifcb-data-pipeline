#' Flag samples dominated by very large particles (air bubbles)
#'
#' Supplements the Bubbles flag in iRfcb::ifcb_psd(), which only fires when the
#' PSD peak is shifted to large sizes. In sparse bubble samples the peak can
#' land in the 0 um bin (ROIs without segmented blobs), so that flag is missed.
#' This rule instead flags samples where a large share of particles are big.
#'
#' @param psd_data Data frame from ifcb_psd()$data (columns `0um`..`199um`
#'   in particles per L, plus `mL_analyzed`).
#' @param min_esd Numeric. ESD (um) at or above which a particle counts as big.
#' @param min_frac Numeric. Minimum fraction of sized particles that must be big.
#' @param min_n Numeric. Minimum number of big particles.
#'
#' @return Tibble with columns `sample`, `n_big` and `frac_big` for flagged samples.
flag_psd_bubbles <- function(psd_data, min_esd = 100, min_frac = 0.2, min_n = 5) {
  um <- grep("um$", names(psd_data), value = TRUE)
  esd <- as.integer(sub("um$", "", um))

  # Convert concentrations (per L) back to particle counts
  counts <- round(as.matrix(psd_data[, um]) * psd_data$mL_analyzed / 1000)
  counts[is.na(counts)] <- 0

  n_big <- rowSums(counts[, esd >= min_esd, drop = FALSE])

  # Exclude the 0 um bin (ROIs without a segmented blob)
  n_valid <- rowSums(counts[, esd > 0, drop = FALSE])

  dplyr::tibble(sample = psd_data$sample,
                n_big = n_big,
                frac_big = n_big / pmax(n_valid, 1)) %>%
    dplyr::filter(n_big >= min_n, frac_big >= min_frac)
}

#' Set flag to "Bubbles" for the given samples, replacing any existing flag
#'
#' @param flags Data frame from ifcb_psd()$flags (may be NULL).
#' @param bubble_samples Character vector of sample names to flag.
#'
#' @return Tibble with columns `sample` and `flag`.
add_bubble_flags <- function(flags, bubble_samples) {
  if (is.null(flags)) {
    flags <- dplyr::tibble(sample = character(), flag = character())
  }

  flags %>%
    dplyr::filter(!sample %in% bubble_samples) %>%
    dplyr::bind_rows(dplyr::tibble(sample = bubble_samples, flag = "Bubbles")) %>%
    dplyr::arrange(sample)
}
