library(tidyverse)
library(worrms)

#' Match class names to WoRMS taxonomy
#'
#' @param class_names Character vector of class names
#' @return A tibble with columns: class_name, cleaned_name, reported_name, aphia_id, flag
clean_taxa <- function(class_names) {
  class_names <- class_names |>
    str_squish() |>
    keep(\(x) nchar(x) > 0)

  classlist <- tibble(class_name = class_names)

  # Non-taxonomic classes to skip
  non_taxa <- c("Air_bubbles", "Beads", "Debris", "mixed")

  # Manual overrides: class_name -> (reported_name, aphia_id)
  manual_overrides <- list(
    Akinete = list(reported_name = "Cyanophyceae", aphia_id = 146542L),
    Heterocyte = list(reported_name = "Cyanophyceae", aphia_id = 146542L),
    Echinoderm_larvae = list(reported_name = "Echinodermata", aphia_id = 1806L),
    Mussel_larvae = list(reported_name = "Bivalvia", aphia_id = 105L),
    Zooplankton = list(reported_name = "Crustacea", aphia_id = 1066L),
    Actinocyclus_spp = list(reported_name = "Actinocyclus_spp", aphia_id = 148944L),
    Unicells = list(reported_name = "Biota", aphia_id = 1L)
  )

  # --- Cleaning helpers ---

  clean_taxon <- function(name) {
    name |>
      str_replace_all("_", " ") |>
      str_remove("\\s+(smaller|larger)\\s+than\\s+\\d+") |>
      str_remove("\\s+(small|large)$") |>
      str_replace("(\\D)(\\d+)$", "\\1") |>
      str_remove("\\s+(CC|CS|Eliptical|pair)$") |>
      str_remove("\\s+(bundle|filament|coil|chain|single cell|single|transformation|heterotropic)$") |>
      str_replace("^Tintinnids$", "Tintinnida") |>
      str_squish()
  }

  detect_flag <- function(cleaned) {
    if (str_detect(cleaned, "\\bgroup$")) return("GRP")
    if (str_detect(cleaned, "\\bcf\\b")) return("CF")
    if (str_detect(cleaned, "-like$")) return("CF")
    if (str_detect(cleaned, "\\bspp$")) return("SPP")
    if (str_detect(cleaned, "\\bsp$")) return("SP")
    NA_character_
  }

  detect_type <- function(cleaned) {
    if (str_detect(cleaned, "\\bgroup$")) return("group")
    if (str_detect(cleaned, "\\bcf\\b")) return("cf")
    if (str_detect(cleaned, "-like$")) return("like")
    parts <- str_split(cleaned, "\\s+")[[1]]
    caps <- which(str_detect(parts, "^[A-Z]"))
    if (length(caps) > 1 && caps[1] == 1) {
      if (any(caps[-1] > 1)) return("multi")
    }
    if (str_detect(cleaned, "^[A-Z][a-z]+-[A-Z][a-z]+$")) return("multi")
    if (str_detect(cleaned, "^[A-Z][a-z]+\\s+[A-Z][a-z]+$")) return("multi")
    "direct"
  }

  extract_queries <- function(cleaned, type) {
    if (type == "cf") {
      return(str_replace(cleaned, "\\s+cf\\s+", " "))
    }
    if (type == "like") {
      return(str_remove(cleaned, "-like$"))
    }
    if (type == "multi") {
      if (str_detect(cleaned, "^[A-Z][a-z]+-[A-Z][a-z]+$")) {
        return(str_split(cleaned, "-")[[1]])
      }
      parts <- str_split(cleaned, "\\s+")[[1]]
      taxa <- list()
      current <- c()
      for (part in parts) {
        if (length(current) > 0 && str_detect(part, "^[A-Z]")) {
          taxa <- c(taxa, list(paste(current, collapse = " ")))
          current <- part
        } else {
          current <- c(current, part)
        }
      }
      taxa <- c(taxa, list(paste(current, collapse = " ")))
      return(unlist(taxa))
    }
    cleaned |>
      str_remove("\\s+spp$") |>
      str_remove("\\s+sp$")
  }

  # --- WoRMS query helpers ---

  worms_cache <- list()


  empty_worms <- list(
    aphia_id = NA_integer_, name = NA_character_, rank = NA_character_,
    worms_kingdom = NA_character_, worms_phylum = NA_character_,
    worms_class = NA_character_, worms_order = NA_character_,
    worms_family = NA_character_, worms_genus = NA_character_
  )

  query_worms_single <- function(name) {
    name <- str_squish(name)
    if (name %in% names(worms_cache)) return(worms_cache[[name]])

    Sys.sleep(0.3)
    cat("  Querying:", name, "\n")

    result <- tryCatch({
      matches <- wm_records_names(name, fuzzy = TRUE, marine_only = FALSE)[[1]]
      if (nrow(matches) > 0) {
        m <- matches[1, ]
        list(
          aphia_id = m$AphiaID, name = m$scientificname, rank = m$rank,
          worms_kingdom = m$kingdom, worms_phylum = m$phylum,
          worms_class = m$class, worms_order = m$order,
          worms_family = m$family, worms_genus = m$genus
        )
      } else {
        empty_worms
      }
    }, error = \(e) {
      cat("    Error:", conditionMessage(e), "\n")
      empty_worms
    })

    worms_cache[[name]] <<- result
    result
  }

  get_parent <- function(aphia_id) {
    if (is.na(aphia_id)) return(list(aphia_id = NA_integer_, name = NA_character_))

    Sys.sleep(0.3)
    tryCatch({
      cls <- wm_classification(aphia_id)
      if (nrow(cls) >= 2) {
        parent_row <- cls[nrow(cls) - 1, ]
        list(aphia_id = parent_row$AphiaID, name = parent_row$scientificname)
      } else {
        list(aphia_id = NA_integer_, name = NA_character_)
      }
    }, error = \(e) {
      cat("    Error getting parent for", aphia_id, ":", conditionMessage(e), "\n")
      list(aphia_id = NA_integer_, name = NA_character_)
    })
  }

  find_common_parent <- function(aphia_ids) {
    aphia_ids <- aphia_ids[!is.na(aphia_ids)]
    if (length(aphia_ids) == 0) return(list(aphia_id = NA_integer_, name = NA_character_))
    if (length(aphia_ids) == 1) return(get_parent(aphia_ids[1]))

    classifications <- map(aphia_ids, \(id) {
      Sys.sleep(0.3)
      tryCatch(wm_classification(id), error = \(e) NULL)
    })
    classifications <- compact(classifications)
    if (length(classifications) == 0) return(list(aphia_id = NA_integer_, name = NA_character_))

    id_lists <- map(classifications, \(cls) cls$AphiaID)
    common <- reduce(id_lists, intersect)

    if (length(common) == 0) return(list(aphia_id = NA_integer_, name = NA_character_))

    first_cls <- classifications[[1]]
    common_rows <- first_cls |> filter(AphiaID %in% common)
    deepest <- common_rows[nrow(common_rows), ]

    list(aphia_id = deepest$AphiaID, name = deepest$scientificname)
  }

  # --- Main processing ---

  as_taxonomy_tibble <- function(w, reported_name = w$name, reported_aphia_id = w$aphia_id) {
    tibble(
      reported_name = reported_name, aphia_id = reported_aphia_id,
      rank = w$rank, worms_kingdom = w$worms_kingdom, worms_phylum = w$worms_phylum,
      worms_class = w$worms_class, worms_order = w$worms_order,
      worms_family = w$worms_family, worms_genus = w$worms_genus
    )
  }

  process_taxon <- function(cleaned_raw, type) {
    queries <- extract_queries(cleaned_raw, type)

    if (type == "direct") {
      w <- query_worms_single(queries)
      return(as_taxonomy_tibble(w))
    }

    if (type %in% c("cf", "like")) {
      w <- query_worms_single(queries)
      if (!is.na(w$aphia_id)) {
        cat("    Getting parent for", w$name, "(", type, ")\n")
        parent <- get_parent(w$aphia_id)
        return(as_taxonomy_tibble(w,
          reported_name = parent$name,
          reported_aphia_id = parent$aphia_id
        ))
      }
      return(as_taxonomy_tibble(empty_worms))
    }

    if (type == "multi") {
      cat("  Multi-taxa:", paste(queries, collapse = " + "), "\n")
      worms_hits <- map(queries, query_worms_single)
      aphia_ids <- map_int(worms_hits, \(w) as.integer(w$aphia_id))

      cat("    Finding common parent...\n")
      parent <- find_common_parent(aphia_ids)
      # Use taxonomy from the first hit
      first_hit <- worms_hits[[1]]
      return(as_taxonomy_tibble(first_hit,
        reported_name = parent$name,
        reported_aphia_id = parent$aphia_id
      ))
    }
  }

  result <- classlist |>
    mutate(
      is_taxon = !class_name %in% non_taxa,
      cleaned_raw = map_chr(class_name, clean_taxon),
      flag = map_chr(cleaned_raw, detect_flag),
      type = map_chr(cleaned_raw, detect_type),
      cleaned = cleaned_raw |>
        str_remove("\\s+group$") |>
        str_replace("\\s+cf\\s+", " ") |>
        str_remove("-like$") |>
        str_remove("\\s+spp$") |>
        str_remove("\\s+sp$") |>
        str_squish()
    )

  cat("Processing", sum(result$is_taxon), "taxonomic entries...\n")
  cat("Types:", paste(table(result$type[result$is_taxon]), collapse = ", "), "\n\n")

  rows <- list()
  for (i in seq_len(nrow(result))) {
    cn <- result$class_name[i]
    if (cn %in% names(manual_overrides)) {
      ov <- manual_overrides[[cn]]
      # Look up taxonomy for the override's aphia_id
      w <- query_worms_single(ov$reported_name)
      rows[[i]] <- as_taxonomy_tibble(w,
        reported_name = ov$reported_name,
        reported_aphia_id = ov$aphia_id
      )
    } else if (!result$is_taxon[i] || result$type[i] == "group") {
      rows[[i]] <- as_taxonomy_tibble(empty_worms)
    } else {
      rows[[i]] <- process_taxon(result$cleaned_raw[i], result$type[i])
    }
  }

  out <- bind_rows(rows)

  final <- bind_cols(result |> select(class_name, cleaned_name = cleaned, is_taxon, flag, type), out) |>
    mutate(
      species_flag_taxon_name = case_when(
        is.na(flag) ~ NA_character_,
        flag == "SPP" ~ paste(cleaned_name, "spp."),
        flag == "GRP" ~ paste(cleaned_name, "group"),
        flag == "CF" & rank == "Species" ~ sub(" ", " cf. ", cleaned_name),
        flag == "CF" & rank == "Genus" ~ paste0(worms_order, ": ", worms_family, " cf. ", worms_genus),
        TRUE ~ cleaned_name
      )
    ) |>
    select(
      class_name, cleaned_name, reported_name, aphia_id, flag, rank,
      worms_kingdom, worms_phylum, worms_class, worms_order,
      worms_family, worms_genus, species_flag_taxon_name
    )

  cat("\nDone!\n")
  cat("Matched:", sum(!is.na(final$aphia_id)), "/", sum(result$is_taxon), "taxonomic entries\n")

  unmatched <- final |> filter(result$is_taxon & is.na(aphia_id))
  if (nrow(unmatched) > 0) {
    cat("\nUnmatched taxa:\n")
    print(unmatched, n = Inf)
  }

  final
}
