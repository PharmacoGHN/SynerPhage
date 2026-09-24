# Reading whatever table the user has, discovering its treatments and folding
# it into the wide (strain x replicate x treatment) layout the engine works on.
#
# Nothing here hardcodes phage or strain names. Treatments are discovered from
# the header: a combination is any column whose name splits into components
# that all exist as single-treatment columns, whatever the separator, whatever
# the number of components.

box::use(
  readr[read_delim],
  readxl[excel_sheets, read_excel],
  stats[ave],
  tidyr[all_of, pivot_wider],
  tools[file_ext],
)
box::use(
  app/logic/postprocessing[clamp_probability, clamp_report],
)

# Simulated example shipped with the app (see app/logic/mock_data.R). It sits in
# app/static so the sidebar can offer it as a plain download link.
#' @export
example_file <- box::file("../static/example_mock.csv")

# Columns that are experimental controls rather than treatments, and columns a
# spreadsheet export leaves behind. Both are only defaults: the user can undo
# either choice in the sidebar.
#' @export
control_pattern <- "^(no.?phages?|control|ctrl|blank|untreated|mock|all.?phages)$"
index_pattern <- "^(|\\.{3}[0-9]+|x|unnamed.*)$"

# --- reading whatever the user has -----------------------------------------

#' @export
file_type <- function(name) tolower(file_ext(name))

clean_frame <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- trimws(names(df))
  df[!grepl(index_pattern, names(df), ignore.case = TRUE)]
}

# Excel workbooks and .RData files hold several tables; the caller offers the
# list to the user and passes the chosen one back as `item`.
#' @export
source_items <- function(path, ext) {
  if (ext %in% c("xlsx", "xls")) {
    excel_sheets(path)
  } else if (ext %in% c("rdata", "rda")) {
    e <- new.env()
    load(path, e)
    ls(e)
  } else {
    character(0)
  }
}

#' @export
read_source <- function(path, ext, item = NULL) {
  pick <- function(x) if (!is.null(item) && nzchar(item) && item %in% x) item else x[[1]]
  df <- switch(ext,
    xlsx = ,
    xls = read_excel(path, sheet = pick(excel_sheets(path))),
    rdata = ,
    rda = {
      e <- new.env()
      load(path, e)
      get(pick(ls(e)), e)
    },
    rds = readRDS(path),
    # csv / tsv / txt: readr guesses the delimiter.
    read_delim(path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal")
  )
  clean_frame(df)
}

#' @export
read_pasted <- function(txt) {
  clean_frame(read_delim(I(txt), show_col_types = FALSE, progress = FALSE,
                         name_repair = "minimal"))
}

parse_values <- function(txt) {
  v <- suppressWarnings(as.numeric(strsplit(trimws(txt), "[,;[:space:]]+")[[1]]))
  v[!is.na(v)]
}

#' @export
guess_col <- function(cols, pattern) {
  hit <- cols[grepl(pattern, cols, ignore.case = TRUE)]
  if (length(hit)) hit[[1]] else ""
}

# --- treatment discovery ---------------------------------------------------

# A column is a combination when its name splits into >= 2 components that are
# all present as single-treatment columns. Controls such as "all_phages" or
# "Nophage" drop out on their own: their parts match nothing.
#' @export
split_treatments <- function(cols, sep) {
  parts <- lapply(strsplit(cols, sep, fixed = TRUE), trimws)
  names(parts) <- cols
  singles <- cols[lengths(parts) == 1L]
  is_combo <- lengths(parts) > 1L &
    vapply(parts, function(p) all(p %in% singles), logical(1))
  list(singles = singles, combos = parts[is_combo])
}

# Pick the separator that resolves the most combinations, so that
# "V1SA9+V1SA12" and "v1sa9_v1sa12" both work without the user being asked.
#' @export
detect_separator <- function(cols) {
  cand <- c("+", "_", "-", ":", "&", "/", "|", ";")
  n <- vapply(cand, function(s) length(split_treatments(cols, s)$combos), integer(1))
  if (max(n) == 0L) "+" else cand[[which.max(n)]]
}

# --- dataset preparation ---------------------------------------------------

# Long / tidy input (one row per measurement) folded into the wide layout the
# rest of the engine works on.
#' @export
widen_long <- function(df, strain_col, rep_col, treat_col, value_col) {
  df <- clean_frame(df)
  for (nm in c(treat_col, value_col)) {
    if (!nzchar(nm) || !nm %in% names(df)) stop("Pick a treatment column and a value column.")
  }
  keep <- unique(c(strain_col, rep_col, treat_col, value_col))
  out <- df[keep[nzchar(keep)]]
  out[[value_col]] <- suppressWarnings(as.numeric(as.character(out[[value_col]])))
  out[[treat_col]] <- trimws(as.character(out[[treat_col]]))
  as.data.frame(
    pivot_wider(out, names_from = all_of(treat_col),
                values_from = all_of(value_col), values_fn = mean),
    check.names = FALSE
  )
}

#' @export
prepare_dataset <- function(df, strain_col = "", rep_col = "", sep = "",
                            ignore = character(), scale = "auto") {
  df <- clean_frame(df)

  # A Shiny input that has not been rendered yet arrives as NULL; fall back to
  # the default rather than failing on a zero-length comparison.
  or_default <- function(x, default) if (length(x) == 1L && !is.na(x)) x else default
  strain_col <- or_default(strain_col, "")
  rep_col <- or_default(rep_col, "")
  sep <- or_default(sep, "")
  scale <- or_default(scale, "auto")

  strain <- if (nzchar(strain_col) && strain_col %in% names(df)) {
    as.character(df[[strain_col]])
  } else {
    rep("sample", nrow(df))
  }
  strain[is.na(strain) | !nzchar(strain)] <- "sample"

  effect_cols <- setdiff(names(df), c(strain_col, rep_col, ignore))
  values <- suppressWarnings(lapply(df[effect_cols], function(x) as.numeric(as.character(x))))
  values <- values[vapply(values, function(v) any(is.finite(v)), logical(1))]
  if (length(values) < 3) {
    stop("Need at least three numeric treatment columns: two single phages and one combination.")
  }

  # Percent-scaled effects are rescaled to the 0-1 support the Bliss model needs.
  top <- suppressWarnings(max(unlist(values), na.rm = TRUE))
  as_percent <- scale == "percent" || (scale == "auto" && is.finite(top) && top > 1.5)
  if (as_percent) values <- lapply(values, function(v) v / 100)

  eff <- as.data.frame(values, check.names = FALSE)
  clamped <- clamp_report(as.matrix(eff))
  eff[] <- clamp_probability(as.matrix(eff))

  # Order rows by strain (stable, so a replicate column only refines the order)
  # and renumber replicates 1..k: the JAGS matrices are filled row by row.
  ord <- if (nzchar(rep_col) && rep_col %in% names(df)) {
    order(strain, df[[rep_col]])
  } else {
    order(strain)
  }
  strain <- strain[ord]
  data <- data.frame(bacteria = strain,
                     replica = ave(seq_along(strain), strain, FUN = seq_along),
                     eff[ord, , drop = FALSE], check.names = FALSE)

  sep <- if (nzchar(sep)) sep else detect_separator(names(eff))
  found <- split_treatments(names(eff), sep)
  if (!length(found$combos)) {
    stop(sprintf("No combination column found with separator '%s'. Expected names such as %s.",
                 sep, paste0("phageA", sep, "phageB")))
  }

  list(data = data, singles = found$singles, combos = found$combos, sep = sep,
       n_rep = max(data$replica), strains = unique(data$bacteria),
       rescaled = as_percent, clamped = clamped)
}

# One row per (strain, phage pair) typed straight into the sidebar. It goes
# through exactly the same engine as an uploaded table.
#' @export
pair_dataset <- function(strain, name_a, name_b, values_a, values_b, values_ab) {
  vals <- list(parse_values(values_a), parse_values(values_b), parse_values(values_ab))
  if (any(lengths(vals) == 0L)) stop("Enter a value in each of the three effect fields.")
  n <- max(lengths(vals))
  if (!all(lengths(vals) %in% c(1L, n))) {
    stop("Give either one value or the same number of replicates in each field.")
  }
  vals <- lapply(vals, rep_len, n)

  safe <- function(x, fallback) {
    x <- trimws(gsub("[+_&|:;/]", "-", x))
    if (nzchar(x)) x else fallback
  }
  nm <- make.unique(c(safe(name_a, "Phage A"), safe(name_b, "Phage B")))

  df <- data.frame(bacteria = rep(safe(strain, "sample"), n), replica = seq_len(n),
                   check.names = FALSE)
  df[[nm[1]]] <- vals[[1]]
  df[[nm[2]]] <- vals[[2]]
  df[[paste(nm, collapse = " + ")]] <- vals[[3]]
  df
}
