#' Stratigraphic Plug Alignment (SPA)
#'
#' Linearly interpolates plug-based measurements (e.g., TOC, porosity, XRD)
#' onto a high-resolution reference depth grid (e.g., XRF). The procedure
#' uses base R's \code{approx()} with \code{rule = 1} to prevent extrapolation
#' beyond the observed depth range, ensuring stratigraphically consistent
#' alignment.
#'
#' @param ref A data.frame containing the reference depth grid and (optionally)
#'   high-resolution variables (e.g., XRF). Must include \code{depth_col}.
#' @param ... One or more named data.frames containing plug-based measurements
#'   to be interpolated (e.g., \code{xrd = xrd_df}, \code{plugs = plug_df}).
#' @param depth_col Character string with the name of the depth column shared
#'   by all input datasets.
#' @param rule Integer passed to \code{approx()} (default = 1).
#' @param add_suffix Logical; TRUE = suffix variable names with dataset name.
#'
#' @return A data.frame with the reference depth grid and interpolated variables.
#' @export
spa_align <- function(ref,
                      ...,
                      depth_col = "Depth_m",
                      rule = 1,
                      add_suffix = TRUE) {
  plug_list <- list(...)
  if (length(plug_list) == 0L) {
    stop("Provide datasets to interpolate, e.g., xrd = xrd_df")
  }

  if (!depth_col %in% names(ref)) stop("Depth column not in ref dataset.")

  ref <- ref[order(ref[[depth_col]]), , drop = FALSE]
  ref_depth <- ref[[depth_col]]
  out <- ref

  for (nm in names(plug_list)) {
    df <- plug_list[[nm]]

    if (!is.data.frame(df)) stop("Dataset '", nm, "' is not a data.frame.")
    if (!depth_col %in% names(df)) stop("Depth column missing in dataset '", nm, "'.")

    df <- df[order(df[[depth_col]]), , drop = FALSE]

    numeric_cols <- vapply(df, is.numeric, logical(1))
    numeric_cols[depth_col] <- FALSE
    vars <- names(df)[numeric_cols]

    for (v in vars) {
      interp_vals <- stats::approx(
        x    = df[[depth_col]],
        y    = df[[v]],
        xout = ref_depth,
        rule = rule
      )$y

      new_name <- if (add_suffix) paste0(v, "_", nm) else v
      out[[new_name]] <- interp_vals
    }
  }

  out
}

