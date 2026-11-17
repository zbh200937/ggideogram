#' Read paired chromosome alignments, structural rearrangements or synteny blocks
#'
#' Imports genome alignment and collinearity outputs directly into paired
#' one-based, closed coordinate tables ready for [geom_chrlink()]:
#' - `"paf"`: `minimap2` / `wfmash` / `MashMap` Pairwise mApping Format (`.paf`
#'   or `.paf.gz`). Target coordinates (`cols 6--9`) become `Chr1/Start1/End1`
#'   and query coordinates (`cols 1--4`) become `Chr2/Start2/End2`, converted
#'   from zero-based half-open `[start, end)` to one-based closed
#'   `[start + 1, end]`.
#' - `"syri"`: `SyRI` `.out` structural rearrangement tables (`SYN`, `INV`,
#'   `TRANS`, `INVTR`, `DUP`, `INVDP`, etc.), preserving one-based closed
#'   coordinates and normalizing inverted classes to `Orientation = "-"`.
#' - `"mcscanx"`: `MCScanX` `.collinearity` files paired with a gene coordinate
#'   table (`genes`, either a 4-column MCScanX `.gff` `[chr, gene, start, end]`
#'   or a [read_chr_features()] GFF3/GTF table). Set `collapse_blocks = TRUE`
#'   to summarize each collinear block into one ribbon interval.
#'
#' @param file Local alignment/synteny file path, optionally `.gz` compressed.
#' @param format Input format; `"auto"` inspects the file extension and header.
#' @param genes Gene coordinate data frame or file path required when
#'   `format = "mcscanx"`. Accepts a 4-column MCScanX `.gff` (`chr`, `gene`,
#'   `start`, `end`) or a GFF3/GTF table/file with `Chr`, `Start`, `End` and an
#'   identifier column (`ID`, `Gene`, or `gene_id`).
#' @param collapse_blocks Logical; for `"mcscanx"`, collapse individual gene
#'   pairs in each collinearity alignment block into block-level intervals.
#' @param min_length Minimum alignment/block span in base pairs (default `0`).
#' @param min_mapq Minimum mapping quality for `"paf"` records (default `0`).
#' @param types Character vector of `SyRI` structural types to retain when
#'   `format = "syri"`.
#' @param chr_map Optional named character vector `c(old_name = "new_name")`
#'   applied to both `Chr1` and `Chr2` before optional `chr_key()` construction.
#' @param target_genome,query_genome Optional genome identifiers; when supplied,
#'   `Chr1` and/or `Chr2` are converted into composite keys via [chr_key()].
#' @param target_assembly,query_assembly Optional assembly identifiers passed to
#'   [chr_key()] alongside `target_genome` / `query_genome`.
#' @param karyotype Optional karyotype or `ideogram_data` used to validate
#'   both endpoint chromosomes and interval bounds.
#'
#' @return A data frame with `Chr1`, `Start1`, `End1`, `Chr2`, `Start2`, `End2`,
#'   `Orientation` (`"+"` or `"-"`) and format-specific metadata columns.
#' @seealso [geom_chrlink()], [chr_key()], [read_chr_features()]
#' @export
read_chr_links <- function(
    file,
    format = c("auto", "paf", "syri", "mcscanx"),
    genes = NULL,
    collapse_blocks = FALSE,
    min_length = 0,
    min_mapq = 0,
    types = c("SYN", "INV", "TRANS", "INVTR", "DUP", "INVDP"),
    chr_map = NULL,
    target_genome = NULL,
    query_genome = NULL,
    target_assembly = NULL,
    query_assembly = NULL,
    karyotype = NULL) {
  format <- match.arg(format)
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    stopf("`file` must be one local file path.")
  }
  if (!file.exists(file)) {
    stopf("Link file `%s` does not exist.", file)
  }
  check_nonnegative_layout(min_length, "min_length")
  check_nonnegative_layout(min_mapq, "min_mapq")
  if (!is.logical(collapse_blocks) || length(collapse_blocks) != 1L || is.na(collapse_blocks)) {
    stopf("`collapse_blocks` must be TRUE or FALSE.")
  }

  con <- if (grepl("\\.gz$", file, ignore.case = TRUE)) gzfile(file, "rt") else base::file(file, "rt")
  on.exit(close(con))
  raw_lines <- readLines(con, warn = FALSE)

  if (format == "auto") {
    format <- detect_link_format(file, raw_lines)
  }

  out <- switch(
    format,
    paf = parse_paf_links(raw_lines, min_length = min_length, min_mapq = min_mapq),
    syri = parse_syri_links(raw_lines, types = types, min_length = min_length),
    mcscanx = parse_mcscanx_links(raw_lines, genes = genes,
      collapse_blocks = collapse_blocks, min_length = min_length)
  )

  if (!is.null(chr_map)) {
    if (!is.character(chr_map) || is.null(names(chr_map)) || anyNA(chr_map) ||
        anyNA(names(chr_map)) || anyDuplicated(names(chr_map)) ||
        any(!nzchar(names(chr_map))) || any(!nzchar(chr_map))) {
      stopf("`chr_map` must be a named character vector with unique non-empty names.")
    }
    if (nrow(out)) {
      idx1 <- match(out$Chr1, names(chr_map))
      out$Chr1[!is.na(idx1)] <- unname(chr_map[idx1[!is.na(idx1)]])
      idx2 <- match(out$Chr2, names(chr_map))
      out$Chr2[!is.na(idx2)] <- unname(chr_map[idx2[!is.na(idx2)]])
    }
  }

  if (!is.null(target_genome)) {
    if (!is.character(target_genome) || length(target_genome) != 1L ||
        is.na(target_genome) || !nzchar(target_genome)) {
      stopf("`target_genome` must be one non-empty string.")
    }
    if (nrow(out)) {
      out$Chr1 <- chr_key(target_genome, out$Chr1, assembly = target_assembly)
    }
  }
  if (!is.null(query_genome)) {
    if (!is.character(query_genome) || length(query_genome) != 1L ||
        is.na(query_genome) || !nzchar(query_genome)) {
      stopf("`query_genome` must be one non-empty string.")
    }
    if (nrow(out)) {
      out$Chr2 <- chr_key(query_genome, out$Chr2, assembly = query_assembly)
    }
  }

  if (!is.null(karyotype) && nrow(out)) {
    k <- as_ideogram_data(karyotype)$karyotype
    i1 <- match(out$Chr1, k$.chr)
    i2 <- match(out$Chr2, k$.chr)
    unknown <- unique(c(out$Chr1[is.na(i1)], out$Chr2[is.na(i2)]))
    if (length(unknown)) {
      stopf("Unknown link chromosomes: %s.", format_chr_rows(unknown))
    }
    if (any(out$Start1 - 1 < k$.start[i1] | out$End1 > k$.end[i1] |
            out$Start2 - 1 < k$.start[i2] | out$End2 > k$.end[i2])) {
      stopf("Link intervals are outside their karyotype chromosome bounds.")
    }
  }

  attr(out, "source_format") <- format
  attr(out, "coordinates") <- "1-based-closed"
  out
}

detect_link_format <- function(file, raw_lines) {
  clean <- sub("\\.gz$", "", file, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(clean))
  if (ext == "paf") return("paf")
  if (ext %in% c("collinearity", "collinear")) return("mcscanx")
  if (ext == "syri") return("syri")

  non_empty <- raw_lines[nzchar(trimws(raw_lines))]
  if (any(grepl("^##\\s*Alignment\\b", non_empty))) return("mcscanx")
  body <- non_empty[!grepl("^\\s*#", non_empty)]
  if (length(body)) {
    fields <- strsplit(body[1], "\t", fixed = TRUE)[[1]]
    if (length(fields) >= 12L && fields[5] %in% c("+", "-") &&
        all(grepl("^[0-9]+$", fields[c(2, 3, 4, 7, 8, 9, 10, 11, 12)]))) {
      return("paf")
    }
    if (length(fields) == 12L) return("syri")
  }
  stopf("Cannot infer link file format for `%s`; supply `format` explicitly.", basename(file))
}

parse_paf_links <- function(raw_lines, min_length = 0, min_mapq = 0) {
  lines <- raw_lines[nzchar(trimws(raw_lines)) & !grepl("^\\s*#", raw_lines)]
  empty_paf <- function() {
    data.frame(
      Chr1 = character(0), Start1 = numeric(0), End1 = numeric(0),
      Chr2 = character(0), Start2 = numeric(0), End2 = numeric(0),
      Orientation = character(0), Matches = numeric(0),
      Block_length = numeric(0), Identity = numeric(0),
      MapQ = numeric(0), Target_length = numeric(0), Query_length = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  if (!length(lines)) return(empty_paf())

  split_rows <- strsplit(lines, "\t", fixed = TRUE)
  widths <- lengths(split_rows)
  if (any(widths < 12L)) {
    stopf("PAF input requires at least 12 tab-separated columns (row %d has %d).",
          which(widths < 12L)[1], widths[widths < 12L][1])
  }
  mat <- do.call(rbind, lapply(split_rows, `[`, 1:12))
  qname <- validate_chr(mat[, 1], "PAF query name")
  tname <- validate_chr(mat[, 6], "PAF target name")
  strand <- mat[, 5]
  if (any(!strand %in% c("+", "-"))) {
    stopf("PAF column 5 (`strand`) must be `+` or `-`.")
  }

  num_cols <- c(2L, 3L, 4L, 7L, 8L, 9L, 10L, 11L, 12L)
  nums <- lapply(num_cols, function(j) suppressWarnings(as.numeric(mat[, j])))
  if (any(vapply(nums, function(v) anyNA(v) || any(!is.finite(v)) || any(v != floor(v)) || any(v < 0), logical(1)))) {
    stopf("PAF coordinate, length, match and MapQ columns must be non-negative integers.")
  }
  qlen <- nums[[1]]; qstart <- nums[[2]]; qend <- nums[[3]]
  tlen <- nums[[4]]; tstart <- nums[[5]]; tend <- nums[[6]]
  nmatch <- nums[[7]]; alen <- nums[[8]]; mapq <- nums[[9]]

  if (any(qstart >= qend | tstart >= tend)) {
    stopf("PAF alignments require `start < end` on both target and query.")
  }
  if (any(qend > qlen | tend > tlen)) {
    stopf("PAF alignment coordinates exceed reported sequence lengths.")
  }

  out <- data.frame(
    Chr1 = tname,
    Start1 = tstart + 1,
    End1 = tend,
    Chr2 = qname,
    Start2 = qstart + 1,
    End2 = qend,
    Orientation = strand,
    Matches = nmatch,
    Block_length = alen,
    Identity = ifelse(alen > 0, nmatch / alen, NA_real_),
    MapQ = mapq,
    Target_length = tlen,
    Query_length = qlen,
    stringsAsFactors = FALSE
  )
  keep <- out$Block_length >= min_length & out$MapQ >= min_mapq
  out <- out[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}

parse_syri_links <- function(raw_lines, types, min_length = 0) {
  if (!is.character(types) || !length(types) || anyNA(types)) {
    stopf("`types` must be a non-empty character vector of SyRI annotation types.")
  }
  lines <- raw_lines[nzchar(trimws(raw_lines)) & !grepl("^\\s*#", raw_lines)]
  empty_syri <- function() {
    data.frame(
      Chr1 = character(0), Start1 = numeric(0), End1 = numeric(0),
      Chr2 = character(0), Start2 = numeric(0), End2 = numeric(0),
      Orientation = character(0), Type = character(0),
      ID = character(0), Parent = character(0), Copy_status = character(0),
      stringsAsFactors = FALSE
    )
  }
  if (!length(lines)) return(empty_syri())

  split_rows <- strsplit(lines, "\t", fixed = TRUE)
  widths <- lengths(split_rows)
  if (any(widths != 12L)) {
    stopf("SyRI `.out` input requires 12 tab-separated columns (row %d has %d).",
          which(widths != 12L)[1], widths[widths != 12L][1])
  }
  mat <- do.call(rbind, split_rows)
  if (nrow(mat) && tolower(mat[1, 1]) %in% c("chra", "chr1", "refchr") &&
      is.na(suppressWarnings(as.numeric(mat[1, 2])))) {
    mat <- mat[-1, , drop = FALSE]
  }
  if (!nrow(mat)) return(empty_syri())

  keep_type <- mat[, 11] %in% types &
    mat[, 1] != "-" & mat[, 6] != "-" &
    mat[, 2] != "-" & mat[, 3] != "-" &
    mat[, 7] != "-" & mat[, 8] != "-"
  mat <- mat[keep_type, , drop = FALSE]
  if (!nrow(mat)) return(empty_syri())

  s1 <- suppressWarnings(as.numeric(mat[, 2]))
  e1 <- suppressWarnings(as.numeric(mat[, 3]))
  s2 <- suppressWarnings(as.numeric(mat[, 7]))
  e2 <- suppressWarnings(as.numeric(mat[, 8]))
  if (anyNA(c(s1, e1, s2, e2)) || any(!is.finite(c(s1, e1, s2, e2))) ||
      any(c(s1, e1, s2, e2) < 1) ||
      any(c(s1, e1, s2, e2) != floor(c(s1, e1, s2, e2)))) {
    stopf("SyRI coordinates for retained structural blocks must be positive integers.")
  }

  inv_types <- c("INV", "INVTR", "INVDP", "INVAL", "INVTRAL", "INVDPAL")
  inverted <- (mat[, 11] %in% inv_types) | ((s1 > e1) != (s2 > e2))
  start1 <- pmin(s1, e1); end1 <- pmax(s1, e1)
  start2 <- pmin(s2, e2); end2 <- pmax(s2, e2)

  out <- data.frame(
    Chr1 = validate_chr(mat[, 1], "SyRI reference chromosome"),
    Start1 = start1,
    End1 = end1,
    Chr2 = validate_chr(mat[, 6], "SyRI query chromosome"),
    Start2 = start2,
    End2 = end2,
    Orientation = ifelse(inverted, "-", "+"),
    Type = mat[, 11],
    ID = mat[, 9],
    Parent = mat[, 10],
    Copy_status = mat[, 12],
    stringsAsFactors = FALSE
  )
  span <- pmin(out$End1 - out$Start1 + 1, out$End2 - out$Start2 + 1)
  out <- out[span >= min_length, , drop = FALSE]
  rownames(out) <- NULL
  out
}

load_mcscanx_genes <- function(genes) {
  if (is.null(genes)) {
    stopf("`genes` is required when `format = \"mcscanx\"` (supply a 4-column MCScanX `.gff` or GFF3/GTF feature table).")
  }
  if (is.character(genes)) {
    if (length(genes) != 1L || is.na(genes) || !file.exists(genes)) {
      stopf("`genes` file path does not exist.")
    }
    con <- if (grepl("\\.gz$", genes, ignore.case = TRUE)) gzfile(genes, "rt") else base::file(genes, "rt")
    on.exit(close(con))
    glines <- readLines(con, warn = FALSE)
    glines <- glines[nzchar(trimws(glines)) & !grepl("^\\s*#", glines)]
    if (!length(glines)) stopf("`genes` file contains no records.")
    first_cols <- length(strsplit(glines[1], "\t", fixed = TRUE)[[1]])
    if (first_cols == 4L) {
      tab <- utils::read.table(text = glines, sep = "\t", header = FALSE,
        quote = "", comment.char = "", colClasses = "character", stringsAsFactors = FALSE)
      s <- suppressWarnings(as.numeric(tab[[3]]))
      e <- suppressWarnings(as.numeric(tab[[4]]))
      if (anyNA(c(s, e)) || any(s < 1 | e < s | s != floor(s) | e != floor(e))) {
        stopf("4-column MCScanX gene file requires positive integer `start <= end` in columns 3 and 4.")
      }
      return(data.frame(
        Gene = validate_chr(tab[[2]], "MCScanX gene ID"),
        Chr = validate_chr(tab[[1]], "MCScanX chromosome"),
        Start = pmin(s, e),
        End = pmax(s, e),
        stringsAsFactors = FALSE
      ))
    }
    genes <- read_chr_features(genes)
  }
  if (!is.data.frame(genes)) {
    stopf("`genes` must be a gene coordinate data frame or file path.")
  }
  if (ncol(genes) == 4L && !all(c("Chr", "Start", "End") %in% names(genes))) {
    s <- validate_coordinate(genes[[3]], "start")
    e <- validate_coordinate(genes[[4]], "end")
    if (any(s < 1 | e < s | s != floor(s) | e != floor(e))) {
      stopf("4-column `genes` table requires positive integer `start <= end` in columns 3 and 4.")
    }
    return(data.frame(
      Gene = validate_chr(genes[[2]], "gene"),
      Chr = validate_chr(genes[[1]], "chr"),
      Start = s,
      End = e,
      stringsAsFactors = FALSE
    ))
  }
  if (!all(c("Chr", "Start", "End") %in% names(genes))) {
    stopf("`genes` table must contain `Chr`, `Start`, `End` and a gene identifier column (`Gene`, `ID` or `gene_id`).")
  }
  sub_genes <- genes
  if ("Type" %in% names(sub_genes)) {
    preferred <- sub_genes$Type %in% c("gene", "mRNA", "transcript", "ncRNA_gene", "pseudogene")
    if (any(preferred)) sub_genes <- sub_genes[preferred, , drop = FALSE]
  }
  id_col <- intersect(c("Gene", "ID", "gene_id", "Transcript", "transcript_id", "Name"), names(sub_genes))
  if (!length(id_col)) {
    stopf("`genes` table must contain an identifier column (`Gene`, `ID` or `gene_id`).")
  }
  raw_id <- as.character(sub_genes[[id_col[1]]])
  keep <- !is.na(raw_id) & nzchar(raw_id)
  sub_genes <- sub_genes[keep, , drop = FALSE]
  raw_id <- raw_id[keep]
  base_df <- data.frame(
    Gene = raw_id,
    Chr = validate_chr(sub_genes$Chr, "Chr"),
    Start = validate_coordinate(sub_genes$Start, "Start"),
    End = validate_coordinate(sub_genes$End, "End"),
    stringsAsFactors = FALSE
  )
  stripped <- sub("^(gene|transcript|mRNA):", "", raw_id)
  extra <- stripped != raw_id & !(stripped %in% raw_id)
  if (any(extra)) {
    alias_df <- base_df[extra, , drop = FALSE]
    alias_df$Gene <- stripped[extra]
    base_df <- rbind(base_df, alias_df)
  }
  base_df[!duplicated(base_df$Gene), , drop = FALSE]
}

parse_mcscanx_links <- function(raw_lines, genes, collapse_blocks = FALSE, min_length = 0) {
  gene_lut <- load_mcscanx_genes(genes)
  lines <- raw_lines[nzchar(trimws(raw_lines))]

  is_hdr <- grepl("^##\\s*Alignment\\s+[0-9]+", lines)
  is_pair <- !grepl("^\\s*#", lines) & grepl("^\\s*[0-9]+\\s*-\\s*[0-9]+:", lines)

  empty_mcscanx <- function() {
    if (isTRUE(collapse_blocks)) {
      data.frame(
        Chr1 = character(0), Start1 = numeric(0), End1 = numeric(0),
        Chr2 = character(0), Start2 = numeric(0), End2 = numeric(0),
        Orientation = character(0), Block = integer(0),
        N = integer(0), Score = numeric(0), E_value = numeric(0),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Chr1 = character(0), Start1 = numeric(0), End1 = numeric(0),
        Chr2 = character(0), Start2 = numeric(0), End2 = numeric(0),
        Orientation = character(0), Block = integer(0),
        Gene1 = character(0), Gene2 = character(0), E_value = numeric(0),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!any(is_pair)) return(empty_mcscanx())

  hdr_idx <- which(is_hdr)
  pair_idx <- which(is_pair)
  block_num <- integer(length(lines))
  block_score <- rep(NA_real_, length(lines))
  block_eval <- rep(NA_real_, length(lines))
  block_orient <- rep(NA_character_, length(lines))

  if (length(hdr_idx)) {
    hdr_text <- lines[hdr_idx]
    b_ids <- as.integer(sub("^##\\s*Alignment\\s+([0-9]+).*", "\\1", hdr_text))
    b_sc <- ifelse(grepl("score=([0-9.eE+-]+)", hdr_text),
      suppressWarnings(as.numeric(sub(".*score=([0-9.eE+-]+).*", "\\1", hdr_text))), NA_real_)
    b_ev <- ifelse(grepl("e_value=([0-9.eE+-]+)", hdr_text),
      suppressWarnings(as.numeric(sub(".*e_value=([0-9.eE+-]+).*", "\\1", hdr_text))), NA_real_)
    b_or <- ifelse(grepl("\\bminus\\s*$", trimws(hdr_text)), "-",
      ifelse(grepl("\\bplus\\s*$", trimws(hdr_text)), "+", NA_character_))
    assign_hdr <- findInterval(pair_idx, hdr_idx)
    valid_hdr <- assign_hdr > 0L
    block_num[pair_idx[valid_hdr]] <- b_ids[assign_hdr[valid_hdr]]
    block_score[pair_idx[valid_hdr]] <- b_sc[assign_hdr[valid_hdr]]
    block_eval[pair_idx[valid_hdr]] <- b_ev[assign_hdr[valid_hdr]]
    block_orient[pair_idx[valid_hdr]] <- b_or[assign_hdr[valid_hdr]]
  }

  pair_lines <- lines[pair_idx]
  if (!length(hdr_idx)) {
    block_num[pair_idx] <- as.integer(sub("^\\s*([0-9]+)\\s*-.*", "\\1", pair_lines))
  }
  parts <- strsplit(pair_lines, "\t", fixed = TRUE)
  if (any(lengths(parts) < 3L)) {
    stopf("MCScanX collinearity pair lines must contain tab-separated gene1 and gene2 columns.")
  }
  g1 <- trimws(vapply(parts, `[[`, character(1), 2L))
  g2 <- trimws(vapply(parts, `[[`, character(1), 3L))
  pair_ev <- vapply(parts, function(p) {
    if (length(p) >= 4L) suppressWarnings(as.numeric(trimws(p[4L]))) else NA_real_
  }, numeric(1))

  m1 <- match(g1, gene_lut$Gene)
  m2 <- match(g2, gene_lut$Gene)
  missing_genes <- unique(c(g1[is.na(m1)], g2[is.na(m2)]))
  if (length(missing_genes)) {
    stopf("MCScanX collinearity references unknown gene ID(s) absent from `genes`: %s.",
          format_chr_rows(utils::head(missing_genes, 5L)))
  }

  pairs_df <- data.frame(
    Chr1 = gene_lut$Chr[m1],
    Start1 = gene_lut$Start[m1],
    End1 = gene_lut$End[m1],
    Chr2 = gene_lut$Chr[m2],
    Start2 = gene_lut$Start[m2],
    End2 = gene_lut$End[m2],
    Orientation = ifelse(is.na(block_orient[pair_idx]), "+", block_orient[pair_idx]),
    Block = block_num[pair_idx],
    Gene1 = g1,
    Gene2 = g2,
    E_value = pair_ev,
    .block_score = block_score[pair_idx],
    .block_eval = block_eval[pair_idx],
    .raw_orient = block_orient[pair_idx],
    stringsAsFactors = FALSE
  )

  if (!isTRUE(collapse_blocks)) {
    pairs_df$.block_score <- NULL
    pairs_df$.block_eval <- NULL
    pairs_df$.raw_orient <- NULL
    if (min_length > 0) {
      span <- pmin(pairs_df$End1 - pairs_df$Start1 + 1, pairs_df$End2 - pairs_df$Start2 + 1)
      pairs_df <- pairs_df[span >= min_length, , drop = FALSE]
    }
    rownames(pairs_df) <- NULL
    return(pairs_df)
  }

  grp_key <- paste(pairs_df$Block, pairs_df$Chr1, pairs_df$Chr2, sep = "\r")
  grp_order <- unique(grp_key)
  groups <- split(pairs_df, factor(grp_key, levels = grp_order))
  blocks <- lapply(groups, function(df) {
    orient <- df$.raw_orient[1]
    if (is.na(orient)) {
      if (nrow(df) >= 2L) {
        mid1 <- (df$Start1 + df$End1) / 2
        mid2 <- (df$Start2 + df$End2) / 2
        ord <- order(mid1)
        slope <- sum(diff(mid1[ord]) * diff(mid2[ord]))
        orient <- if (slope < 0) "-" else "+"
      } else {
        orient <- "+"
      }
    }
    data.frame(
      Chr1 = df$Chr1[1],
      Start1 = min(df$Start1),
      End1 = max(df$End1),
      Chr2 = df$Chr2[1],
      Start2 = min(df$Start2),
      End2 = max(df$End2),
      Orientation = orient,
      Block = df$Block[1],
      N = nrow(df),
      Score = df$.block_score[1],
      E_value = df$.block_eval[1],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, blocks)
  span <- pmin(out$End1 - out$Start1 + 1, out$End2 - out$Start2 + 1)
  out <- out[span >= min_length, , drop = FALSE]
  rownames(out) <- NULL
  out
}
