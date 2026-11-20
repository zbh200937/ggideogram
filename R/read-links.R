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
#' @param unknown_chr Policy when link records reference chromosomes absent from
#'   `karyotype` (or genes absent from `genes` in `"mcscanx"`): `"error"`
#'   (default) rejects unknown identifiers; `"drop"` explicitly removes those
#'   records while still validating bounds on retained chromosomes.
#'
#' @return A data frame with `Chr1`, `Start1`, `End1`, `Chr2`, `Start2`, `End2`,
#'   `Orientation` (`"+"` or `"-"`) and format-specific metadata columns.
#' @seealso [merge_chr_links()], [geom_chrlink()], [chr_key()], [read_chr_features()]
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
    karyotype = NULL,
    unknown_chr = c("error", "drop")) {
  format <- match.arg(format)
  unknown_chr <- match.arg(unknown_chr)
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
      collapse_blocks = collapse_blocks, min_length = min_length,
      drop_unknown_genes = identical(unknown_chr, "drop"))
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
    if (identical(unknown_chr, "drop")) {
      keep <- !is.na(i1) & !is.na(i2)
      out <- out[keep, , drop = FALSE]
      i1 <- i1[keep]; i2 <- i2[keep]
      rownames(out) <- NULL
    } else {
      unknown <- unique(c(out$Chr1[is.na(i1)], out$Chr2[is.na(i2)]))
      if (length(unknown)) {
        stopf("Unknown link chromosomes: %s.", format_chr_rows(unknown))
      }
    }
    if (nrow(out) && any(out$Start1 - 1 < k$.start[i1] | out$End1 > k$.end[i1] |
            out$Start2 - 1 < k$.start[i2] | out$End2 > k$.end[i2])) {
      stopf("Link intervals are outside their karyotype chromosome bounds.")
    }
  }

  attr(out, "source_format") <- format
  attr(out, "coordinates") <- "1-based-closed"
  out
}

escape_regex_literal <- function(x) {
  gsub("([][{}()+*^$.|\\\\?])", "\\\\\\1", x, perl = TRUE)
}

detect_link_format <- function(file, raw_lines) {
  clean <- sub("\\.gz$", "", file, ignore.case = TRUE)
  ext <- tolower(tools::file_ext(clean))
  if (ext == "paf") return("paf")
  if (ext %in% c("collinearity", "collinear")) return("mcscanx")
  if (ext == "syri") return("syri")

  head_lines <- utils::head(raw_lines, 200L)
  non_empty <- head_lines[nzchar(head_lines) & !grepl("^\\s*$", head_lines, perl = TRUE)]
  if (any(grepl("^##\\s*Alignment\\b", non_empty, perl = TRUE))) return("mcscanx")
  body <- non_empty[!grepl("^\\s*#", non_empty, perl = TRUE)]
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

strip_comment_lines <- function(raw_lines) {
  lines <- raw_lines[nzchar(raw_lines) & !startsWith(raw_lines, "#")]
  if (!length(lines)) return(lines)
  has_lead_ws <- startsWith(lines, " ") | startsWith(lines, "\t")
  if (any(has_lead_ws)) {
    drop <- has_lead_ws & grepl("^\\s*(?:#|$)", lines, perl = TRUE)
    if (any(drop)) lines <- lines[!drop]
  }
  lines
}

parse_paf_links <- function(raw_lines, min_length = 0, min_mapq = 0) {
  lines <- strip_comment_lines(raw_lines)
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

  # Use C-level scan(flush = TRUE) to read only the 12 core columns and discard
  # arbitrarily long trailing SAM/CIGAR tags (e.g. cg:Z:) without copying them.
  cols <- scan(
    text = lines,
    what = rep(list(character(0)), 12L),
    sep = "\t", quote = "", comment.char = "",
    fill = TRUE, flush = TRUE, na.strings = character(0), quiet = TRUE
  )
  short <- !nzchar(cols[[12L]])
  if (any(short)) {
    bad_i <- which(short)[1L]
    n_col <- length(strsplit(lines[bad_i], "\t", fixed = TRUE)[[1L]])
    if (n_col < 12L) {
      stopf("PAF input requires at least 12 tab-separated columns (row %d has %d).", bad_i, n_col)
    }
  }
  rm(lines)

  qname <- validate_chr(cols[[1L]], "PAF query name")
  tname <- validate_chr(cols[[6L]], "PAF target name")
  strand <- cols[[5L]]
  if (any(!strand %in% c("+", "-"))) {
    stopf("PAF column 5 (`strand`) must be `+` or `-`.")
  }

  num_idx <- c(2L, 3L, 4L, 7L, 8L, 9L, 10L, 11L, 12L)
  nums <- lapply(num_idx, function(j) suppressWarnings(as.numeric(cols[[j]])))
  if (any(vapply(nums, function(v) anyNA(v) || any(!is.finite(v)) || any(v != floor(v)) || any(v < 0), logical(1)))) {
    stopf("PAF coordinate, length, match and MapQ columns must be non-negative integers.")
  }
  qlen <- nums[[1L]]; qstart <- nums[[2L]]; qend <- nums[[3L]]
  tlen <- nums[[4L]]; tstart <- nums[[5L]]; tend <- nums[[6L]]
  nmatch <- nums[[7L]]; alen <- nums[[8L]]; mapq <- nums[[9L]]

  if (any(qstart >= qend | tstart >= tend)) {
    stopf("PAF alignments require `start < end` on both target and query.")
  }
  if (any(qend > qlen | tend > tlen)) {
    stopf("PAF alignment coordinates exceed reported sequence lengths.")
  }

  keep <- alen >= min_length & mapq >= min_mapq
  if (!any(keep)) return(empty_paf())
  out <- data.frame(
    Chr1 = tname[keep],
    Start1 = tstart[keep] + 1,
    End1 = tend[keep],
    Chr2 = qname[keep],
    Start2 = qstart[keep] + 1,
    End2 = qend[keep],
    Orientation = strand[keep],
    Matches = nmatch[keep],
    Block_length = alen[keep],
    Identity = ifelse(alen[keep] > 0, nmatch[keep] / alen[keep], NA_real_),
    MapQ = mapq[keep],
    Target_length = tlen[keep],
    Query_length = qlen[keep],
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

parse_syri_links <- function(raw_lines, types, min_length = 0) {
  if (!is.character(types) || !length(types) || anyNA(types)) {
    stopf("`types` must be a non-empty character vector of SyRI annotation types.")
  }
  lines <- strip_comment_lines(raw_lines)
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

  first_fields <- strsplit(lines[1L], "\t", fixed = TRUE)[[1L]]
  if (length(first_fields) != 12L) {
    stopf("SyRI `.out` input requires 12 tab-separated columns (row 1 has %d).", length(first_fields))
  }
  if (tolower(first_fields[1L]) %in% c("chra", "chr1", "refchr") &&
      is.na(suppressWarnings(as.numeric(first_fields[2L])))) {
    lines <- lines[-1L]
  }
  if (!length(lines)) return(empty_syri())

  # For moderate files verify every line's 12-column schema; for massive files
  # (>50k rows dominated by SNPs/INDELs) pre-filter by target `types` first and
  # validate every retained structural line cleanly.
  checked_all <- length(lines) <= 50000L
  if (checked_all) {
    valid12 <- grepl("^(?:[^\t]*\t){11}[^\t]*$", lines, perl = TRUE)
    if (!all(valid12)) {
      bad_i <- which(!valid12)[1L]
      n_col <- length(strsplit(lines[bad_i], "\t", fixed = TRUE)[[1L]])
      stopf("SyRI `.out` input requires 12 tab-separated columns (row %d has %d).", bad_i, n_col)
    }
    rm(valid12)
  }

  type_pat <- paste0("\t(?:", paste(escape_regex_literal(types), collapse = "|"), ")\t[^\t]*$")
  lines <- lines[grepl(type_pat, lines, perl = TRUE)]
  if (!length(lines)) return(empty_syri())
  if (!checked_all) {
    valid12 <- grepl("^(?:[^\t]*\t){11}[^\t]*$", lines, perl = TRUE)
    if (!all(valid12)) {
      bad_i <- which(!valid12)[1L]
      n_col <- length(strsplit(lines[bad_i], "\t", fixed = TRUE)[[1L]])
      stopf("SyRI `.out` input requires 12 tab-separated columns (row %d has %d).", bad_i, n_col)
    }
    rm(valid12)
  }

  cols <- scan(
    text = lines,
    what = rep(list(character(0)), 12L),
    sep = "\t", quote = "", comment.char = "",
    na.strings = character(0), quiet = TRUE
  )
  rm(lines)

  keep_type <- cols[[11L]] %in% types &
    cols[[1L]] != "-" & cols[[6L]] != "-" &
    cols[[2L]] != "-" & cols[[3L]] != "-" &
    cols[[7L]] != "-" & cols[[8L]] != "-"
  if (!any(keep_type)) return(empty_syri())
  if (!all(keep_type)) cols <- lapply(cols, `[`, keep_type)

  s1 <- suppressWarnings(as.numeric(cols[[2L]]))
  e1 <- suppressWarnings(as.numeric(cols[[3L]]))
  s2 <- suppressWarnings(as.numeric(cols[[7L]]))
  e2 <- suppressWarnings(as.numeric(cols[[8L]]))
  if (anyNA(c(s1, e1, s2, e2)) || any(!is.finite(c(s1, e1, s2, e2))) ||
      any(c(s1, e1, s2, e2) < 1) ||
      any(c(s1, e1, s2, e2) != floor(c(s1, e1, s2, e2)))) {
    stopf("SyRI coordinates for retained structural blocks must be positive integers.")
  }

  inv_types <- c("INV", "INVTR", "INVDP", "INVAL", "INVTRAL", "INVDPAL")
  inverted <- (cols[[11L]] %in% inv_types) | ((s1 > e1) != (s2 > e2))
  start1 <- pmin(s1, e1); end1 <- pmax(s1, e1)
  start2 <- pmin(s2, e2); end2 <- pmax(s2, e2)

  span <- pmin(end1 - start1 + 1, end2 - start2 + 1)
  keep_len <- span >= min_length
  if (!any(keep_len)) return(empty_syri())

  out <- data.frame(
    Chr1 = validate_chr(cols[[1L]][keep_len], "SyRI reference chromosome"),
    Start1 = start1[keep_len],
    End1 = end1[keep_len],
    Chr2 = validate_chr(cols[[6L]][keep_len], "SyRI query chromosome"),
    Start2 = start2[keep_len],
    End2 = end2[keep_len],
    Orientation = ifelse(inverted[keep_len], "-", "+"),
    Type = cols[[11L]][keep_len],
    ID = cols[[9L]][keep_len],
    Parent = cols[[10L]][keep_len],
    Copy_status = cols[[12L]][keep_len],
    stringsAsFactors = FALSE
  )
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
    glines <- strip_comment_lines(readLines(con, warn = FALSE))
    if (!length(glines)) stopf("`genes` file contains no records.")
    first_cols <- length(strsplit(glines[1], "\t", fixed = TRUE)[[1]])
    if (first_cols == 4L) {
      tab <- scan(
        text = glines,
        what = list(character(0), character(0), character(0), character(0)),
        sep = "\t", quote = "", comment.char = "",
        na.strings = character(0), quiet = TRUE
      )
      s <- suppressWarnings(as.numeric(tab[[3L]]))
      e <- suppressWarnings(as.numeric(tab[[4L]]))
      if (anyNA(c(s, e)) || any(s < 1 | e < s | s != floor(s) | e != floor(e))) {
        stopf("4-column MCScanX gene file requires positive integer `start <= end` in columns 3 and 4.")
      }
      return(data.frame(
        Gene = validate_chr(tab[[2L]], "MCScanX gene ID"),
        Chr = validate_chr(tab[[1L]], "MCScanX chromosome"),
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

parse_mcscanx_links <- function(raw_lines, genes, collapse_blocks = FALSE,
                                min_length = 0, drop_unknown_genes = FALSE) {
  gene_lut <- load_mcscanx_genes(genes)
  lines <- raw_lines[nzchar(raw_lines) & !grepl("^\\s*$", raw_lines, perl = TRUE)]

  is_hdr <- grepl("^##\\s*Alignment\\s+[0-9]+", lines, perl = TRUE)
  is_pair <- !grepl("^\\s*#", lines, perl = TRUE) & grepl("^\\s*[0-9]+\\s*-\\s*[0-9]+:", lines, perl = TRUE)

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
  n_pairs <- length(pair_idx)
  b_num <- integer(n_pairs)
  b_score_v <- rep(NA_real_, n_pairs)
  b_eval_v <- rep(NA_real_, n_pairs)
  b_orient_v <- rep(NA_character_, n_pairs)

  if (length(hdr_idx)) {
    hdr_text <- lines[hdr_idx]
    b_ids <- as.integer(sub("^##\\s*Alignment\\s+([0-9]+).*", "\\1", hdr_text, perl = TRUE))
    b_sc <- ifelse(grepl("score=([0-9.eE+-]+)", hdr_text, perl = TRUE),
      suppressWarnings(as.numeric(sub(".*score=([0-9.eE+-]+).*", "\\1", hdr_text, perl = TRUE))), NA_real_)
    b_ev <- ifelse(grepl("e_value=([0-9.eE+-]+)", hdr_text, perl = TRUE),
      suppressWarnings(as.numeric(sub(".*e_value=([0-9.eE+-]+).*", "\\1", hdr_text, perl = TRUE))), NA_real_)
    b_or <- ifelse(grepl("\\bminus\\s*$", hdr_text, perl = TRUE), "-",
      ifelse(grepl("\\bplus\\s*$", hdr_text, perl = TRUE), "+", NA_character_))
    assign_hdr <- findInterval(pair_idx, hdr_idx)
    valid_hdr <- assign_hdr > 0L
    b_num[valid_hdr] <- b_ids[assign_hdr[valid_hdr]]
    b_score_v[valid_hdr] <- b_sc[assign_hdr[valid_hdr]]
    b_eval_v[valid_hdr] <- b_ev[assign_hdr[valid_hdr]]
    b_orient_v[valid_hdr] <- b_or[assign_hdr[valid_hdr]]
  }

  pair_lines <- lines[pair_idx]
  if (!length(hdr_idx)) {
    b_num <- as.integer(sub("^\\s*([0-9]+)\\s*-.*", "\\1", pair_lines, perl = TRUE))
  }
  if (!all(grepl("^[^\t]*\t[^\t]*\t[^\t]*", pair_lines, perl = TRUE))) {
    stopf("MCScanX collinearity pair lines must contain tab-separated gene1 and gene2 columns.")
  }
  pcols <- scan(
    text = pair_lines,
    what = list(character(0), character(0), character(0), character(0)),
    sep = "\t", quote = "", comment.char = "",
    fill = TRUE, flush = TRUE, na.strings = character(0), quiet = TRUE
  )
  g1 <- trimws(pcols[[2L]])
  g2 <- trimws(pcols[[3L]])
  ev_raw <- trimws(pcols[[4L]])
  pair_ev <- ifelse(nzchar(ev_raw), suppressWarnings(as.numeric(ev_raw)), NA_real_)

  m1 <- match(g1, gene_lut$Gene)
  m2 <- match(g2, gene_lut$Gene)
  if (isTRUE(drop_unknown_genes)) {
    keep_g <- !is.na(m1) & !is.na(m2)
    if (!any(keep_g)) return(empty_mcscanx())
    m1 <- m1[keep_g]; m2 <- m2[keep_g]
    g1 <- g1[keep_g]; g2 <- g2[keep_g]
    pair_ev <- pair_ev[keep_g]
    b_num <- b_num[keep_g]; b_score_v <- b_score_v[keep_g]
    b_eval_v <- b_eval_v[keep_g]; b_orient_v <- b_orient_v[keep_g]
  } else {
    missing_genes <- unique(c(g1[is.na(m1)], g2[is.na(m2)]))
    if (length(missing_genes)) {
      stopf("MCScanX collinearity references unknown gene ID(s) absent from `genes`: %s.",
            format_chr_rows(utils::head(missing_genes, 5L)))
    }
  }

  chr1 <- gene_lut$Chr[m1]; start1 <- gene_lut$Start[m1]; end1 <- gene_lut$End[m1]
  chr2 <- gene_lut$Chr[m2]; start2 <- gene_lut$Start[m2]; end2 <- gene_lut$End[m2]

  if (!isTRUE(collapse_blocks)) {
    keep_len <- if (min_length > 0) pmin(end1 - start1 + 1, end2 - start2 + 1) >= min_length else rep(TRUE, length(m1))
    if (!any(keep_len)) return(empty_mcscanx())
    out <- data.frame(
      Chr1 = chr1[keep_len],
      Start1 = start1[keep_len],
      End1 = end1[keep_len],
      Chr2 = chr2[keep_len],
      Start2 = start2[keep_len],
      End2 = end2[keep_len],
      Orientation = ifelse(is.na(b_orient_v[keep_len]), "+", b_orient_v[keep_len]),
      Block = b_num[keep_len],
      Gene1 = g1[keep_len],
      Gene2 = g2[keep_len],
      E_value = pair_ev[keep_len],
      stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    return(out)
  }

  # Vectorized block aggregation (no split/lapply/rbind over data.frame rows).
  grp_key <- paste(b_num, chr1, chr2, sep = "\r")
  grp_levels <- unique(grp_key)
  grp_id <- match(grp_key, grp_levels)
  n_grp <- length(grp_levels)

  ord <- order(grp_id, (start1 + end1) / 2)
  g_ord <- grp_id[ord]
  first_pos <- which(c(TRUE, g_ord[-1L] != g_ord[-length(g_ord)]))
  last_pos <- c(first_pos[-1L] - 1L, length(g_ord))
  n_per <- last_pos - first_pos + 1L
  rep_idx <- ord[first_pos]

  s1_min <- as.numeric(tapply(start1, grp_id, min))
  e1_max <- as.numeric(tapply(end1, grp_id, max))
  s2_min <- as.numeric(tapply(start2, grp_id, min))
  e2_max <- as.numeric(tapply(end2, grp_id, max))

  orient <- b_orient_v[rep_idx]
  miss_or <- which(is.na(orient))
  if (length(miss_or)) {
    mid1_ord <- (start1[ord] + end1[ord]) / 2
    mid2_ord <- (start2[ord] + end2[ord]) / 2
    for (gi in miss_or) {
      if (n_per[gi] >= 2L) {
        rng <- first_pos[gi]:last_pos[gi]
        slope <- sum(diff(mid1_ord[rng]) * diff(mid2_ord[rng]))
        orient[gi] <- if (slope < 0) "-" else "+"
      } else {
        orient[gi] <- "+"
      }
    }
  }

  span <- pmin(e1_max - s1_min + 1, e2_max - s2_min + 1)
  keep_b <- span >= min_length
  if (!any(keep_b)) return(empty_mcscanx())

  out <- data.frame(
    Chr1 = chr1[rep_idx[keep_b]],
    Start1 = s1_min[keep_b],
    End1 = e1_max[keep_b],
    Chr2 = chr2[rep_idx[keep_b]],
    Start2 = s2_min[keep_b],
    End2 = e2_max[keep_b],
    Orientation = orient[keep_b],
    Block = b_num[rep_idx[keep_b]],
    N = as.integer(n_per[keep_b]),
    Score = b_score_v[rep_idx[keep_b]],
    E_value = b_eval_v[rep_idx[keep_b]],
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

#' Merge adjacent collinear alignment fragments into synteny blocks
#'
#' Chains dense pairwise alignment fragments (such as raw `minimap2` PAF or
#' fragmented `SyRI` blocks) into clean, contiguous one-based closed synteny
#' blocks ready for [geom_chrlink()]. Adjacent fragments on the same
#' `(Chr1, Chr2, Orientation)` (and optional `group` columns) are merged when
#' both endpoints progress consistently with the block orientation (`"+"`
#' co-directional, `"-"` inverted) and each end's gap does not exceed
#' `max_gap` base pairs.
#'
#' @param data A paired interval data frame such as the result of
#'   [read_chr_links()].
#' @param max_gap Non-negative maximum gap in base pairs allowed between
#'   consecutive fragments on both target and query ends (default `1e5`).
#' @param min_length Minimum merged block span (minimum of target and query
#'   spans) in base pairs (default `0`).
#' @param min_blocks Minimum number of merged input fragments required to keep
#'   a block (default `1L`).
#' @param group Optional character vector of extra grouping column names
#'   (for example `"Type"` in SyRI tables) that must match to merge fragments.
#' @param chr1,start1,end1,chr2,start2,end2,orientation Column names in `data`.
#'
#' @return A data frame with `Chr1`, `Start1`, `End1`, `Chr2`, `Start2`, `End2`,
#'   `Orientation`, any requested `group` columns, plus `N` (merged fragment
#'   count), `Span1`, `Span2`, `Aligned_bp1`, and `Aligned_bp2`.
#' @seealso [read_chr_links()], [geom_chrlink()]
#' @export
merge_chr_links <- function(
    data,
    max_gap = 1e5,
    min_length = 0,
    min_blocks = 1L,
    group = NULL,
    chr1 = "Chr1",
    start1 = "Start1",
    end1 = "End1",
    chr2 = "Chr2",
    start2 = "Start2",
    end2 = "End2",
    orientation = "Orientation") {
  if (!is.data.frame(data)) stopf("`data` must be a paired-interval data frame.")
  check_nonnegative_layout(max_gap, "max_gap")
  check_nonnegative_layout(min_length, "min_length")
  if (!is.numeric(min_blocks) || length(min_blocks) != 1L || is.na(min_blocks) ||
      !is.finite(min_blocks) || min_blocks < 1 || min_blocks != floor(min_blocks)) {
    stopf("`min_blocks` must be one positive integer.")
  }
  min_blocks <- as.integer(min_blocks)

  req <- c(chr1, start1, end1, chr2, start2, end2)
  miss <- setdiff(req, names(data))
  if (length(miss)) {
    stopf("`data` is missing required column(s): %s.", paste(miss, collapse = ", "))
  }
  if (!is.null(group)) {
    if (!is.character(group) || anyNA(group) || any(!group %in% names(data)) || anyDuplicated(group)) {
      stopf("`group` must name distinct existing columns in `data`.")
    }
  }

  empty_merged <- function() {
    res <- data.frame(
      Chr1 = character(0), Start1 = numeric(0), End1 = numeric(0),
      Chr2 = character(0), Start2 = numeric(0), End2 = numeric(0),
      Orientation = character(0), stringsAsFactors = FALSE
    )
    for (gcol in group) res[[gcol]] <- data[[gcol]][FALSE]
    res$N <- integer(0)
    res$Span1 <- numeric(0); res$Span2 <- numeric(0)
    res$Aligned_bp1 <- numeric(0); res$Aligned_bp2 <- numeric(0)
    attr(res, "coordinates") <- "1-based-closed"
    res
  }
  if (!nrow(data)) return(empty_merged())

  c1 <- validate_chr(data[[chr1]], chr1)
  c2 <- validate_chr(data[[chr2]], chr2)
  s1 <- validate_coordinate(data[[start1]], start1)
  e1 <- validate_coordinate(data[[end1]], end1)
  s2 <- validate_coordinate(data[[start2]], start2)
  e2 <- validate_coordinate(data[[end2]], end2)
  if (any(s1 < 1 | e1 < s1 | s1 != floor(s1) | e1 != floor(e1) |
          s2 < 1 | e2 < s2 | s2 != floor(s2) | e2 != floor(e2))) {
    stopf("Link intervals require one-based closed integer coordinates with start <= end.")
  }
  or_v <- if (orientation %in% names(data)) as.character(data[[orientation]]) else rep("+", nrow(data))
  if (anyNA(or_v) || any(!or_v %in% c("+", "-"))) {
    stopf("`%s` must contain only `+` or `-`.", orientation)
  }

  key_cols <- c(list(c1, c2, or_v), if (length(group)) lapply(group, function(g) encode_group_token(data[[g]])) else list())
  grp_key <- do.call(paste, c(key_cols, list(sep = "\r")))
  ord <- order(grp_key, s1, e1, ifelse(or_v == "+", s2, -e2))
  n <- length(ord)

  gk <- grp_key[ord]
  c1 <- c1[ord]; c2 <- c2[ord]; or_v <- or_v[ord]
  s1 <- s1[ord]; e1 <- e1[ord]; s2 <- s2[ord]; e2 <- e2[ord]
  is_plus <- or_v == "+"

  # Single-pass monotonic chain resolution tracking exact running block bounds.
  chain_id <- integer(n)
  out_e1 <- numeric(n)
  out_s2 <- numeric(n)
  out_e2 <- numeric(n)
  cid <- 1L
  chain_id[1L] <- 1L
  cur_e1 <- e1[1L]
  cur_s2 <- s2[1L]; cur_e2 <- e2[1L]
  cur_last_s2 <- s2[1L]; cur_last_e2 <- e2[1L]

  if (n > 1L) {
    same_g <- gk[-1L] == gk[-n]
    for (i in 2L:n) {
      mergeable <- FALSE
      if (same_g[i - 1L] && s1[i] <= cur_e1 + max_gap + 1) {
        if (is_plus[i]) {
          mergeable <- s2[i] <= cur_e2 + max_gap + 1 &&
            s2[i] >= cur_last_s2 - max_gap &&
            e2[i] >= cur_s2 - max_gap
        } else {
          mergeable <- e2[i] >= cur_s2 - max_gap - 1 &&
            e2[i] <= cur_last_e2 + max_gap &&
            s2[i] <= cur_e2 + max_gap
        }
      }
      if (mergeable) {
        if (e1[i] > cur_e1) cur_e1 <- e1[i]
        if (s2[i] < cur_s2) cur_s2 <- s2[i]
        if (e2[i] > cur_e2) cur_e2 <- e2[i]
        cur_last_s2 <- s2[i]
        cur_last_e2 <- e2[i]
      } else {
        out_e1[cid] <- cur_e1
        out_s2[cid] <- cur_s2
        out_e2[cid] <- cur_e2
        cid <- cid + 1L
        cur_e1 <- e1[i]
        cur_s2 <- s2[i]; cur_e2 <- e2[i]
        cur_last_s2 <- s2[i]; cur_last_e2 <- e2[i]
      }
      chain_id[i] <- cid
    }
  }
  out_e1[cid] <- cur_e1
  out_s2[cid] <- cur_s2
  out_e2[cid] <- cur_e2
  out_e1 <- out_e1[seq_len(cid)]
  out_s2 <- out_s2[seq_len(cid)]
  out_e2 <- out_e2[seq_len(cid)]

  first_i <- which(c(TRUE, chain_id[-1L] != chain_id[-n]))
  last_i <- c(first_i[-1L] - 1L, n)
  counts <- last_i - first_i + 1L

  out_s1 <- s1[first_i]

  span1 <- out_e1 - out_s1 + 1
  span2 <- out_e2 - out_s2 + 1
  keep <- counts >= min_blocks & pmin(span1, span2) >= min_length
  if (!any(keep)) return(empty_merged())

  kept_chains <- which(keep)
  n_out <- length(kept_chains)
  aln1 <- numeric(n_out)
  aln2 <- numeric(n_out)
  for (k in seq_len(n_out)) {
    ci <- kept_chains[k]
    a <- first_i[ci]; b <- last_i[ci]
    if (a == b) {
      aln1[k] <- e1[a] - s1[a] + 1
      aln2[k] <- e2[a] - s2[a] + 1
    } else {
      rng <- a:b
      aln1[k] <- union_covered_bp(s1[rng], e1[rng])
      aln2[k] <- union_covered_bp(s2[rng], e2[rng])
    }
  }

  rep_first <- first_i[kept_chains]
  res <- data.frame(
    Chr1 = c1[rep_first],
    Start1 = out_s1[kept_chains],
    End1 = out_e1[kept_chains],
    Chr2 = c2[rep_first],
    Start2 = out_s2[kept_chains],
    End2 = out_e2[kept_chains],
    Orientation = or_v[rep_first],
    stringsAsFactors = FALSE
  )
  orig_first <- ord[rep_first]
  for (gcol in group) res[[gcol]] <- data[[gcol]][orig_first]
  res$N <- as.integer(counts[kept_chains])
  res$Span1 <- span1[kept_chains]
  res$Span2 <- span2[kept_chains]
  res$Aligned_bp1 <- aln1
  res$Aligned_bp2 <- aln2
  rownames(res) <- NULL
  attr(res, "coordinates") <- "1-based-closed"
  res
}
