# ggideogram 0.6.2

- Added `chr_layout = c("grid", "serial")` and `row_align = c("left", "center")`
  to `ggideogram()` and `ideogram_layout()` for end-to-end single-genome tracks
  (e.g., GWAS Manhattan plots) and multi-species macrosynteny stacks, with
  automatic bold genome row titles, uncluttered exterior/inside-arm chromosome
  labels, and adaptive end-cap clamping for short chromosomes/scaffolds.
- Added `scale_length = "per_genome"` (for both linear serial and circular
  layouts) so genomes with orders-of-magnitude size contrast share equal total
  display span while preserving uniform bp scaling within each genome.
- Added `curve_style = c("auto", "straight", "quadratic", "bezier")` to
  `geom_chrlink()`, providing perpendicular-exit cubic S-curve Bezier ribbons
  (enabled automatically for cross-row links in `chr_layout = "serial"`).
- Added `merge_chr_links()` to chain dense co-directional and inverted pairwise
  alignment fragments (PAF/SyRI) into contiguous synteny blocks with union
  aligned base-pair summaries.
- Added explicit `unknown_chr = c("error", "drop")` to `read_chr_links()`,
  `read_chr_features()`, `as_chr_features()`, and `geom_chrlink()` to filter
  unplaced scaffolds cleanly when requested while preserving strict bounds
  checking on retained chromosomes.
- Automatically scaled down default circular sector gaps (`gap_angle = NULL`)
  when chromosome counts exceed 90 (capping total inter-sector gaps at 25% of
  available ring span), eliminating failures on assemblies with 171+ scaffolds.
- Accelerated high-throughput parsing, binning and rendering: Cartesian fast
  draw for pre-sampled circular `geom_chrlink()` layers (10–130x faster grob/PNG
  rendering with ~10x lower memory on 10k links), pre-filtered SyRI / C-level
  `scan()` PAF / vectorized MCScanX block collapsing, PCRE-vectorized GFF3/GTF
  attribute parsing (40–70x faster), vectorized `pmin`/`pmax` track rectangle
  projection, `tabulate()` + pre-bucketed `bin_genome()`, vectorized circular
  label leaders, and sweep-line conflict pruning in `geom_text_repel`.

# ggideogram 0.6.1

- Decoupled linear chromosome primary names and track titles from cross-track
  start-margin push: chromosome names stay snug beside their chromosome bodies
  (retreating outward only on actual transverse bounding-box collision), and
  each track title clears only its own start-side axis or channel labels
  without stacking non-overlapping widths in `plot.margin`.
- Extended `geom_track(label = ...)` with `angle` (degrees; automatically `90`
  for `ylab`-style vertical titles on wide horizontal discrete-`y` tracks with
  start-side category labels) and `position = c("side", "top")` (left-aligned
  top-corner header in horizontal layouts with zero left-margin overhang).
- Added `pack = NULL` (`TRUE`/`FALSE`/`NULL`) and `label_position = c("auto",
  "above", "beside", "start")` to `geom_genemodel()`, packing non-overlapping
  genes (accounting for both exon intervals and label widths) into compact
  lanes with local labels directly above or beside each gene in linear
  `mode = "gene"` tracks while preserving `"start"` channel alignment for
  `mode = "transcript"`, explicit `lane_order` and circular layouts.

# ggideogram 0.6.0

- Removed redundant legacy wrappers `ideogram()`, `GFFex()`, `save_ideogram()`
  and `as_ideogram_grob()` in favor of `ggideogram()`, `read_chr_features()` +
  `bin_genome()`, `ggplot2::ggsave()` and `ggplot2::ggplotGrob()`.
- Added `resolve_gene_models()` to expand GFF3 (`ID`/`Parent`, including
  multi-parent shared exons, non-coding transcripts and two-level `gene ->
  exon/CDS` structures) and GTF (`gene_id`/`transcript_id`) tables directly for
  `geom_genemodel()`.
- Added `read_chr_links()` to parse `minimap2` PAF (`.paf`/`.paf.gz`), `SyRI`
  structural rearrangement outputs (`.out`) and `MCScanX` `.collinearity`
  alignments (gene pairs or `collapse_blocks = TRUE` block ribbons) into
  1-based closed coordinate tables with optional `chr_key()` genome keys.
- Extended `read_chr_features()` with `format = "bedgraph"` (`.bedgraph`/`.bg`).
- Extended `geom_track()` with native discrete/categorical `y` channels
  (`factor` or `character`, automatic `1..K` sub-lanes and category axis
  labels), per-chromosome `geom_smooth()` (`StatSmooth`), track-edge
  `geom_rug()`, `geom_hline()`, `geom_vline()`, interval-only `geom_rect()` /
  `geom_linerange()` / `geom_errorbar()`, and `log10`-safe `geom_col()` /
  `geom_area()` baselines.
- Extended `axis_units` and `geom_chr(component = "axis", units = ...)` with
  `"cM"` for genetic linkage maps and `"none"` for unitless numeric ticks.
- Generalized horizontal microsynteny stacks to four or more genomes (`N >= 2`)
  with shared corridor spacing, diamond middle-genome locus markers and a
  shared bottom distance ruler.
- Eliminated hidden default graphics device side-effects (`Rplots.pdf`) during
  lazy text and label layout construction, fixed outer tick-length clearance in
  axis title reservation, and aligned gene-model track titles and chromosome
  names without over-displacing non-collinear names.
- Vectorized circular/curved `geom_chrlink()` Bezier and arc sampling,
  `geom_genemodel()` interval merging and $O(1)$ lane lookup, and deduplicated
  lazy text-dimension measurements.
- Added the end-to-end `ggideogram-workflows` vignette and multi-domain
  pangenome, QTL fine-mapping, T2T SyRI, breeding BinMap and Marey map gallery
  scripts.

# ggideogram 0.5.0

- `min_gap` sets a minimum clearance for beside tracks. The bp-axis
  `axis_gap` and the track value-axis `axis = list(gap = ...)` control
  independent distances.
- `axis_angle` rotates bp tick labels. `axis_check_overlap` controls label
  suppression at draw time, including for explicit breaks. Tick marks remain
  visible, and rotated label bounds are used for collision checks.
- Connections automatically draw translucent bands for intervals and lines for
  point endpoints. Explicit connection types still take precedence.
- Forward and reversed bands between chromosomes use the same default edges.
  Linear layouts use straight edges and circular layouts use curved edges.
- Track value-axis text adapts to the rendered axis length and label spacing.
  An explicit `axis$size` fixes its physical size.
  Track titles at shared starts in linear layouts align and clear nearby axes.
  Chromosome names join the same alignment when the label boxes fit.

- Automatic track width follows the final side and contents, and `width = NULL`
  or `gap = NULL` restores automatic sizing. Track edits retain explicit
  dimensions; chromosome selections are checked against the source karyotype.
- Automatic value axes skip chromosomes without a resolved range. Blank
  axis labels hide text without removing ticks, including in circular gaps.
- Point and interval annotations in a track share its source and styling
  without inheriting unrelated value mappings. Equivalent starting angles
  produce the same circular layout and label positions.
- Ordered labels minimize displacement along the chromosome or sector when
  there is room. Linear synteny bands expand the source table once and sample
  only curved bands.
- Separate overlay tracks cannot overlap on the same chromosome, and untracked
  markers have their own space before beside tracks. Local synteny reuses space
  for tracks assigned to different chromosomes and keeps automatic ribbons
  between adjacent chromosomes out of their tracks.

- Insets check physical and mixed-unit viewports in the final panel, including
  collisions across layers and with neighboring tracks. Reusing a grob at a
  different size recalculates the viewport bounds.
- Empty or out-of-view text annotations no longer fail or allocate an unused
  label track. Spread and repel labels share source validation and selection;
  label columns and inset checks reuse the same track boundaries.
- Centered track titles, tick labels and chromosome names reserve transverse
  text overhang even when their anchors lie just inside the panel boundary.

- Track titles clear value-axis labels at the same endpoint. Removing or
  shortening annotations releases their automatic margins and heading space,
  while explicitly supplied theme spacing is retained.
- Untracked locus markers reserve space between neighboring chromosomes.
  Changing their gap rebuilds dependent components through the shared layout.
- Physical text spacing is maintained in one module; adding components no
  longer changes the name-spacing parameters of an existing plot.

- Chromosome axes retain default and named labels when ggalign reorders or
  groups observations. Complete ggalign layouts, compositions and nested
  insets enter `geom_locus_inset()` as standard grobs via `ggalignGrob()`.
- `save_ideogram()` delegates object support to `ggplot2::ggsave()`, including
  grobs and compositions with a standard drawing method.
- Added optional ggalign interoperability examples and regression coverage.

- Native tile footprints and nudge/stack/fill results train automatic track
  ranges. Completely missing layers can use explicit or shared finite limits.
- Chromosome-axis labels preserve plotmath expressions. Gene-model labels
  expose independent colour, family and font face, including measured italic
  text clearance after track updates.
- Circular coordinate distances retain missing path separators produced by
  native area layers on ggplot2 3.5.0.
- Track rebuilding is isolated from public geometry constructors. Compatibility
  checks cover minimum and current ggplot2 versions with matching extensions.
- Repelled labels require ggrepel 0.9.6 or newer to retain valid circular labels
  through the native coordinate transformation.
- Gallery figures use compact typography and separate captions; the composition
  example compares gene and LTR density per Mb using actual window widths.
  Local gene models, counts and exon coverage share the protein-coding scope.

- Cytoband tables consistently use 1-based closed intervals; UCSC imports,
  centromere boundaries and single-base bands share the same bp projection.
- Locus layers reject unsupported statistics and honor explicit native style
  and position overrides. Chromosome bodies accept alpha and linetype;
  unsupported component parameters are reported.
- GFF3 examples expand multiple Parent relationships. Local coverage is drawn
  by window boundaries, and gallery figures include source context.
- Coverage and weighted-mean summaries scan sorted overlapping intervals.
  Track ranges, distribution statistics and value axes have separate modules.

- Window statistics share a metadata contract across FUN, explicit methods
  and GFFex. count_position selects midpoint, start or end assignments;
  GFFex retains start counts through the same engine. Summary units are not
  inferred from value-column names. Metadata reads only the exact na.rm
  argument, preserving lazy evaluation of other arguments passed to FUN.
- Package and gallery defaults share a 6 mm legend gap. Circular component
  updates preserve explicitly chosen legend gaps.

- Ordinary beside tracks now default to 3 chromosome body widths. Gene-model
  and overlay defaults remain 6 and 1; explicit widths take priority.

- Identity track layers retain missing observations for native path breaks.
  Occupied windows preserve missing summaries independently of empty-window fill.
- GFFex counts feature starts with the shared window engine, respecting nonzero
  source bounds, complete interval validation and mapped chromosome keys.
  Interval insets validate closed integer coordinates at start - 1 and end.
- Gene-density examples normalize actual window widths; RNA examples identify
  the original random annotation sample. Multi-genome REST snapshots reproduce
  selected alignment records and verify regenerated table checksums.
- Gallery track-axis sizes and legend spacing remain independently adjustable;
  grid planning uses the same constructor parameters as the final plot.

- Local bp axes preserve fractional unit offsets. Vertical bp axes receive
  wider default chromosome spacing; explicit spacing remains available.
- Gene-model tracks default to a wider layout, with measured label margins
  and chromosome-name clearance. Track rebuilds retain value-axis margins.
- Track titles and value axes inherit base_family unless explicitly styled.
  Component colours, line widths and font sizes share one set of defaults.

- Native boxplot and violin statistics and positions run before chromosome
  projection, preserving bp widths, value bandwidths and nonlinear summaries.
  Automatic limits include untrimmed density tails and enabled notches.
- All bin_genome methods validate complete source intervals and begin windows
  at the karyotype's declared Start boundary.
- Track rebuilds retain interleaved draw order, chromosome metadata, locus
  attachments and effective native aesthetic prototypes across multiple scales.
- Native column footprints stop at visible chromosome boundaries while source
  observations and widths remain unchanged.
- Standard locus Positions retain their displayed adjustments; statistical
  tracks defer after_stat and after_scale mappings to the native pipeline.
- Closed interval annotations retain single-base spans and exact local clipping.
  Gene and transcript grouping uses collision-free chromosome/model keys.
- Circular opening_angle controls the closing opening independently of the
  chromosome gaps. The default opening is 20 degrees; explicit gap_angle values
  retain their closing value unless opening_angle is supplied.
- Circular track titles and values retain their track radii with bounded
  tangential spacing. Crowded openings report insufficient space; numerical
  leaders are short native segments.
- Rice comparison examples use the karyotype's display label for species names
  on the shared outer name ring.
- Object-based public layers consolidate chromosome selection, gene models,
  loci, labels and chromosome relations into a small ggplot-style grammar.
- Locus text uses geom_locus(geom = "text") with position = "identity",
  "spread" or "repel". Text spacing, font and leaders belong to the annotation;
  optional label_width reserves its display space within the host layout.
- Circular free text layouts keep native text boxes within their allocated
  sector when a display translation fits; font size and source anchors persist.
- Shared-coordinate tracks own raw bp/value ranges; positioned complete plots
  own their coordinates, theme, guides and viewport dimensions. Closed-interval
  inset anchors use the same midpoint as source nodes and interval geometry.
- Named geom_track declarations own their data, mapping, native layers, title
  and value axis. Left/right and inner/outer use the same entry point.
- The same track identifier supports geometry, range, data and chromosome
  selection edits. Omitted settings persist; title/axis lists update individual
  settings, and explicit NULL/FALSE removes titles or axes. replace = TRUE
  replaces contained layers while retaining the track's other settings.
- Track titles default to 2.8 mm and accept native font styling and em clearance.
  Automatic axes follow the track's chromosome ownership; overlay tiles and
  native locus layers retain stable Geom inheritance during track preparation.
- Locus annotations accept a relative track_position, so node markers can follow
  a track boundary and stay attached to circular relation ports after relayout.
- Adding a track rebuilds all prior source components together, including gene
  structures, node connections, labels and complete-plot insets.
- Native area and ribbon layers complete raw-coordinate statistics and position
  adjustments before projection. Automatic limits include interpolated stacks.
- chr_nodes resolves point and interval connections through source feature IDs;
  local same-chromosome and reversed interval relations default to smooth arcs.
- Ordered dense labels retain true genomic anchors and use device text dimensions.
  Near text edges share a straight baseline or circular arc, with readable normal
  or radial facing and a 5 mm connection zone. Chromosome names share a name ring.
- Chromosome names own the complete track group and use physical text clearance.
  Default chromosome, feature, bp-axis and value-axis sizes form a common hierarchy.
- Circular layouts reserve a closing gap above and to the right, with a vertical
  left boundary. Shared value axes and track titles occupy this gap; explicit
  angles and chromosome-specific axes remain available. Numerical labels use
  measured text clearance and native segment leaders while ticks retain their
  true radial positions and fonts retain physical sizes.
- Selected chromosome bp axes retain their selection when loci or tracks are added.
- Public projection helpers accept a plot directly. Updated real-data examples
  cover shared local coordinates, multiple tracks, nodes and external composition.

# ggideogram 0.4.1

- Chromosome metadata attaches by source and display keys independently after
  genome ordering or local-window selection.
- Generic native tracks project mapped and constant segment/rectangle endpoints
  from original bp, and include endpoint values in automatic track limits.
- Native rectangle widths retain their drawing geometry across local-window
  boundaries and short final bins in circular layouts.
- Identity geom constructors that fix their Stat internally retain native
  columns and position adjustments with ggplot2 3.5.2 and 4.0.3.
- Text-track adapters project standard positions through the final layout;
  third-party repel labels preserve true anchors and requested nudges without
  being pushed to the panel edge.
- Boxplot and violin position adjustments run independently on each chromosome
  while retaining ordinary within-chromosome group dodging.
- Gene and transcript lane labels anchor to the displayed start of reversed
  chromosomes and remain outside circular structures.
- Circular sector names reserve their complete multiline height; annular label
  columns check the entire text box, including concave boundaries.
- Clipped complete-plot insets mask their full viewport, preventing child guides
  and axes from covering the host chromosome.
- Plot headings reserve physical clearance above names and local axes. Gallery
  PDFs use Cairo to retain Unicode text.

# ggideogram 0.4.0

- Circular ideograms allocate chromosome sectors by displayed bp length, with
  explicit radius, start angle, named gaps, clockwise order and reversed loci.
- Inner/outer tracks share existing value scales and native ggplot2 geoms;
  internal intervals and tiles follow the annular body and centromere outline.
- Circular connections use native segments along quadratic Bezier curves;
  synteny ribbons retain both interval arcs and explicit orientation. Automatic
  endpoints attach inside the innermost declared track.
- Sector names, genomic/value axes, physical label columns, gene structures,
  radial insets, external composition and ggiraph layers use the same layout.
- Added three reproducible circular PNG/PDF examples using human gene windows,
  Arabidopsis MCScanX blocks and all 21 wheat A/B/D chromosomes.

# ggideogram 0.3.0

- Multi-genome karyotypes share explicit genome/assembly/chromosome keys across
  layers and guides, with genome/homolog ordering and reversible display axes.
- Annotation assembly checks use matched chromosomes and mapped assembly fields;
  closed intervals respect nonzero source plotting boundaries.
- Gene/transcript structures share lane spacing and block thickness across
  chromosomes with different visible model counts.
- Comparison examples place numerical tracks and bp axes outside the connection
  region, with corresponding chromosome-specific sides.
- Optional ggiraph widgets support locus tooltips, gene links, category filters,
  selection and original-interval TSV downloads at a declared physical size.
- Added 21 reproducible PNG/PDF examples and three corresponding interactive
  examples, including real poplar haplotypes, wheat subgenomes and rice comparisons.

- Physical label columns adapt ggrepel's device measurements, retain true loci,
  and draw leaders with native segments; side, priority and omission are explicit.
- Explicit window methods provide grouped midpoint counts, union coverage and
  overlap-weighted means, with observed/missing counts and actual window widths.
- Paired loci and interval blocks connect chromosomes through native segments
  and polygons; original endpoints, relative orientation and clipped views are preserved.
- Internal overlay tracks support explicit offsets and masked interval fills via
  `geom_chr(component = "fill")`, sharing chromosome coordinates with native point/text layers.

- Add gene/transcript structures using native rectangles and segments, with strand arrows and window clipping.
- Place chromosome names at the aligned start in both orientations.

- Added `read_karyotype()` for chromosome sizes and FAI indexes, including
  gzip input, explicit chromosome selection and preserved sequence names.
- Added BED/GFF3/GTF and GRanges annotation adapters with explicit coordinate
  conversion, plus full-length karyotypes from GRanges/Seqinfo.
- Added `chr_view()` and `view_chr_data()` for local genomic views with
  original bp axes, clipped intervals and flat truncated chromosome ends.

# ggideogram 0.2.0

- Fixed chromosome rows to align genomic starts (normally 0 bp) in both
  orientations, keeping markers, tracks and base-pair axes on the same origin.
- Rebuilt chromosome ideograms as composable ggplot2 objects with one
  dimensionless layout and chromosome-aware coordinate system.
- Added native chromosome x/y axis guides for ordinary ggplot2 geoms.
- Added standard ggplot2 markers, declarative bp-aligned tracks, complete-plot
  insets, and patchwork/cowplot interoperability.
- Added physical axis-label clearance and chromosome-silhouette clipping for
  internal heatmap tiles.
- Added a reproducible six-figure gallery based on the original RIdeogram data.
- Documented RIdeogram provenance, attribution, and citation in Chinese and
  English.
