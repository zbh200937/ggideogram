# ggideogram 0.5.0

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
