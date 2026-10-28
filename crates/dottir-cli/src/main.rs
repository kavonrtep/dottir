//! `dottir` CLI binary — batch mode (spec §4.4.8).
//!
//! Computes a dotplot from a (query, subject) sequence pair and writes
//! a greyscale PNG + a `.params.toml` sidecar. Mirrors C dotter CLI
//! option names where reasonable.
//!
//! Each positional input is either a FASTA file or a GFF3 file that
//! carries its own sequences after a `##FASTA` directive (see
//! `dottir_io::input`), so an annotated sequence can be supplied as a
//! single self-contained file.

mod align_export;
mod find_peaks;
mod genome_batch;
mod large_batch;
mod periodogram;

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand, ValueEnum};

use dottir_core::genome::SortAxis;
use dottir_core::{
    compute_dotplot, pick_auto_zoom, snap_zoom_to_period_divisor, BlastMode, PlotConfig,
    ScoreMatrix, Strand, Triangle, PIXELMAP_FORMAT_VERSION,
};
use dottir_io::{
    params::{
        sha256_bytes, DottirInfo, HostInfo, InputInfo, KarlinInfo, ParamsSidecar, PlotParamsInfo,
    },
    png_export,
};

#[derive(Parser, Debug)]
#[command(
    name = "dottir",
    version,
    about = "Modern Rust reimplementation of Dotter (Sonnhammer & Durbin 1995)"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Compute a dotplot and write a PNG + params sidecar.
    Batch(Box<BatchArgs>),
    /// Compute a self-comparison periodogram per FASTA record and
    /// write a single TSV. Sums the dotplot's windowed, greyramp-
    /// shaped pixels along each anti-diagonal — peaks correspond to
    /// the visible repeat periodicities in the dotplot view.
    Periodogram(periodogram::PeriodogramArgs),
    /// Classify periodogram or FFT peaks into fundamentals,
    /// harmonics, and (optionally) sub-repeats. Takes the TSV
    /// produced by `dottir periodogram` (or its `--fft` companion).
    FindPeaks(find_peaks::FindPeaksArgs),
}

#[derive(clap::Args, Debug)]
#[command(version)]
struct BatchArgs {
    /// Query sequence, horizontal axis. FASTA, or a GFF3 with an
    /// embedded `##FASTA` section.
    #[arg(value_name = "QUERY")]
    query: PathBuf,
    /// Subject sequence, vertical axis. FASTA or GFF3-with-FASTA.
    /// Omit for a self-comparison (subject = query, mirror enabled).
    #[arg(value_name = "SUBJECT")]
    subject: Option<PathBuf>,
    /// Output PNG path. The params sidecar is written to
    /// `<output>.params.toml` alongside it. Optional with
    /// `--export-alignment-at`.
    #[arg(
        short = 'o',
        long,
        value_name = "PATH",
        required_unless_present = "export_alignment_at"
    )]
    output: Option<PathBuf>,

    #[command(flatten)]
    align_export: align_export::AlignExportArgs,

    /// BLAST mode.
    #[arg(long, value_enum, default_value_t = ModeArg::Blastn)]
    mode: ModeArg,

    /// Built-in score matrix name. Default: BLOSUM62 for protein modes,
    /// DNA+5/-4 for BLASTN.
    #[arg(long)]
    matrix: Option<String>,

    /// Window size W. Default: Karlin/Altschul estimate.
    #[arg(short = 'W', long, value_name = "INT")]
    window: Option<u32>,

    /// Zoom factor (pixels per matrix block).
    #[arg(short = 'z', long, default_value_t = 1)]
    zoom: u32,

    /// Pixel factor (multiplier in `min(255, score * pixel_fac / W)`).
    /// `0` = auto-derive from Karlin's expected residue score (dotter's
    /// default behaviour). The resolved value is reported in the
    /// params sidecar.
    #[arg(long, default_value_t = 0)]
    pixel_fac: u32,

    /// Strand selection. BLASTP ignores this.
    #[arg(long, value_enum, default_value_t = StrandArg::Both)]
    strand: StrandArg,

    /// Compute only the forward (watson) strand. Equivalent to
    /// `--strand forward`. Original Dotter flag: `--watson-only`.
    #[arg(long, default_value_t = false, conflicts_with_all = ["crick_only"])]
    watson_only: bool,

    /// Compute only the reverse (crick) strand. Equivalent to
    /// `--strand reverse`. Original Dotter flag: `--crick-only`.
    #[arg(long, default_value_t = false)]
    crick_only: bool,

    /// Pre-process: reverse-complement the query before computation.
    /// Original Dotter flag: `-r`.
    #[arg(short = 'r', long, default_value_t = false)]
    reverse_query: bool,

    /// Pre-process: reverse-complement the subject before computation.
    /// Original Dotter flag: `-v`.
    #[arg(short = 'v', long, default_value_t = false)]
    reverse_subject: bool,

    /// Self-comparison triangle.
    #[arg(long, value_enum, default_value_t = TriangleArg::Both)]
    triangle: TriangleArg,

    /// Disable the self-comparison mirror post-process.
    #[arg(long, default_value_t = false)]
    disable_mirror: bool,

    /// Memory limit for the pixelmap, in bytes.
    #[arg(long, default_value_t = 512 * 1024 * 1024)]
    memory_limit_bytes: u64,

    /// If set, automatically pick `zoom` so the largest output dimension
    /// is at most this. Avoids surprise OOMs on large inputs
    /// (spec §4.4.8). Multi-record repeat arrays may snap to a nearby
    /// divisor zoom so identical records tile uniformly.
    #[arg(long, value_name = "PIXELS")]
    auto_zoom: Option<u32>,

    /// Render PNG/SVG from a boundary-aware normalized raster. This computes
    /// a zoom=1 source pixelmap, then downsamples each FASTA record tile in
    /// tile-local coordinates. The `.dot` export remains the authoritative
    /// Dotter-style pixelmap at `--zoom` / `--auto-zoom`.
    #[arg(long, default_value_t = false)]
    repeat_uniform_rendering: bool,

    /// Skip writing the `.params.toml` sidecar.
    #[arg(long, default_value_t = false)]
    no_sidecar: bool,

    /// Exact PNG width in pixels. Height is derived to preserve
    /// the input aspect ratio. The pixelmap is nearest-neighbour
    /// resized to fit (upscale gives pixel-perfect blocks;
    /// downscale drops pixels — prefer to omit this flag if you
    /// want a smaller image). Capped at 16384 px. Pass 0 to skip
    /// resizing and write the pixelmap at its native resolution.
    /// Has no effect on the SVG (SVG renderers scale themselves).
    #[arg(long, value_name = "PX", default_value_t = 2000)]
    width: u32,

    /// Also write an SVG version of the plot (with embedded PNG +
    /// axis labels) to this path.
    #[arg(long, value_name = "PATH")]
    svg: Option<PathBuf>,

    /// Also write a one-page PDF of the plot (the PNG image with axes, at
    /// the pixelmap's resolution).
    #[arg(long, value_name = "PATH")]
    pdf: Option<PathBuf>,

    /// Draw the alignments of a PAF file over the plot (HSPs, spec
    /// §4.2.3). Records not plotted are skipped; either orientation of
    /// the PAF works. The PNG and SVG become RGB.
    #[arg(long, value_name = "PAF")]
    hsp_paf: Option<PathBuf>,

    /// How `--hsp-paf` is drawn: red lines, lines coloured by identity, or
    /// grey by identity instead of the dots.
    #[arg(long, value_name = "MODE", default_value = "line",
          value_parser = ["line", "score", "greyscale"], requires = "hsp_paf")]
    hsp_mode: String,

    /// Also write a C-dotter compatible `.dot` binary to this path.
    /// Format 3 (matches the current C-dotter savePlot output).
    #[arg(long, value_name = "PATH")]
    dot: Option<PathBuf>,

    #[command(flatten)]
    large: large_batch::LargeArgs,

    #[command(flatten)]
    genome: genome_batch::GenomeArgs,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum ModeArg {
    Blastn,
    Blastp,
    Blastx,
}

impl ModeArg {
    fn to_core(self) -> BlastMode {
        match self {
            ModeArg::Blastn => BlastMode::Blastn,
            ModeArg::Blastp => BlastMode::Blastp,
            ModeArg::Blastx => BlastMode::Blastx,
        }
    }
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum StrandArg {
    Forward,
    Reverse,
    Both,
}

impl StrandArg {
    fn to_core(self) -> Strand {
        match self {
            StrandArg::Forward => Strand::Forward,
            StrandArg::Reverse => Strand::Reverse,
            StrandArg::Both => Strand::Both,
        }
    }
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum TriangleArg {
    Both,
    Upper,
    Lower,
}

impl TriangleArg {
    fn to_core(self) -> Triangle {
        match self {
            TriangleArg::Both => Triangle::Both,
            TriangleArg::Upper => Triangle::Upper,
            TriangleArg::Lower => Triangle::Lower,
        }
    }
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();
    let cli = Cli::parse();
    match cli.command {
        Command::Batch(args) => run_batch(*args),
        Command::Periodogram(args) => periodogram::run(args),
        Command::FindPeaks(args) => find_peaks::run(args),
    }
}

fn run_batch(args: BatchArgs) -> Result<()> {
    // Genome mode never loads the sequences, so it branches off first.
    if genome_batch::wanted(&args)? {
        return genome_batch::run(&args);
    }
    let mode = args.mode.to_core();
    let matrix = pick_matrix(args.matrix.as_deref(), mode)?;

    // One positional → self-comparison. Two → pairwise.
    let (subject_path, is_self_comparison) = match args.subject {
        Some(p) => (p, false),
        None => (args.query.clone(), true),
    };

    tracing::info!("reading query  {}", args.query.display());
    let query_input = dottir_io::load_sequence_input(&args.query)
        .with_context(|| format!("reading query {}", args.query.display()))?;
    report_input("query", &query_input);

    let subject_input = if is_self_comparison {
        tracing::info!("self-comparison (subject = query)");
        query_input.clone()
    } else {
        tracing::info!("reading subject {}", subject_path.display());
        let s = dottir_io::load_sequence_input(&subject_path)
            .with_context(|| format!("reading subject {}", subject_path.display()))?;
        report_input("subject", &s);
        s
    };
    let (query_input, subject_input, record_sort) = match args.large.sort_records {
        None => (query_input, subject_input, None),
        Some(axis) => {
            if is_self_comparison {
                anyhow::bail!("--sort-records needs two inputs");
            }
            if mode != BlastMode::Blastn {
                anyhow::bail!("--sort-records needs DNA (blastn) inputs");
            }
            let (q, s, info) = large_batch::sort_records(
                &args.large,
                axis.axis(),
                &matrix,
                query_input,
                subject_input,
            )
            .context("sorting records")?;
            (q, s, Some(info))
        }
    };
    if let Some(at) = &args.align_export.export_alignment_at {
        if mode == BlastMode::Blastx {
            anyhow::bail!("--export-alignment-at: not available for BLASTX (DNA vs protein)");
        }
        let text = align_export::export(
            &args.align_export,
            at,
            &query_input.sequence,
            &subject_input.sequence,
        )?;
        print!("{text}");
    }
    let Some(output) = args.output.clone() else {
        return Ok(());
    };
    let query = query_input.sequence.clone();
    let subject = subject_input.sequence.clone();

    // auto-zoom: pick zoom so max(qlen, slen) / zoom <= auto_zoom.
    // For repeat arrays, snap to a nearby record-period divisor so
    // identical records land on the same pixel phase.
    let use_large = args.large.use_large(mode, query.len(), subject.len())?;
    if use_large && args.repeat_uniform_rendering {
        anyhow::bail!("--repeat-uniform-rendering is not supported by the large engine");
    }
    // The large engine targets inputs whose zoom-1 pixelmap could never be
    // allocated; without an explicit zoom, fit the plot to 4000 px.
    let auto_zoom = match args.auto_zoom {
        None if use_large && args.zoom == 1 && query.len().max(subject.len()) > 4000 => {
            tracing::info!("large engine: no --zoom/--auto-zoom given, using --auto-zoom 4000");
            Some(4000)
        }
        a => a,
    };
    let period_hints = zoom_period_hints(&query, &subject, is_self_comparison);
    let zoom = match auto_zoom {
        Some(target) => {
            let base = pick_auto_zoom(query.len(), subject.len(), target);
            let snapped = snap_zoom_to_period_divisor(base, &period_hints, 2.0);
            if snapped != base {
                tracing::info!(
                    "auto_zoom snapped zoom {base} -> {snapped} using period hints {:?}",
                    period_hints
                );
            }
            snapped
        }
        None => args.zoom,
    };
    if zoom != args.zoom {
        tracing::info!("auto_zoom chose zoom={zoom} (--zoom={} ignored)", args.zoom);
    }

    // --watson-only / --crick-only override --strand. Mutually
    // exclusive (clap's `conflicts_with` already enforced that).
    let strand = if args.watson_only {
        Strand::Forward
    } else if args.crick_only {
        Strand::Reverse
    } else {
        args.strand.to_core()
    };
    let mut cfg = PlotConfig {
        mode,
        matrix,
        window_size: args.window,
        zoom,
        pixel_fac: args.pixel_fac,
        strand,
        self_comparison: is_self_comparison,
        triangle: args.triangle.to_core(),
        disable_mirror: args.disable_mirror,
        memory_limit_bytes: args.memory_limit_bytes,
        separate_strand_channels: false,
        reverse_query: args.reverse_query,
        reverse_subject: args.reverse_subject,
    };
    // BLASTP only supports Forward.
    if mode == BlastMode::Blastp {
        cfg.strand = Strand::Forward;
    }

    tracing::info!(
        "computing {}×{} dotplot, mode={:?}, W={:?}, zoom={}",
        query.len(),
        subject.len(),
        mode,
        args.window,
        zoom
    );
    let (plot, large_info, chains) = if use_large {
        let out = large_batch::run(
            &args.large,
            &cfg,
            &query_input,
            &subject_input,
            is_self_comparison,
        )?;
        (out.plot, Some(out.info), out.chains)
    } else {
        let plot = compute_dotplot(query.bytes(), subject.bytes(), &cfg)
            .context("compute_dotplot failed")?;
        (plot, None, Vec::new())
    };
    tracing::info!(
        "wrote pixelmap {}x{} ({} bytes)",
        plot.width,
        plot.height,
        plot.pixels.len()
    );

    // Provenance text chunks for the PNG (spec §4.4.4).
    let dottir_version = env!("CARGO_PKG_VERSION");
    let resolved = format!(
        "mode={:?};matrix={};W={};zoom={};pixel_fac={};strand={:?}",
        cfg.mode,
        cfg.matrix.name,
        plot.params.window_size,
        cfg.zoom,
        plot.params.pixel_fac, // resolved value (may differ from cfg if auto)
        cfg.strand
    );
    let pixelmap_format_version = PIXELMAP_FORMAT_VERSION.to_string();
    let mut text_chunks = vec![
        ("dottir-version".to_string(), dottir_version.to_string()),
        (
            "dottir-pixelmap-format-version".to_string(),
            pixelmap_format_version,
        ),
        ("dottir-parameters".to_string(), resolved),
    ];
    if let Some(info) = &large_info {
        text_chunks.push(("dottir-engine".to_string(), "large".to_string()));
        text_chunks.push((
            "dottir-large-parameters".to_string(),
            format!(
                "format_version={};seed={};two_hit_dist={};max_occ={};extend={};depth={};classes={:?};shown={}",
                info.format_version,
                info.seed,
                info.two_hit_dist,
                info.max_occ,
                info.extend,
                info.depth,
                info.class_thresholds,
                info.shown_classes.join(",")
            ),
        ));
    }

    let (render_pixels, render_w, render_h) = if args.repeat_uniform_rendering {
        let mut render_cfg = cfg.clone();
        render_cfg.zoom = 1;
        tracing::info!("repeat-uniform rendering: computing zoom=1 source raster");
        let source_plot = if plot.params.zoom == 1 {
            plot.clone()
        } else {
            compute_dotplot(query.bytes(), subject.bytes(), &render_cfg)
                .context("repeat-uniform zoom=1 compute_dotplot failed")?
        };
        match dottir_io::text_overlay::downsample_record_tiles_max(
            &source_plot.pixels,
            source_plot.width,
            source_plot.height,
            cfg.zoom,
            &query.records,
            &subject.records,
        ) {
            Some((pixels, w, h)) => {
                tracing::info!(
                    "repeat-uniform rendering downsampled {}x{} -> {}x{} at zoom {}",
                    source_plot.width,
                    source_plot.height,
                    w,
                    h,
                    cfg.zoom
                );
                text_chunks.push((
                    "dottir-rendering".to_string(),
                    "repeat-uniform-boundary-aware".to_string(),
                ));
                (pixels, w, h)
            }
            None => {
                tracing::warn!(
                    "repeat-uniform rendering unavailable; using authoritative pixelmap"
                );
                (plot.pixels.clone(), plot.width, plot.height)
            }
        }
    } else {
        (plot.pixels.clone(), plot.width, plot.height)
    };

    // Resize the render pixelmap to the requested PNG width (preserving
    // aspect ratio). `--width 0` is the opt-out. Cap at a generous
    // upper bound so a stray digit can't OOM. SVG keeps the native
    // resolution (SVG renderers scale themselves).
    const MAX_PNG_WIDTH: u32 = 16384;
    let requested_w = args.width.min(MAX_PNG_WIDTH);
    if requested_w != args.width {
        tracing::warn!("PNG width {} clamped to {MAX_PNG_WIDTH}", args.width);
    }
    // HSPs: the displayed (inverted) image in RGB with the alignments
    // drawn in, at the pixelmap's resolution.
    let hsp_rgb: Option<Vec<u8>> = match &args.hsp_paf {
        None => None,
        Some(path) => {
            if mode == BlastMode::Blastx {
                anyhow::bail!("--hsp-paf: not available for BLASTX (codon columns)");
            }
            let set = dottir_io::hsp::load(path, &query, &subject)
                .with_context(|| format!("reading HSPs {}", path.display()))?;
            tracing::info!(
                "HSPs: {} alignments, {} drawn, {} for records not plotted",
                set.alignments,
                set.hsps.len(),
                set.skipped
            );
            let hsp_mode = dottir_io::hsp::HspMode::parse(&args.hsp_mode).expect("clap-checked");
            let mut rgba: Vec<u8> = render_pixels
                .iter()
                .flat_map(|&v| [255 - v, 255 - v, 255 - v, 255])
                .collect();
            dottir_io::hsp::draw_rgba(
                &mut rgba,
                render_w as usize,
                render_h as usize,
                &set.hsps,
                hsp_mode,
                0,
                0,
                plot.params.zoom as usize,
            );
            text_chunks.push(("hsp-paf".into(), path.display().to_string()));
            text_chunks.push(("hsp-mode".into(), args.hsp_mode.clone()));
            Some(
                rgba.chunks_exact(4)
                    .flat_map(|p| [p[0], p[1], p[2]])
                    .collect(),
            )
        }
    };
    let text_chunk_refs: Vec<(&str, &str)> = text_chunks
        .iter()
        .map(|(k, v)| (k.as_str(), v.as_str()))
        .collect();
    let (png_pixels, png_w, png_h) = if requested_w > 0 {
        dottir_io::text_overlay::resize_nearest(&render_pixels, render_w, render_h, requested_w)
    } else {
        (render_pixels.clone(), render_w, render_h)
    };
    if png_w != render_w {
        tracing::info!(
            "PNG resized {}x{} -> {}x{} (--width {})",
            render_w,
            render_h,
            png_w,
            png_h,
            args.width
        );
    }
    // Coord-space dims for axis labels are the sequence lengths
    // (residues), regardless of any pixelmap upscaling. For zoom > 1
    // the pre-upscale pixelmap already represents `seq_len / zoom`
    // pixels, but the axes still label in residue coordinates.
    //
    // For multi-record FASTAs we also pass per-axis record metadata
    // so each record's id is drawn as a segment label above the top
    // axis / beside the left axis (H1).
    let axis_records_x: Vec<dottir_io::text_overlay::AxisRecord> = query
        .records
        .iter()
        .map(|r| {
            dottir_io::text_overlay::AxisRecord::new(
                r.id.clone(),
                r.range.start as u32,
                r.range.end as u32,
            )
        })
        .collect();
    let axis_records_y: Vec<dottir_io::text_overlay::AxisRecord> = subject
        .records
        .iter()
        .map(|r| {
            dottir_io::text_overlay::AxisRecord::new(
                r.id.clone(),
                r.range.start as u32,
                r.range.end as u32,
            )
        })
        .collect();
    // Larger margin when we have record names to draw, so the extra
    // label row above the ticks has breathing room.
    let png_margin = if query.records.len() > 1 || subject.records.len() > 1 {
        70
    } else {
        50
    };
    match &hsp_rgb {
        None => png_export::write_grayscale_png_with_axes(
            &output,
            png_w,
            png_h,
            query.len() as u32,
            subject.len() as u32,
            &png_pixels,
            png_margin,
            &axis_records_x,
            &axis_records_y,
            &text_chunk_refs,
            true, // raw kernel output → invert for "white bg, dark hits"
        )?,
        Some(rgb) => {
            // Resize each channel like the greyscale image.
            let channel =
                |c: usize| -> Vec<u8> { rgb.iter().skip(c).step_by(3).copied().collect() };
            let resized: Vec<Vec<u8>> = (0..3)
                .map(|c| {
                    if png_w == render_w {
                        channel(c)
                    } else {
                        dottir_io::text_overlay::resize_nearest(
                            &channel(c),
                            render_w,
                            render_h,
                            png_w,
                        )
                        .0
                    }
                })
                .collect();
            let rgb_png: Vec<u8> = (0..resized[0].len())
                .flat_map(|i| [resized[0][i], resized[1][i], resized[2][i]])
                .collect();
            png_export::write_rgb_png_with_axes(
                &output,
                png_w,
                png_h,
                query.len() as u32,
                subject.len() as u32,
                &rgb_png,
                png_margin,
                &axis_records_x,
                &axis_records_y,
                &text_chunk_refs,
            )?
        }
    }
    tracing::info!("wrote {}", output.display());

    if let Some(svg_path) = &args.svg {
        // Chains as vectors over the image, mirrored for a
        // self-comparison like the raster.
        let chain_lines = dottir_io::large::chain_svg_lines(
            &chains,
            plot.params.zoom.max(1) as usize,
            is_self_comparison,
        );
        let rgb = match &hsp_rgb {
            Some(rgb) => Some(rgb.clone()),
            None if !chain_lines.is_empty() => Some(
                render_pixels
                    .iter()
                    .flat_map(|&v| [255 - v, 255 - v, 255 - v])
                    .collect(),
            ),
            None => None,
        };
        match &rgb {
            None => dottir_io::svg_export::write_svg(
                svg_path,
                render_w,
                render_h,
                query.len() as u32,
                subject.len() as u32,
                &render_pixels,
                png_margin,
                &axis_records_x,
                &axis_records_y,
                &text_chunk_refs,
                true, // raw kernel output → invert for "white bg, dark hits"
            ),
            Some(rgb) => dottir_io::svg_export::write_rgb_svg_with_lines(
                svg_path,
                render_w,
                render_h,
                (0, 0),
                query.len() as u32,
                subject.len() as u32,
                rgb,
                png_margin,
                &axis_records_x,
                &axis_records_y,
                &text_chunk_refs,
                &chain_lines,
            ),
        }
        .context("SVG export failed")?;
        tracing::info!("wrote {}", svg_path.display());
    }

    if let Some(pdf_path) = &args.pdf {
        use dottir_io::pdf_export::{write_pdf, PdfPixels};
        use dottir_io::text_overlay::{
            compose_image_with_axes, compose_rgb_image_with_axes, inverted,
        };
        let (qlen, slen) = (query.len() as u32, subject.len() as u32);
        match &hsp_rgb {
            None => {
                let (img, w, h) = compose_image_with_axes(
                    &inverted(&render_pixels),
                    render_w,
                    render_h,
                    qlen,
                    slen,
                    png_margin,
                    &axis_records_x,
                    &axis_records_y,
                );
                write_pdf(pdf_path, w, h, PdfPixels::Grey(&img), &text_chunk_refs)
            }
            Some(rgb) => {
                let (img, w, h) = compose_rgb_image_with_axes(
                    rgb,
                    render_w,
                    render_h,
                    qlen,
                    slen,
                    png_margin,
                    &axis_records_x,
                    &axis_records_y,
                );
                write_pdf(pdf_path, w, h, PdfPixels::Rgb(&img), &text_chunk_refs)
            }
        }
        .with_context(|| format!("writing {}", pdf_path.display()))?;
        tracing::info!("wrote {}", pdf_path.display());
    }

    if let (Some(_), Some(_)) = (&args.dot, &large_info) {
        tracing::warn!(
            ".dot stores the shown classes merged; repeat-class channels are not preserved"
        );
    }
    if let Some(dot_path) = &args.dot {
        // Construct a format-3 .dot. Matrix is only 24×24 for proteins;
        // for DNA we pad the 4×4 matrix into the top-left of a 24×24
        // block (the rest stays zero — C-dotter doesn't care about
        // those cells when reading DNA `.dot` files because the inner
        // loop only indexes the first 4 rows/cols).
        let mut matrix24 = vec![0_i32; 24 * 24];
        let n = cfg.matrix.size();
        for r in 0..n {
            for c in 0..n {
                matrix24[r * 24 + c] = cfg.matrix.get(r, c);
            }
        }
        let dot = dottir_io::dot_format::DotFile {
            format: 3,
            zoom: cfg.zoom as f64,
            width: plot.width,
            height: plot.height,
            pixel_fac: Some(plot.params.pixel_fac as i32),
            sliding_window_size: Some(plot.params.window_size as i32),
            matrix_name: Some(cfg.matrix.name.clone()),
            matrix: Some(matrix24),
            pixels: plot.pixels.clone(),
        };
        dot.write_format3(dot_path).context(".dot export failed")?;
        tracing::info!("wrote {}", dot_path.display());
    }

    if !args.no_sidecar {
        let sidecar_path = sidecar_path(&output);
        let sidecar = ParamsSidecar {
            dottir: DottirInfo {
                version: dottir_version.to_string(),
                git_sha: option_env!("GIT_SHA").map(|s| s.to_string()),
                pixelmap_format_version: PIXELMAP_FORMAT_VERSION,
            },
            query: input_info(
                &args.query,
                &query_input.bytes,
                &query.records,
                query.bytes(),
            )?,
            subject: input_info(
                &subject_path,
                &subject_input.bytes,
                &subject.records,
                subject.bytes(),
            )?,
            plot: PlotParamsInfo {
                mode: format!("{:?}", cfg.mode),
                matrix: cfg.matrix.name.clone(),
                strand: format!("{:?}", cfg.strand),
                window_size: plot.params.window_size,
                zoom: cfg.zoom,
                pixel_fac: plot.params.pixel_fac,
                self_comparison: cfg.self_comparison,
                karlin: plot.params.karlin.map(|k| KarlinInfo {
                    lambda: k.lambda,
                    k: k.k,
                    h: k.h,
                    expected_residue_score: k.expected_residue_score,
                    expected_msp_score: k.expected_msp_score,
                    predicted_msp_length: k.predicted_msp_length,
                }),
                width: plot.width,
                height: plot.height,
            },
            host: HostInfo {
                hostname: dottir_io::params::hostname(),
                os: std::env::consts::OS.to_string(),
                timestamp_utc: now_iso_utc(),
            },
            large: large_info.clone(),
            genome: None,
            record_sort: record_sort.map(|r| dottir_io::params::RecordSortInfo {
                axis: match r.axis {
                    SortAxis::Vertical => "vertical",
                    SortAxis::Horizontal => "horizontal",
                }
                .to_string(),
                chains: r.chains,
                order: r.order,
            }),
        };
        std::fs::write(&sidecar_path, sidecar.to_toml()?)?;
        tracing::info!("wrote sidecar {}", sidecar_path.display());
    }

    Ok(())
}

fn pick_matrix(name: Option<&str>, mode: BlastMode) -> Result<ScoreMatrix> {
    match (name, mode) {
        (None, BlastMode::Blastn) => Ok(ScoreMatrix::dna_identity()),
        (None, _) => Ok(ScoreMatrix::blosum62()),
        (Some(n), BlastMode::Blastn) => {
            if n.eq_ignore_ascii_case("DNA") || n.eq_ignore_ascii_case("DNA+5/-4") {
                Ok(ScoreMatrix::dna_identity())
            } else {
                anyhow::bail!("matrix {n:?} not valid for BLASTN — use DNA+5/-4")
            }
        }
        (Some(n), _) => {
            ScoreMatrix::by_name(n).ok_or_else(|| anyhow::anyhow!("unknown matrix {n:?}"))
        }
    }
}

/// Log what a positional input turned out to be. Batch mode has no
/// annotation overlay, so features found in a GFF3 input are counted
/// and dropped — saying so beats silence when a user passes an
/// annotated GFF3 and sees a plain dotplot.
fn report_input(role: &str, input: &dottir_io::SeqInput) {
    match input.kind {
        dottir_io::InputKind::Fasta => tracing::info!(
            "{role}: FASTA, {} records, {} residues",
            input.sequence.records.len(),
            input.sequence.len(),
        ),
        dottir_io::InputKind::Gff3WithFasta => {
            let n_features = input.annotations.as_ref().map_or(0, |a| a.features.len());
            tracing::info!(
                "{role}: GFF3 with embedded FASTA, {} records, {} residues, \
                 {n_features} features (not rendered in batch mode)",
                input.sequence.records.len(),
                input.sequence.len(),
            );
        }
    }
}

fn input_info(
    path: &Path,
    bytes: &[u8],
    recs: &[dottir_io::RecordSpan],
    sequence: &[u8],
) -> Result<InputInfo> {
    Ok(InputInfo {
        path: path.display().to_string(),
        sha256: sha256_bytes(bytes),
        size_bytes: bytes.len() as u64,
        n_records: recs.len(),
        total_residues: sequence.len(),
    })
}

fn zoom_period_hints(
    query: &dottir_io::Sequence,
    subject: &dottir_io::Sequence,
    is_self_comparison: bool,
) -> Vec<usize> {
    if is_self_comparison {
        return query.record_period_hint().into_iter().collect();
    }
    match (query.record_period_hint(), subject.record_period_hint()) {
        (Some(q), Some(s)) => vec![q, s],
        _ => Vec::new(),
    }
}

fn sidecar_path(output: &Path) -> PathBuf {
    let mut s = output.as_os_str().to_owned();
    s.push(".params.toml");
    PathBuf::from(s)
}

fn now_iso_utc() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Cheap and correct enough for a sidecar (a few-second drift is fine).
    // Avoids pulling in chrono.
    iso_from_unix(secs as i64)
}

fn iso_from_unix(secs: i64) -> String {
    // Days from epoch to Y-M-D using a public-domain Howard Hinnant
    // algorithm (date.h). Civil time = UTC; no timezone handling.
    let z = secs.div_euclid(86_400);
    let secs_of_day = secs.rem_euclid(86_400);
    let (h, m, s) = (
        (secs_of_day / 3600) as u32,
        ((secs_of_day / 60) % 60) as u32,
        (secs_of_day % 60) as u32,
    );
    let z = z + 719_468; // shift epoch to 0000-03-01
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m_civil = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
    let y_civil = y + if m_civil <= 2 { 1 } else { 0 };
    format!("{y_civil:04}-{m_civil:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}
