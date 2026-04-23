# B9: render the Quarto report.
#
# Drives `quarto::quarto_render()` on `report/template.qmd` with the active
# `config.yaml` path (and `project_root`) passed through `execute_params`.
# Quarto emits `template.html` next to the source; we copy it to the contracted
# Snakemake output path. PDF rendering is attempted when the active config
# requests it but fails gracefully if the LaTeX toolchain is absent.

log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

suppressPackageStartupMessages({
  library(quarto)
  library(yaml)
})

template_path <- snakemake@input[["template"]]
html_out      <- snakemake@output[["html"]]
formats       <- snakemake@params[["formats"]]

if (is.null(formats) || length(formats) == 0) formats <- c("html")

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
template_src <- normalizePath(template_path, winslash = "/", mustWork = TRUE)

# Quarto writes its sidecar (template.html, .quarto/, template_files/, etc.)
# next to the input. In Docker the baked template lives under /app/report,
# which is owned by mambauser — after the entrypoint drops to the host UID
# it becomes read-only, so rendering in place fails. Stage the report dir
# (and the sibling VERSION file, which the template reads via ../VERSION)
# into a writable location under the output tree before rendering.
template_src_dir <- dirname(template_src)
stage_root       <- file.path(dirname(snakemake@output[["html"]]), "_render")
stage_report     <- file.path(stage_root, basename(template_src_dir))
unlink(stage_root, recursive = TRUE, force = TRUE)
dir.create(stage_report, showWarnings = FALSE, recursive = TRUE)
file.copy(list.files(template_src_dir, full.names = TRUE, all.files = TRUE,
                     no.. = TRUE),
          stage_report, recursive = TRUE)
version_src <- file.path(dirname(template_src_dir), "VERSION")
if (file.exists(version_src)) {
  file.copy(version_src, file.path(stage_root, "VERSION"), overwrite = TRUE)
}
template_abs <- normalizePath(file.path(stage_report, basename(template_src)),
                              winslash = "/", mustWork = TRUE)

# Snakemake's `@config` is the merged configuration (base configfile + any
# `--configfile` overrides + `--config` flags). Write it to a temp YAML so the
# Quarto template sees the exact config that drove the run, regardless of
# which `--configfile` was passed.
dir.create(dirname(html_out), showWarnings = FALSE, recursive = TRUE)
merged_cfg_path <- file.path(dirname(html_out), "_merged_config.yaml")
yaml::write_yaml(as.list(snakemake@config), merged_cfg_path)
config_abs <- normalizePath(merged_cfg_path, winslash = "/", mustWork = TRUE)

message(sprintf("render_report: project_root=%s", project_root))
message(sprintf("render_report: template=%s", template_abs))
message(sprintf("render_report: config=%s", config_abs))
message(sprintf("render_report: formats=%s", paste(formats, collapse = ", ")))

dir.create(dirname(html_out), showWarnings = FALSE, recursive = TRUE)

execute_params <- list(
  config_path  = config_abs,
  project_root = project_root
)

ensure_latex_wrappers <- function() {
  make_wrapper <- function(wrapper, target) {
    wrapper_path <- Sys.which(wrapper)
    target_path <- Sys.which(target)
    if (nzchar(wrapper_path) || !nzchar(target_path)) return(invisible(FALSE))

    bin_dir <- dirname(target_path)
    dest <- file.path(bin_dir, wrapper)
    ok <- tryCatch({
      if (!file.exists(dest)) file.symlink(basename(target_path), dest)
      file.exists(dest)
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok) {
      message(sprintf("render_report: created LaTeX wrapper %s -> %s",
                      dest, basename(target_path)))
    } else {
      message(sprintf("render_report: could not create LaTeX wrapper %s", dest))
    }
    invisible(ok)
  }

  make_wrapper("xelatex", "xetex")
  make_wrapper("pdflatex", "pdftex")
}

ensure_latex_wrappers()

render_fmt <- function(fmt, sidecar_ext) {
  message(sprintf("render_report: rendering format=%s", fmt))
  quarto::quarto_render(
    input          = template_abs,
    output_format  = fmt,
    execute_params = execute_params,
    quiet          = FALSE
  )
  sidecar <- sub("\\.qmd$", paste0(".", sidecar_ext), template_abs)
  if (!file.exists(sidecar)) {
    stop(sprintf("Expected Quarto output not found: %s", sidecar))
  }
  sidecar
}

# HTML is required (Snakemake output contract).
html_sidecar <- render_fmt("html", "html")
file.copy(html_sidecar, html_out, overwrite = TRUE)
unlink(html_sidecar)
message(sprintf("render_report: wrote %s", html_out))

# PDF is optional.
if ("pdf" %in% formats) {
  tryCatch({
    pdf_sidecar <- render_fmt("pdf", "pdf")
    pdf_out <- sub("\\.html$", ".pdf", html_out)
    file.copy(pdf_sidecar, pdf_out, overwrite = TRUE)
    unlink(pdf_sidecar)
    message(sprintf("render_report: wrote %s", pdf_out))
  }, error = function(e) {
    message("render_report: PDF render failed, continuing with HTML only. Reason: ",
            conditionMessage(e))
  })
}

unlink(stage_root, recursive = TRUE, force = TRUE)

invisible(NULL)
