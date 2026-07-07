# B2m CD40 analysis

This repository contains the computational analysis (single-cell RNA-seq and TCGA bulk RNA-seq) supporting a manuscript currently under review. It is built as a [Quarto book](https://quarto.org/docs/books/); the rendered chapters walk through the analysis end to end.

## Structure

- `book/` — analysis chapters, in the order they appear in the book (see `_quarto.yml`):
  - **Single cell analysis**: Seurat object setup, cell type annotation, T/myeloid/NK/γδT subclustering, doublet removal, pseudobulk analysis, and the core B2m-KO-vs-WT and CD40-agonist-vs-isotype comparisons, ending in `figures.qmd` which generates the manuscript figures.
  - **TCGA analysis**: derivation of an inflammatory monocyte signature and its application to TCGA bulk tumor data (signature scoring, survival analysis).
  - **Miscellaneous**: supporting exports (e.g. cell-type proportions).
- `lib/` — shared setup and helper R code (Seurat helpers, plotting themes/palettes, Single Cell Portal export utilities) sourced by the book chapters.
- `renv.lock` — pinned R package versions (R 4.2.1) for reproducing the environment via [renv](https://rstudio.github.io/renv/).

## Reproducing the analysis

1. Clone the repository and open `b2m-cd40-quarto.Rproj` in RStudio (or launch R in this directory — `.Rprofile` will activate `renv` automatically).
2. Restore the R environment: `renv::restore()`.
3. Download the source data:
   - Single-cell sequencing data (10x CellRanger output) from NCBI GEO.
   - Processed single-cell data is also available on the Broad Institute [Single Cell Portal](https://singlecell.broadinstitute.org/).
   - See `book/setup_sc.qmd` for details on expected data layout (`data/cellranger_output_seq1`, `data/cellranger_output_seq2`).
   - TCGA data used in `book/15_prelim_tcga_analysis.qmd` and `book/16_full_tcga_analysis.qmd` is public bulk RNA-seq/clinical data, downloaded separately (see those chapters for details).
4. Render the book: `quarto render`.

Data and intermediate results (`data/`, `results/`) are not tracked in this repository — they must be downloaded/generated locally per the steps above.

## Contact

Please contact a corresponding author for any questions, comments, or concerns regarding the paper.

## License

Apache License 2.0 — see [`LICENSE`](LICENSE).
