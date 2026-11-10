# Hair Follicle Regeneration Atlas

A static gene expression browser for paired low- and high-density mouse hair follicle regeneration sites. Search an official mouse gene symbol and switch between whole skin and the hair follicle subset. The L/H view uses the same expression scale for both groups.

The viewer uses saved UMAP coordinates and final annotations. RNA counts are normalized to 10,000 counts per cell and log1p-transformed. For distribution, expression is rounded to the nearest 1/4096 of a log-normalized unit, and coordinates are stored as 32-bit floats. The display caps each gene's color scale at the 99th percentile of positive values in the selected atlas.

## Local preview

Serve the published directory using any static HTTP server. For example, run `python -m http.server 8765 --directory docs` and open `http://localhost:8765`. The browser needs native gzip decompression support.

## Data format

Each atlas has a gene index, cell type palette, anonymous UMAP coordinates with cell-type and L/H codes, and gzip-compressed expression shards of up to 24 genes. Only the selected gene's shard is fetched. The source RDS objects and raw sequencing files are separate from this repository.

`tools/export.R` exports the two frozen analysis objects and verifies annotation counts against the figure source tables. Run `node tools/verify.mjs` to check the published payloads.
