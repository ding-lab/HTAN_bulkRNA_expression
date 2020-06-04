## Prepare annotation
The workflow to create a Ensembl v100 (GENCODE v34) annotation compatible with GDC genome reference (GRCh38.d1.vd1).

1. Download all upstream sources by `snakemake all_upstream_sources`
2. Run `analysis/03_make_gene_info.Rmd` to generate the gene info table