## Prepare annotation
The workflow to create a Ensembl v100 (GENCODE v34) annotation compatible with GDC genome reference (GRCh38.d1.vd1).


### Download all upstream sources

    mkdir upstream_sources
    curl -L -o upstream_sources/gencode.v34.annotation.gtf.gz \
        https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz
    curl -L -o upstream_sources/EnsDb.Hsapiens.v100.sqlite \
        http://s3.amazonaws.com/annotationhub/AHEnsDbs/v100/EnsDb.Hsapiens.v100.sqlite


### Run R notebooks
Run `analysis/03_make_gene_info.Rmd` to generate the gene info table.
