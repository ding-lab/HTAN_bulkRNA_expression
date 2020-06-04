import csv
from collections import namedtuple
from pathlib import Path
import logging

# List of sample names to run the pipeline
SAMPLE_LIST_PTH = config['sample_list']
# The mapping of sample name to other information.
BAM_MAP_PTH = config['bam_map']
GENE_GTF_PTH = config['gtf']  # gencode.v34.annotation.gtf
GENE_INFO_PTH = config['gene_info']  # gencode.gene.info.v34.tsv
WORKFLOW_ROOT = config['workflow_root']  # Path to this repository

_logger = logging.getLogger(__name__)

# Read all the cases to process
with open(SAMPLE_LIST_PTH) as f:
    SAMPLES = f.read().splitlines()
if len(SAMPLES) != len(set(SAMPLES)):
    _logger.error('There are duplicated samples in the sample list!')
SAMPLES = set(SAMPLES)


# Select all the available samples of the selected cases.
# SAMPLE_INFO is a mapping of:
#   sample_name -> SampleInfo(case, sample_type, disease, UUID, BAM_src_path)
SampleInfo = namedtuple('SampleInfo', 'case, disease, uuid, bam_pth')
SAMPLE_INFO = {}
with open(BAM_MAP_PTH) as f:
    reader = csv.DictReader(f, dialect='excel-tab')
    for row in reader:
        sample_name = row['# sample_name']
        if sample_name not in SAMPLES:
            continue
        is_genomic_rna_bam = (
            row['experimental_strategy'] == 'RNA-Seq' and
            row['data_format'] == 'BAM' and
            row['reference'] == 'hg38' and
            row['data_path'].endswith('.rna_seq.genomic.gdc_realn.bam')
        )
        if not is_genomic_rna_bam:
            _logger.warning(f'{sample_name} is not a genomic RNA-Seq BAM')

        bam_pth = Path(row['data_path'])
        SAMPLE_INFO[sample_name] = SampleInfo(
            row['case'], row['disease'], row['UUID'], bam_pth
        )


def find_sample_fq_path(wildcards):
    """Find the FASTQ file paths of a given sample."""
    # FIXME: implement this
    bam_pth = str(SAMPLE_INFO[wildcards.sample].bam_pth)
    return {
        'r1_fq': bam_pth,
        'r2_fq': bam_pth + '.bai',
    }


rule star_align:
    """STAR align one sample."""
    output:
        # unsorted_bam='star/{sample}/Aligned.out.bam',
        sorted_bam='star/{sample}/Aligned.sortedByCoord.out.bam',
        chimeric_sam='star/{sample}/Chimeric.out.sam',
        chimeric_junction='star/{sample}/Chimeric.out.junction',
        quant_tx_bam='star/{sample}/Aligned.toTranscriptome.out.bam',
        quant_gene_count_tab='star/{sample}/ReadsPerGene.out.tab',
        sj_count_tab='star/{sample}/SJ.out.tab',
    input: unpack(find_sample_fq_path)        
    params:
        star_ix=STAR_INDEX_FOLDER,
        star_gtf=STAR_GTF,
        out_folder='star/{sample}/'
    log: 'logs/star/{sample}.log'
    threads: 6
    shell:
        "STAR "
        "--readFilesIn {input.r1_fq} {input.r2_fq} "
        "--alignIntronMax 1000000 "
        "--alignIntronMin 20 "
        "--alignMatesGapMax 1000000 "
        "--alignSJDBoverhangMin 1 "
        "--alignSJoverhangMin 8 "
        "--alignSoftClipAtReferenceEnds Yes "
        "--chimJunctionOverhangMin 15 "
        "--chimMainSegmentMultNmax 1 "
        "--chimOutType Junctions SeparateSAMold WithinBAM SoftClip "
        "--chimOutJunctionFormat 1 "
        "--chimSegmentMin 15 "
        "--genomeDir {params.star_ix} "
        "--genomeLoad NoSharedMemory "
        "--limitBAMsortRAM 0 "
        "--limitSjdbInsertNsj 1200000 "
        "--outFileNamePrefix {params.out_folder} "
        "--outFilterIntronMotifs None "
        "--outFilterMatchNminOverLread 0.33 "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverLmax 0.1 "
        "--outFilterMultimapNmax 20 "
        "--outFilterScoreMinOverLread 0.33 "
        "--outFilterType BySJout "
        "--outSAMattributes NH HI AS nM NM ch "
        "--outSAMstrandField intronMotif "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMunmapped Within "
        "--quantMode TranscriptomeSAM GeneCounts "
        "--readFilesCommand zcat "
        "--runThreadN {threads} "
        "--twopassMode Basic "
        "> {log}"


rule samtools_index_bam:
    """Index a sorted BAM by samtools."""
    output: '{name}.bam.bai'
    input: '{name}.bam'
    resources:
        io_heavy=1
    shell: 'samtools index {input} {output}'


rule star_align_all_samples:
    """Align all RNA-seq samples."""
    input:
        all_sorted_bams=expand(rules.star_align.output.sorted_bam, sample=SAMPLES),
        all_sorted_bam_bais=expand(rules.star_align.output.sorted_bam + '.bai', sample=SAMPLES)


def find_sample_bam_path(wildcards):
    """Find the BAM and BAM index file paths of a given sample."""
    bam_pth = str(SAMPLE_INFO[wildcards.sample].bam_pth)
    return {
        'bam': bam_pth,
        'bai': bam_pth + '.bai',
    }


rule featurecounts_unstranded_readcount:
    """Readcount by featureCounts (unstranded)."""
    output: count_tsv=temp('featurecounts_unstranded_readcount/{sample}.tsv')
    input: unpack(find_sample_bam_path)
    log: 'logs/featurecounts_unstranded/{sample}.log'
    params:
        gtf=GENE_GTF_PTH
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 16000 + 16000 * (attempt - 1)
    threads: 16
    group: "featurecounts"
    shell:
        'featureCounts '
        '-g gene_id '  # feature id (-i in htseq)
        '-t exon '  # feature type (-t in htseq)
        '-T {threads} '
        '-Q 10 '  # htseq set this minimal mapping quality by default
        '-p '  # pair-end reads are considered one fragment; default HTSeq behavior
        '-B '  # both reads of a read pair need to be mapped
        # '-s 2 '  # reversely stranded
        '-a {params.gtf} '
        '-o {output.count_tsv} {input.bam} 2> {log}'


rule compress_featurecounts:
    """Shrink and compress featureCounts output."""
    output: 'featurecounts_unstranded_readcount/{sample}.tsv.gz'
    input: rules.featurecounts_unstranded_readcount.output.count_tsv
    threads: 2
    group: "featurecounts"
    shell: 'python {WORKFLOW_ROOT}/scripts/shrink_featurecounts.py {input} | gzip -9 -c > {output}'


rule generate_fpkm:
    """Generate FPKM and FPKM-UQ from the readcount."""
    output: fpkm='readcount_and_fpkm/{sample}.tsv.gz'
    input: rc=rules.compress_featurecounts.output[0],
           gene_info=GENE_INFO_PTH
    shell: 'python {WORKFLOW_ROOT}/scripts/gen_fpkm.py {input.gene_info} {input.rc} {output.fpkm}'


rule all_featurecounts_stranded_readcount:
    input:
        counts=expand(rules.compress_featurecounts.output[0], sample=SAMPLES)


rule all_fpkms:
    input: fpkms=expand(rules.generate_fpkm.output.fpkm, sample=SAMPLES)


rule make_analysis_summary:
    """Generate the analysis summary table."""
    input: rules.all_fpkms.input.fpkms
    output: analysis_summary='analysis_summary.dat'
    run:
        with open(output.analysis_summary, 'w') as f:
            writer = csv.writer(f, dialect='excel-tab', lineterminator='\n')
            # Write column header
            cols = ['# case', 'disease',
                    'data_path', 'file_format',
                    'sample_name', 'sample_uuid']
            writer.writerow(cols)

            for sample, info in SAMPLE_INFO.items():
                count_tsv_pth = Path(
                    rules.generate_fpkm.output.fpkm.format(sample=sample)
                ).resolve(strict=True)
                writer.writerow([
                    info.case, info.disease,
                    str(count_tsv_pth), 'TSV',
                    sample, info.uuid
                ])
