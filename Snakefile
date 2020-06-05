from collections import namedtuple, defaultdict
import csv
from pathlib import Path
import re

from snakemake.logging import logger

# List of sample names to run the pipeline
SAMPLE_LIST_PTH = config['sample_list']
# The mapping of sample name to other information.
FILE_MAP_PTH = config['file_map']
GENE_GTF_PTH = config['gtf']  # gencode.v34.annotation.gtf
GENE_INFO_PTH = config['gene_info']  # gencode.gene.info.v34.tsv
WORKFLOW_ROOT = config['workflow_root']  # Path to this repository
STAR_INDEX_FOLDER = config['star_index']  # Path to the STAR index

# Read all the sample names (HTAN_Specimen_ID) to process
logger.info('Reading sample list and file map ...')
with open(SAMPLE_LIST_PTH) as f:
    SAMPLES = f.read().splitlines()
if len(SAMPLES) != len(set(SAMPLES)):
    logger.error('There are duplicated samples in the sample list!')
SAMPLES = set(SAMPLES)

# FASTQ map structcure
# {
#     sample1_id: {
#         read_group1_id: {R1: fq_pth, R2: fq_pth},
#         read_group2_id: {R1: fq_pth, R2: fq_pth},
#         ...
#     },
#     sample2_id: { ... },
#     ...
# }
FASTQ_MAP = defaultdict(lambda: defaultdict(lambda: {'R1': None, 'R2': None}))

# SAMPLE_INFO structure:
# {sample_id: SampleInfo(patient_id, cancer_type, sample_type, tissue)}
SampleInfo = namedtuple('SampleInfo', 'patient_id, cancer_type, sample_type, tissue')
SAMPLE_INFO = {}
with open(FILE_MAP_PTH) as f:
    reader = csv.DictReader(f, dialect='excel-tab')
    for row in reader:
        sample_id = row['HTAN_Specimen_ID']
        is_bulk_rna_fastq = (
            row['Experiment_Type'] == 'Bulk_RNA-Seq' and
            row['Data_Format'] == 'FASTQ'
        )
        if sample_id not in SAMPLES or not is_bulk_rna_fastq:
            continue

        # Update FASTQ_MAP
        fq_pth = Path(row['Path'])
        read_group_id = re.match(r'^(\S+)_R[12].fastq.gz$', fq_pth.name).group(1)
        read_strand = row['Data_Format_Details']
        if read_strand not in ['R1', 'R2']:
            logger.error(f'This row of FASTQ entry is malformatted: {row}')
            raise ValueError(f'Invalid FASTQ format in the file map for {sample_id}')
        FASTQ_MAP[sample_id][read_group_id][read_strand] = fq_pth

        # Upate SAMPLE_INFO
        SAMPLE_INFO[sample_id] = SampleInfo(
            row['Patient_ID'], row['Cancer_Type'], row['Sample_Type'], row['Tissue']
        )

# Validate the FASTQ_MAP
for sample, read_groups in FASTQ_MAP.items():
    if not read_groups:
        raise ValueError(f'{sample} has no read groups')
    for rg, fqs in read_groups.items():
        if any(fq_pth is None for fq_pth in fqs.values()):
            raise ValueError(f'{rg} of {sample} has unpaired FASTQs')

# Samples without FASTQ
samples_without_fq = set(SAMPLES) - set(FASTQ_MAP)
if samples_without_fq:
    raise ValueError(f"These samples have no FASTQs: {' '.join(samples_without_fq)}")

# Summarize the FASTQ_MAP
fastq_map_summary = '\n'.join(
    f'\t{sample}: {len(read_groups)} read groups'
    for sample, read_groups in FASTQ_MAP.items()
)
logger.info(
    f'Built FASTQ mapping of {len(FASTQ_MAP):,d} samples:\n'
    f'{fastq_map_summary}'
)
logger.info('Sample list and FASTQ mapping pass all checks!\n')



def find_sample_fq_path(wildcards):
    """Find the FASTQ file paths of a given sample."""
    read_groups = FASTQ_MAP[wildcards.sample]
    return {
        'r1_fq': [str(rg['R1']) for rg in read_groups.values()],
        'r2_fq': [str(rg['R2']) for rg in read_groups.values()],
    }


def merge_input_fqs(wildcards, input):
    fqs = find_sample_fq_path(wildcards)
    r1_fq_str = ','.join(input.r1_fq)
    r2_fq_str = ','.join(fqs['r2_fq'])
    return f'{r1_fq_str} {r2_fq_str}'


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
        star_gtf=GENE_GTF_PTH,
        input_str=merge_input_fqs,
        out_folder='star/{sample}/'
    log: 'logs/star/{sample}.log'
    threads: 16
    resources:
        io_heavy=1,
        mem_mb=lambda wildcards, attempt: 32000 + 16000 * (attempt - 1)
    shell:
        "STAR "
        "--readFilesIn {params.input_str} "
        # Most parameters follow GDC
        "--alignIntronMax 1000000 "
        "--alignIntronMin 20 "
        "--alignMatesGapMax 1000000 "
        "--alignSJDBoverhangMin 1 "
        "--alignSJoverhangMin 8 "
        "--alignSoftClipAtReferenceEnds Yes "

        # Follow arriba's recommendation regarding chimera parameters
        # Ref: https://arriba.readthedocs.io/en/latest/workflow/
        "--chimJunctionOverhangMin 10 "
        "--chimMainSegmentMultNmax 1 "
        "--chimOutType Junctions SeparateSAMold WithinBAM SoftClip "
        "--chimOutJunctionFormat 1 "
        "--chimSegmentMin 10 "
        "--chimScoreMin 1" 
        "--chimScoreDropMax 30 "
        "--chimScoreJunctionNonGTAG 0 "
        "--chimScoreSeparation 1 "
        "--alignSJstitchMismatchNmax 5 -1 5 5 "
        "--chimSegmentReadGapMax 3 "
    
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


rule featurecounts_unstranded_readcount:
    """Readcount by featureCounts (unstranded)."""
    output: count_tsv=temp('featurecounts_unstranded_readcount/{sample}.tsv')
    input:
        bam=rules.star_align.output.sorted_bam,
        bai=rules.star_align.output.sorted_bam + '.bai'
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
    """featureCounts of all samples."""
    input: counts=expand(rules.compress_featurecounts.output[0], sample=SAMPLES)


rule all_fpkms:
    """FPKM TSVs of all samples."""
    input: fpkms=expand(rules.generate_fpkm.output.fpkm, sample=SAMPLES)


rule make_analysis_summary:
    """Generate the analysis summary table."""
    input: 
        tsvs=rules.all_fpkms.input.fpkms,
        genomic_bams=rules.star_align_all_samples.input.sorted_bam
    output: analysis_summary='analysis_summary.dat'
    run:
        with open(output.analysis_summary, 'w') as f:
            writer = csv.writer(f, dialect='excel-tab', lineterminator='\n')
            # Write column header
            cols = ['# case', 'disease',
                    'data_path', 'file_format',
                    'sample_id', 'sample_uuid']
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
