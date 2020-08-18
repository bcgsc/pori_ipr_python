# Variant Input Files

The reference documentation for `genomic_report.inputs` contains the function signatures. This
document is focused on users building the files. All variant files for the genomic report
use a tab delimited format. The column names, types and example inputs are shown below.

## Expression Variant Data

### Required Columns

| Column          | Type   | Description                                                            |
| --------------- | ------ | ---------------------------------------------------------------------- |
| gene            | string | the gene name (or source identifier)                                   |
| kbCategory      | string | the graphkb expression variant vocabulary term this variant belongs to |
| expressionState | string | the variant used for display in the report                             |


### Optional Columns


| Column                | Type  | Description                                                                                        |
| --------------------- | ----- | -------------------------------------------------------------------------------------------------- |
| rnaReads              | float |                                                                                                    |
| rpkm                  | float | reads per kilobase of transcript, per million mapped reads                                         |
| tpm                   | float | transcript per million                                                                             |
| diseasePercentile     | float | the percentile with respect to the disease expression comparator cohort                            |
| diseaseKIQR           | float | the kIQR with respect to the disease expression comparator cohort                                  |
| diseaseZScore         | float | the zscore with respect to the disease expression comparator cohort                                |
| diseaseFoldChange     | float | the fold change with respect to the median of the disease expression comparator cohort             |
| diseaseQC             | float |                                                                                                    |
| primarySitePercentile | float | the percentile with respect to the normal primary site expression comparator cohort                |
| primarySiteKIQR       | float | the kIQR with respect to the normal primary site expression comparator cohort                      |
| primarySiteZScore     | float | the zscore with respect to the normal primary site expression comparator cohort                    |
| primarySiteFoldChange | float | the fold change with respect to the median of the normal primary site expression comparator cohort |
| primarySiteQC         | float |                                                                                                    |
| biopsySitePercentile  | float | the percentile with respect to the normal biopsy site expression comparator cohort                 |
| biopsySiteKIQR        | float | the kIQR with respect to the normal biopsy site expression comparator cohort                       |
| biopsySiteZScore      | float | the zscore with respect to the normal biopsy site expression comparator cohort                     |
| biopsySiteFoldChange  | float | the fold change with respect to the median of the normal biopsy site expression comparator cohort  |
| biopsySiteQC          | float |                                                                                                    |

## Small Mutation Data

Small mutations are composed of indels and single nucleotide variants.

### Required Columns

| Column        | Type    | Example     | Description                                |
| ------------- | ------- | ----------- | ------------------------------------------ |
| chromosome    | string  | X           | the chromosome                             |
| startPostion  | integer | 1234        | the genomic start position of this variant |
| endPosition   | integer | 1234        | the genomic end position of this variant   |
| refSeq        | string  | A           | the reference sequence                     |
| altSeq        | string  | C           | the alternate/new sequence                 |
| gene          | string  | KRAS        | the gene name                              |
| proteinChange | string  | p.G12D      | the HGVS protein notation                  |
| transcript    | string  | ENST00001.2 | the transcript name                        |


### Optional Columns

| Column         | Type    | Example                | Description                                                                |
| -------------- | ------- | ---------------------- | -------------------------------------------------------------------------- |
| zygosity       | string  | het                    |                                                                            |
| tumourRefCount | integer | 1                      | the number of reference reads in the tumour genome                         |
| tumourAltCount | integer | 1                      | the number of alternate reads in the tumour genome supporting the mutation |
| tumourDepth    | integer | 1                      | the total number of reads at this position in the tumour genome            |
| rnaRefCount    | integer | 1                      | the number of reference reads in the rna                                   |
| rnaAltCount    | integer | 1                      | the number of alternate reads in the rna  supporting the mutation          |
| rnaDepth       | integer | 1                      | the total number of reads at this position in the rna                      |
| normalRefCount | integer | 1                      | the number of reference reads in the normal genome                         |
| normalAltCount | integer | 1                      | the number of alternate reads in the normal genome supporting the mutation |
| normalDepth    | integer | 1                      | the total number of reads at this position in the normal genome            |
| detectedIn     | string  | DNA/RNA                | the sample types this variant was detected in                              |
| hgvsCds        | string  | `ENST0001:c.1234+3A>G` | HGVS coding sequence notation for this variant                             |
| hgvsGenomic    | string  | `1:g.1234A>G`          | HGVS genomic notation for this variant                                     |
| hgvsProtein    | string  | `KRAS:p.G12D`          | HGVS protein notation for this variant                                     |
| ncbiBuild      | string  | GRCh37                 | the genome reference assembly build version                                |

## Copy Variant Data


### Required Columns

| Column     | Type   | Example | Description                                                      |
| ---------- | ------ | ------- | ---------------------------------------------------------------- |
| gene       | string | KRAS    | the gene name                                                    |
| kbCategory | string |         | the graphkb copy variant vocabulary term this variant belongs to |


### Optional Columns

| Column         | Type    | Example | Description                                                                          |
| -------------- | ------- | ------- | ------------------------------------------------------------------------------------ |
| copyChange     | integer | -2      | the ploidy corrected copy change                                                     |
| lohState       | string  | HET     | the loss-of-heterozygosity category for this gene region                             |
| chromosomeBand | string  | X:p12.2 |                                                                                      |
| start          | integer |         | the genomic start position of the copy segment this gene copy number was called from |
| end            | integer |         | the genomic end position of the copy segment this gene copy number was called from   |
| cna            | float   | 1.22    | The copy number alteration (CNA) ratio                                               |
| log2cna        | float   |         |                                                                                      |


## Structural Variant (Fusion) Data

### Required Columns

| Column          | Type    | Example               | Description                                                        |
| --------------- | ------- | --------------------- | ------------------------------------------------------------------ |
| eventType       | string  | deletion              | the type of underlying structural variant                          |
| breakpoint      | string  | 12:123456\|14:1244662 | description of the breakpoints involved in this structural variant |
| gene1           | string  | EWSR1                 | the 5' (n-terminal) gene name                                      |
| gene2           | string  | FLI1                  | the 3' (c-terminal) gene name                                      |
| exon1           | integer | 1                     | the 5' (n-terminal) exon                                           |
| exon2           | integer | 2                     | the 3' (c-terminal) exon                                           |
| nTermTranscript | string  | ENST0001.2            | the 3' transcript name                                             |
| ctermTranscript | string  | ENST004.5             | the 5' transcript name                                             |

### Optional Columns


| Column           | Type    | Example | Description                                                             |
| ---------------- | ------- | ------- | ----------------------------------------------------------------------- |
| detectedIn       | string  | DNA     | the sample type(s) this SV was detected in                              |
| conventionalName | string  |         | cytogenic descriptor                                                    |
| svg              | string  |         | svg image file for this SV                                              |
| svgTitle         | string  |         | title for the svg                                                       |
| omicSupport      | boolean |         | flag to indicate this SV has support from both genome and transcriptome |
