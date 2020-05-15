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


| Column        | Type    | Description                                                      |
| ------------- | ------- | ---------------------------------------------------------------- |
| rnaReads      | float   |                                                                  |
| rpkm          | float   |                                                                  |
| foldChange    | float   |                                                                  |
| tcgaPerc      | float   | percentile value when compared to the primary disease comparator |
| tcgaPercCol   | string  | primary disease comparator                                       |
| tcgakIQR      | float   | kIQR value when compared to the primary disease comparator       |
| tcgaQC        | float   |                                                                  |
| tcgaQCCol     | string  |                                                                  |
| tcgaAvgPerc   | float   |                                                                  |
| tcgaAvgkIQR   | float   |                                                                  |
| tcgaAvgQC     | float   |                                                                  |
| tcgaAvgQCCol  | string  |                                                                  |
| tcgaNormPerc  | float   |                                                                  |
| tcgaNormkIQR  | float   |                                                                  |
| ptxPerc       | float   |                                                                  |
| ptxkIQR       | float   |                                                                  |
| ptxQC         | float   |                                                                  |
| ptxPercCol    | string  |                                                                  |
| ptxTotSampObs | integer |                                                                  |
| ptxPogPerc    | float   |                                                                  |
| gtexComp      | string  |                                                                  |
| gtexPerc      | float   |                                                                  |
| gtexFC        | float   |                                                                  |
| gtexkIQR      | float   |                                                                  |
| gtexAvgPerc   | float   |                                                                  |
| gtexAvgFC     | float   |                                                                  |
| gtexAvgkIQR   | float   |                                                                  |

## Small Mutation Data

Small mutations are composed of indels and single nucleotide variants.

### Required Columns

| Column        | Type   | Example     | Description                          |
| ------------- | ------ | ----------- | ------------------------------------ |
| location      | string | 10:123456   | the genomic position of this variant |
| refAlt        | string | A>T         | the genomic sequence change          |
| gene          | string | KRAS        | the gene name                        |
| proteinChange | string | p.G12D      | the HGVS protein notation            |
| transcript    | string | ENST00001.2 | the transcript name                  |


### Optional Columns

| Column      | Type   | Example | Description                                                                                |
| ----------- | ------ | ------- | ------------------------------------------------------------------------------------------ |
| zygosity    | string | het     |                                                                                            |
| tumourReads | string | 4/8     | the reference and alternate (supporting) read counts at this position in the genome        |
| rnaReads    | string | 4/8     | the reference and alternate (supporting) read counts at this position in the transcriptome |
| detectedIn  | string | DNA/RNA | the sample types this variant was detected in                                              |

## Copy Variant Data


### Required Columns

| Column     | Type   | Example | Description                                                      |
| ---------- | ------ | ------- | ---------------------------------------------------------------- |
| gene       | string | KRAS    | the gene name                                                    |
| kbCategory | string |         | the graphkb copy variant vocabulary term this variant belongs to |


### Optional Columns

| Column              | Type    | Example | Description                                                                          |
| ------------------- | ------- | ------- | ------------------------------------------------------------------------------------ |
| ploidyCorrCpyChange | integer | -2      | the ploidy corrected copy change                                                     |
| lohState            | string  | HET     | the loss-of-heterozygosity category for this gene region                             |
| chromosomeBand      | string  |         |                                                                                      |
| start               | integer |         | the genomic start position of the copy segment this gene copy number was called from |
| end                 | integer |         | the genomic end position of the copy segment this gene copy number was called from   |


## Structural Variant (Fusion) Data

### Required Columns

| Column          | Type    | Example               | Description                                                        |
| --------------- | ------- | --------------------- | ------------------------------------------------------------------ |
| eventType       | string  | deletion              | the type of underlying structural variant                          |
|                 |
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
