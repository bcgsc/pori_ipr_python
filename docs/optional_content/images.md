# Images

There are a number of images that can be uploaded to the report. Images require a path to the image
and the key for the image. The image key is used to tell IPR how to place the image in the report.
The following are examples with their expected format and image key.

```json
{
    "images": [
        {
            "key": "string",
            "path": "/path/to/image/file",
            "title": "",
            "caption": ""
        }
    ]
}
```

## Mutation Signature Plots

key: `mutSignature.barplot.sbs`

![mutation signature plot](../images/mutSignature.barplot.sbs.png)

key: `mutSignature.barplot.dbs`

![mutation signature plot](../images/mutSignature.barplot.dbs.png)

key: `mutSignature.barplot.indels`

![mutation signature plot](../images/mutSignature.barplot.indels.png)


## Expression Correlation Subtyping Plots

key: `subtypePlot\.\S+`

## Copy Number Circos Plot

key: `cnvLoh.circos`

## Mutation Burden Plots

full pattern: `mutationBurden\.(barplot|density|legend)_(sv|snv|indel)\.(primary|secondary|tertiary|quaternary)`

### SNVs

key: `mutationBurden\.density_snv\.(primary|secondary|tertiary|quaternary)`

![snv density plot](../images/mutationBurden.density_snv.primary.png)

key: `mutationBurden\.barplot_snv\.(primary|secondary|tertiary|quaternary)`

![snv bar plot](../images/mutationBurden.barplot_snv.primary.png)

key: `mutationBurden\.density_indel\.(primary|secondary|tertiary|quaternary)`

### Indels

![indel density plot](../images/mutationBurden.density_indel.primary.png)

key: `mutationBurden\.barplot_indel\.(primary|secondary|tertiary|quaternary)`

![indel bar plot](../images/mutationBurden.barplot_indel.primary.png)

### Structural Variants

key: `mutationBurden\.density_sv\.(primary|secondary|tertiary|quaternary)`

![sv density plot](../images/mutationBurden.density_sv.primary.png)

## Expression Density Plots

key: `expDensity\.(\S+)`

In the above the pattern is expected to be `expDensity.<gene name>` where the gene name
matches the gene name(s) used for the expression variant definitions

![expression density plot](../images/expression_density.png)

## Expression Correlation Plot

keys: `expression.chart`, `expression.legend`

![expression correlation plot](../images/expression_correlation.png)

## Structural Variant Circos Plots

keys: `circosSv.genome`, `circosSv.transcriptome`

## Microbial Integration Circos Plots

keys: `microbial.circos.genome`, `microbial.circos.transcriptome`

## Immune Related Plots

### Cibersort

[Cibersort Publication](https://www.nature.com/articles/nmeth.3337)

key: `cibersort.cd8_positive_t-cell_scatter`

![cd8+ T-cell scatter plot](../images/cibersort.cd8_positive_t-cell_scatter.png)

key: `cibersort.combined_t-cell_scatter`

![combined T-cell scatter plot](../images/cibersort.combined_t-cell_scatter.png)

### MiXCR

[MiXCR Publication](https://pubmed.ncbi.nlm.nih.gov/25924071/)

key: `mixcr.circos_trb_vj_gene_usage`

![mixcr circos](../images/mixcr.circos_trb_vj_gene_usage.png)

key: `mixcr.dominance_vs_alpha_beta_t-cells_scatter`

![mixcr dominance](../images/mixcr.dominance_vs_alpha_beta_t-cells_scatter.png)
