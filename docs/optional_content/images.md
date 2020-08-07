# Images

There are a number of images that can be uploaded to the report. Images require a path to the image
and the key for the image. The image key is used to tell IPR how to place the image in the report.
The following are examples with their expected format and image key.

```json
{
    "images": [
        {
            "key": "string",
            "path": "/path/to/image/file"
        }
    ]
}
```

## Mutation Signature Plots

key: `mutSignature.barplot.sbs`

![mutation signature plot](../images/mutation_signature.sbs.png)

key: `mutSignature.barplot.dbs`

![mutation signature plot](../images/mutation_signature.dbs.png)

key: `mutSignature.barplot.indels`

![mutation signature plot](../images/mutation_signature.indels.png)


## Expression Correlation Subtyping Plots

key: `subtypePlot\.\S+`

## Copy Number Circos Plot

key: `cnvLoh.circos`

## Mutation Burden Plots

key: `mutation_summary\.(barplot|density|legend)_(sv|snv|indel)(\.\w+)?`

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
