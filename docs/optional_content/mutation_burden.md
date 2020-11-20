# Mutation Burden

## Data

Measure of the relative counts of various types of mutations compared to those from a given cohort
of samples. The "role" determines with comparator these values are calculated in reference to

```json
{
    "mutationBurden": [
        {
            "snv": 4,
            "snvTruncating": 0,
            "indels": 0,
            "indelsFrameshift": 0,
            "sv": 92,
            "svExpressed": 40,
            "snvPercentile": 4,
            "indelPercentile": 10,
            "svPercentile": 60,
            "role": "primary"
        }
    ],
}
```

## Comparators and Roles

The role used by the mutation burden data is linked to the mutation burden images in the report
that is created by the comparators input.

The images use the following key: `mutationBurden\.(barplot|density|legend)_(sv|snv|indel)\.(primary|secondary|tertiary|quaternary)`

which maps to the comparators

- mutation burden (primary)
- mutation burden (secondary)
- mutation burden (tertiary)
- mutation burden (quaternary)

Therefore for a complete mutation burden entry one might see the following

```json
{
    "mutationBurden": [
        {
            "snv": 4,
            "snvTruncating": 0,
            "indels": 0,
            "indelsFrameshift": 0,
            "sv": 92,
            "svExpressed": 40,
            "snvPercentile": 4,
            "indelPercentile": 10,
            "svPercentile": 60,
            "role": "primary"
        }
    ],
    "comparators": [
        {
            "analysisRole": "mutation burden (primary)",
            "name": "TCGA COAD"
        }
    ],
    "images": [
        {
            "key": "mutationBurden.barplot_snv.primary",
            "path": "/path/to/image/file.png"
        },
        {
            "key": "mutationBurden.legend_snv.primary",
            "path": "/path/to/other/image/file.png"
        },
        {
            "key": "mutationBurden.density_snv.primary",
            "path": "/path/to/other-other/image/file.png"
        }
    ]
}
```

Note that the above only loaded images for the SNVs but similar images would be expected for SVs and
indels.

## SV-specific Comparator Overloading

There are special overload comparators for stuctural variants since they often have a different
comparator dataset from the SNV and Indels. Specifying an SV specific comparator will overload the
more general comparator for the SV types

- mutation burden SV (primary)
- mutation burden SV (secondary)
- mutation burden SV (tertiary)
- mutation burden SV (quaternary)

This means that when the report client lists the comparator that corresponds to the images and data
it will use the SV-specific value where it exists instead of the more general comparator name.
