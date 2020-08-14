# Pairwise RNA Expression Correlation

Provide a list of the most similar other samples with respect to the RNA expression profile.


```json
{
     "pairwiseExpressionCorrelation": [
        {
            "patientId": "UPLOADPAT02",
            "library": "LIB0002",
            "correlation": 0.99,
            "tumourType": "pancreatic cancer",
            "tissueType": "liver",
            "tumourContent": 15
        }
    ]
}
```

All values expect the correlation are attributes of the sample being compared to and not the
sample the report is being generated for
