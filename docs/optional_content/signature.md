# Mutation Signatures

These are the scores of individual mutation signatures. These can include cosmic and/or custom
signatures. The nnls field is the non-negative least squares contribution of the given signature

```json
{
    "mutationSignature": [
        {
            "signature": "SBS1",
            "nnls": 0.344,
            "associations": "Tobacco chewing",
            "features": "D,T",
            "numCancerTypes": 1,
            "cancerTypes": "stomach cancer",
            "selected": false,
            "kbCategory": "strong signature"
        }
    ]
}
```

The `selected` field indicates if this signature should be shown on the front page of the report or
not
