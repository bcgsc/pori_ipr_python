# graphkb/constants

## DEFAULT_LIMIT

```python
DEFAULT_LIMIT = 1000
```

## GKB_BASE_URL

```python
GKB_BASE_URL = "https://graphkb-api.bcgsc.ca/api"
```

## GKB_STAGING_URL

```python
GKB_STAGING_URL = "https://graphkbstaging-api.bcgsc.ca/api"
```

## GKB_DEV_URL

```python
GKB_DEV_URL = "https://graphkbdev-api.bcgsc.ca/api"
```

## DEFAULT_URL

```python
DEFAULT_URL = GKB_BASE_URL
```

## PREFERRED_GENE_SOURCE

```python
PREFERRED_GENE_SOURCE = "#39:5"  # HGNC
```

## BASE_RETURN_PROPERTIES

```python
BASE_RETURN_PROPERTIES = ['@rid', '@class']
```

## GENERIC_RETURN_PROPERTIES

```python
GENERIC_RETURN_PROPERTIES = [
    'name',
    'sourceId',
    'sourceIdVersion',
    'source.name',
    'source.@rid',
    'displayName',
    'deprecated',
] + BASE_RETURN_PROPERTIES
```

## GENE_RETURN_PROPERTIES

```python
GENE_RETURN_PROPERTIES = ['biotype'] + GENERIC_RETURN_PROPERTIES
```

## VARIANT_RETURN_PROPERTIES

```python
VARIANT_RETURN_PROPERTIES = (
    BASE_RETURN_PROPERTIES
    + [f'type.{p}' for p in GENERIC_RETURN_PROPERTIES]
    + [f'reference1.{p}' for p in GENE_RETURN_PROPERTIES]
    + [f'reference2.{p}' for p in GENE_RETURN_PROPERTIES]
    + ['zygosity', 'germline', 'displayName']
```

## POS_VARIANT_RETURN_PROPERTIES

```python
POS_VARIANT_RETURN_PROPERTIES = VARIANT_RETURN_PROPERTIES + [
    'break1Start',
    'break1End',
    'break2Start',
    'break2End',
    'break1Repr',
    'break2Repr',
    'refSeq',
    'untemplatedSeq',
    'untemplatedSeqSize',
    'truncation',
    'assembly',
```

## STATEMENT_RETURN_PROPERTIES

```python
STATEMENT_RETURN_PROPERTIES = (
    BASE_RETURN_PROPERTIES
    + ['displayNameTemplate', 'sourceId', 'source.name', 'source.displayName']
    + [f'conditions.{p}' for p in GENERIC_RETURN_PROPERTIES]
    + [f'subject.{p}' for p in GENERIC_RETURN_PROPERTIES]
    + [f'evidence.{p}' for p in GENERIC_RETURN_PROPERTIES]
    + [f'relevance.{p}' for p in GENERIC_RETURN_PROPERTIES]
    + [f'evidenceLevel.{p}' for p in GENERIC_RETURN_PROPERTIES]
    + ['reviewStatus']
```

## ONCOKB_SOURCE_NAME

```python
ONCOKB_SOURCE_NAME = 'oncokb'
```

## ONCOGENE

```python
ONCOGENE = 'oncogenic'
```

## TUMOUR_SUPPRESSIVE

```python
TUMOUR_SUPPRESSIVE = 'tumour suppressive'
```

## FUSION_NAMES

```python
FUSION_NAMES = ['structural variant', 'fusion']
```

## PHARMACOGENOMIC_SOURCE_EXCLUDE_LIST

```python
PHARMACOGENOMIC_SOURCE_EXCLUDE_LIST = ["cancer genome interpreter", "civic"]
```

## BASE_THERAPEUTIC_TERMS

```python
BASE_THERAPEUTIC_TERMS = ['therapeutic efficacy', 'eligibility']
```

## RELEVANCE_BASE_TERMS

```python
RELEVANCE_BASE_TERMS: CategoryBaseTermMapping = [
    ('therapeutic', BASE_THERAPEUTIC_TERMS),
    ('diagnostic', ['diagnostic indicator']),
    ('prognostic', ['prognostic indicator']),
    ('pharmacogenomic', ['metabolism', 'toxicity', 'dosage']),
    ('cancer predisposition', ['pathogenic']),
    ('biological', ['functional effect', 'tumourigenesis', 'predisposing']),
]
```

## FAILED_REVIEW_STATUS

```python
FAILED_REVIEW_STATUS = 'failed'
```

## CHROMOSOMES_HG38

```python
CHROMOSOMES_HG38 = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
```

## CHROMOSOMES_HG19

```python
CHROMOSOMES_HG19 = [str(i) for i in range(1, 23)] + ['x', 'y', 'mt']
```

## CHROMOSOMES

```python
CHROMOSOMES = CHROMOSOMES_HG38 + CHROMOSOMES_HG19
```

## AMBIGUOUS_AA

```python
AMBIGUOUS_AA = ['x', '?', 'X']
```

## AA_3to1_MAPPING

```python
AA_3to1_MAPPING = {
    'Ala': 'A',
    'Arg': 'R',
    'Asn': 'N',
    'Asp': 'D',
    'Asx': 'B',
    'Cys': 'C',
    'Glu': 'E',
    'Gln': 'Q',
    'Glx': 'Z',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Leu': 'L',
    'Lys': 'K',
    'Met': 'M',
    'Phe': 'F',
    'Pro': 'P',
    'Ser': 'S',
    'Thr': 'T',
    'Trp': 'W',
    'Tyr': 'Y',
    'Val': 'V',
    'Ter': '*',
```

## INPUT_COPY_CATEGORIES

```python
INPUT_COPY_CATEGORIES = IterableNamespace(
    AMP='amplification',
    ANY_GAIN='copy gain',
    ANY_LOSS='copy loss',
    DEEP='deep deletion',
    GAIN='low level copy gain',
    LOSS='shallow deletion',
)
```

## INPUT_EXPRESSION_CATEGORIES

```python
INPUT_EXPRESSION_CATEGORIES = IterableNamespace(
    UP='increased expression', DOWN='reduced expression'
)
```

## class IterableNamespace

**inherits** `argparse.Namespace`






