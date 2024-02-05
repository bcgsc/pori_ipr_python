# graphkb/genes

Methods for retrieving gene annotation lists from GraphKB


## get\_oncokb\_oncogenes()

Gets the list of oncogenes stored in GraphKB derived from OncoKB.

```python
def get_oncokb_oncogenes(conn: GraphKBConnection) -> List[Ontology]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection)): the graphkb connection object

**Returns**

- List\[[Ontology](../types/#class-ontology)\]: gene (Feature) records

## get\_oncokb\_tumour\_supressors()

Gets the list of tumour supressor genes stored in GraphKB derived from OncoKB.

```python
def get_oncokb_tumour_supressors(conn: GraphKBConnection) -> List[Ontology]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection)): the graphkb connection object

**Returns**

- List\[[Ontology](../types/#class-ontology)\]: gene (Feature) records

## get\_therapeutic\_associated\_genes()

Genes related to a cancer-associated statement in Graphkb.

```python
def get_therapeutic_associated_genes(graphkb_conn: GraphKBConnection) -> List[Ontology]:
```

**Args**

- graphkb_conn ([GraphKBConnection](../util/#class-graphkbconnection))

**Returns**

- List\[[Ontology](../types/#class-ontology)\]

## get\_genes\_from\_variant\_types()

Retrieve a list of Genes which are found in variants on the given types

```python
def get_genes_from_variant_types(
    conn: GraphKBConnection,
    types: List[str],
    source_record_ids: List[str] = [],
    ignore_cache: bool = False,
) -> List[Ontology]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection)): the graphkb connection object
- types (`List[str]`): list of names of variant types
- source_record_ids (`List[str]`): list of sources ids to filter genes by
- ignore_cache (`bool`)

**Returns**

- List\[[Ontology](../types/#class-ontology)\]: gene (Feature) records

## get\_preferred\_gene\_name()

Preferred gene symbol of a gene or transcript.

```python
def get_preferred_gene_name(
    conn: GraphKBConnection, gene_name: str, source: str = PREFERRED_GENE_SOURCE
) -> str:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection))
- gene_name (`str`): the gene name to search features by
- source (`str`): id of the preferred gene symbol source

**Returns**

- `str`: preferred displayName symbol.

**Examples**

```python
return KRAS for get_preferred_gene_name(conn, 'NM_033360')
return KRAS for get_preferred_gene_name(conn, 'ENSG00000133703.11')
```


## get\_cancer\_predisposition\_info()

Return two lists from GraphKB, one of cancer predisposition genes and one of associated variants.

GERO-272 - criteria for what counts as a "cancer predisposition" variant

In short:
* Statement 'source' is 'CGL'
* Statement 'relevance' is 'pathogenic'
* gene is gotten from any associated 'PositionalVariant' records

Example: https://graphkb.bcgsc.ca/view/Statement/155:11616

```python
def get_cancer_predisposition_info(conn: GraphKBConnection) -> Tuple[List[str], Dict[str, str]]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection))

**Returns**

- `Tuple[List[str], Dict[str, str]]`: list of cancer predisposition genes variants: dictionary mapping pharmacogenomic variant IDs to variant display names

## get\_pharmacogenomic\_info()

Return two lists from GraphKB, one of pharmacogenomic genes and one of associated variants.

SDEV-2733 - criteria for what counts as a "pharmacogenomic" variant

In short:
* Statement 'source' is not 'CGI' or 'CIViC'
* Statement 'relevance' is 'increased toxicity' or 'decreased toxicity'
* gene is gotten from any associated 'PositionalVariant' records

Example: https://graphkb.bcgsc.ca/view/Statement/154:9574

```python
def get_pharmacogenomic_info(conn: GraphKBConnection) -> Tuple[List[str], Dict[str, str]]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection))

**Returns**

- `Tuple[List[str], Dict[str, str]]`: list of pharmacogenomic genes variants: dictionary mapping pharmacogenomic variant IDs to variant display names
