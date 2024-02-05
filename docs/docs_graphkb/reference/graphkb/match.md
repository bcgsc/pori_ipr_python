# graphkb/match

Functions which return Variants from GraphKB which match some input variant definition

## FEATURES_CACHE

```python
FEATURES_CACHE: Set[str] = set()
```

## get\_equivalent\_features()

Match an equivalent list of features given some input feature name (or ID).

```python
def get_equivalent_features(
    conn: GraphKBConnection,
    gene_name: str,
    ignore_cache: bool = False,
    is_source_id: bool = False,
    source: str = '',
    source_id_version: str = '',
) -> List[Ontology]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection))
- gene_name (`str`): the gene name to search features by
- ignore_cache (`bool`): bypass the cache to always force a new request
- is_source_id (`bool`): treat the gene_name as the gene ID from the source database (ex. ENSG001)
- source (`str`): the name of the source database the gene definition is from (ex. ensembl)
- source_id_version (`str`): the version of the source_id

**Returns**

- List\[[Ontology](../types/#class-ontology)\]: equivalent feature records

**Examples**

```python
get_equivalent_features(conn, 'KRAS')
```

```python
get_equivalent_features(conn, 'ENSG001', source='ensembl', is_source_id=True)
```

```python
get_equivalent_features(conn, 'ENSG001', source='ensembl', source_id_version='1')
```

```python
get_equivalent_features(conn, '#3:44')
```


## cache\_missing\_features()

Create a cache of features that exist to avoid repeatedly querying
for missing features

```python
def cache_missing_features(conn: GraphKBConnection) -> None:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection))

## match\_category\_variant()

Returns a list of variants matching the input variant

```python
def match_category_variant(
    conn: GraphKBConnection,
    gene_name: str,
    category: str,
    root_exclude_term: str = '',
    gene_source: str = '',
    gene_is_source_id: bool = False,
    ignore_cache: bool = False,
) -> List[Variant]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection)): the graphkb connection object
- gene_name (`str`): the name of the gene the variant is in reference to
- category (`str`): the variant category (ex. copy loss)
- root_exclude_term (`str`)
- gene_source (`str`): The source database the gene is defined by (ex. ensembl)
- gene_is_source_id (`bool`): Indicates the gene name(s) input should be treated as sourceIds not names
- ignore_cache (`bool`)

**Returns**

- List\[[Variant](../types/#class-variant)\]: List of variant records from GraphKB which match the input

**Raises**

- [FeatureNotFoundError](../util/#class-featurenotfounderror): The gene could not be found in GraphKB

## match\_copy\_variant()

Returns a list of variants matching the input variant

```python
def match_copy_variant(
    conn: GraphKBConnection, gene_name: str, category: str, drop_homozygous: bool = False, **kwargs
) -> List[Variant]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection)): the graphkb connection object
- gene_name (`str`): the name of the gene the variant is in reference to
- category (`str`): the variant category (ex. copy loss)
- drop_homozygous (`bool`): Drop homozygous matches from the result when true

**Returns**

- List\[[Variant](../types/#class-variant)\]: List of variant records from GraphKB which match the input

**Raises**

- `ValueError`: The input copy category is not recognized


## positions\_overlap()

Check if 2 Position records from GraphKB indicate an overlap

```python
def positions_overlap(
    pos_record: BasicPosition, range_start: BasicPosition, range_end: Optional[BasicPosition] = None
) -> bool:
```

**Args**

- pos_record ([BasicPosition](../types/#class-basicposition)): the record to compare
- range_start ([BasicPosition](../types/#class-basicposition)): the position record indicating the start of an uncertainty range
- range_end (Optional\[[BasicPosition](../types/#class-basicposition)\]): the position record indicating the end of an uncertainty range

**Returns**

- `bool`: True if the positions overlap

**Raises**

- `NotImplementedError`: if a cytoband type position is given

!!! note
	null values indicate not-specified or any

## compare\_positional\_variants()

Compare 2 variant records from GraphKB to determine if they are equivalent

```python
def compare_positional_variants(
    variant: Union[PositionalVariant, ParsedVariant],
    reference_variant: Union[PositionalVariant, ParsedVariant],
) -> bool:
```

**Args**

- variant (Union\[[PositionalVariant](../types/#class-positionalvariant), [ParsedVariant](../types/#class-parsedvariant)\]): the input variant
- reference_variant (Union\[[PositionalVariant](../types/#class-positionalvariant), [ParsedVariant](../types/#class-parsedvariant)\]): the reference (matched) variant record

**Returns**

- `bool`: True if the records are equivalent

## match\_positional\_variant()

Given the HGVS+ representation of some positional variant, parse it and match it to
annotations in GraphKB

```python
def match_positional_variant(
    conn: GraphKBConnection,
    variant_string: str,
    reference1: Optional[str] = None,
    reference2: Optional[str] = None,
    gene_is_source_id: bool = False,
    gene_source: str = '',
    ignore_cache: bool = False,
) -> List[Variant]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection))
- variant_string (`str`): the HGVS+ annotation string
- reference1 (`Optional[str]`): Explicitly specify the first reference link record (gene1)
- reference2 (`Optional[str]`): Explicitly specify the second reference link record (gene2)
- gene_is_source_id (`bool`): Indicates the gene name(s) input should be treated as sourceIds not names
- gene_source (`str`): The source database the gene is defined by (ex. ensembl)
- ignore_cache (`bool`)

**Returns**

- List\[[Variant](../types/#class-variant)\]: A list of matched statement records

**Raises**

- `NotImplementedError`: thrown for uncertain position input (ranges)
- [FeatureNotFoundError](../util/#class-featurenotfounderror): One of the genes does not exist in GraphKB
- `ValueError`: the gene names were given both in the variant_string and explicitly

**Examples**

```python
match_positional_variant(conn, '(EWSR1,FLI1):fusion(e.1,e.2)')
```

```python
match_positional_variant(conn, 'fusion(e.1,e.2)', 'EWSR1', 'FLI1')
```

```python
match_positional_variant(conn, 'fusion(e.1,e.2)', '#3:4', '#4:5')
```

```python
match_positional_variant(conn, 'fusion(e.1,e.2)', '123', '456', gene_is_source_id=True, gene_source='entrez gene')
```

```python
match_positional_variant(conn, 'KRAS:p.G12D')
```

```python
match_positional_variant(conn, 'p.G12D', 'KRAS')
```

