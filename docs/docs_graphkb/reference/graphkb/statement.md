# graphkb/statement

## categorize\_relevance()

Given the record ID of some relevance term, return the higher level categorization

```python
def categorize_relevance(
    graphkb_conn: GraphKBConnection,
    relevance_rid: str,
    category_base_terms: CategoryBaseTermMapping = RELEVANCE_BASE_TERMS,
) -> str:
```

**Args**

- graphkb_conn ([GraphKBConnection](../util/#class-graphkbconnection))
- relevance_rid (`str`)
- category_base_terms ([CategoryBaseTermMapping](../types/#categorybasetermmapping))

**Returns**

- `str`

## get\_statements\_from\_variants()

Given a list of variant records from GraphKB, return related statements.

```python
def get_statements_from_variants(
    graphkb_conn: GraphKBConnection,
    variants: List[Variant],
    failed_review: bool = False,
) -> List[Statement]:
```

**Args**

- graphkb_conn ([GraphKBConnection](../util/#class-graphkbconnection)): the graphkb api connection object
- variants (List\[[Variant](../types/#class-variant)\]): list of variant records. (Have @rid property.)
- failed_review (`bool`): Include statements that failed review

**Returns**

- List\[[Statement](../types/#class-statement)\]: list of Statement records from graphkb
