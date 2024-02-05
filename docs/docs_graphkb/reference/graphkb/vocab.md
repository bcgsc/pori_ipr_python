# graphkb/vocab


## get\_equivalent\_terms()

Get a list of terms equivalent to the current term up to the root term

```python
def get_equivalent_terms(
    conn: GraphKBConnection,
    base_term_name: str,
    root_exclude_term: str = '',
    ontology_class: str = 'Vocabulary',
    ignore_cache: bool = False,
    build_base_query: Callable = query_by_name,
) -> List[Ontology]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection))
- base_term_name (`str`): the name to get superclasses of
- root_exclude_term (`str`): the parent term to exlcude along with all of its parent terms
- ontology_class (`str`)
- ignore_cache (`bool`)
- build_base_query (`Callable`)

**Returns**

- List\[[Ontology](../types/#class-ontology)\]

## get\_term\_tree()

Get terms equivalent to the base term by traversing the subclassOf tree and expanding related
alias and cross reference edges

```python
def get_term_tree(
    conn: GraphKBConnection,
    base_term_name: str,
    root_exclude_term: str = '',
    ontology_class: str = 'Vocabulary',
    include_superclasses: bool = True,
    ignore_cache: bool = False,
    build_base_query: Callable = query_by_name,
) -> List[Ontology]:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection)): the graphkb connection object
- base_term_name (`str`): the term to use as the base of the subclass tree
- root_exclude_term (`str`)
- ontology_class (`str`): the default class to query. Defaults to 'Vocabulary'
- include_superclasses (`bool`): when True the query will include superclasses of the current term
- ignore_cache (`bool`)
- build_base_query (`Callable`)

**Returns**

- List\[[Ontology](../types/#class-ontology)\]: GraphKB records Note: this must be done in 2 calls to avoid going up and down the tree in a single query (exclude adjacent siblings)

## get\_term\_by\_name()

Retrieve a vocaulary term by name

```python
def get_term_by_name(
    conn: GraphKBConnection,
    name: str,
    ontology_class: str = 'Vocabulary',
    ignore_cache: bool = False,
    **kwargs,
) -> Ontology:
```

**Args**

- conn ([GraphKBConnection](../util/#class-graphkbconnection)): the graphkb connection object
- name (`str`): the name of the Vocabulary term to retrieve
- ontology_class (`str`)
- ignore_cache (`bool`)

**Returns**

- [Ontology](../types/#class-ontology): Vocabulary record

**Raises**

- `AssertionError`: more than one term or no terms with that name were found
- `AssertionError`: if the term was not found or more than 1 match was found (expected to be unique)

## get\_terms\_set()

Get a set of vocabulary rids given some base/parent term names.

```python
def get_terms_set(
    graphkb_conn: GraphKBConnection, base_terms: Iterable[str], ignore_cache: bool = False
) -> Set[str]:
```

**Args**

- graphkb_conn ([GraphKBConnection](../util/#class-graphkbconnection))
- base_terms (`Iterable[str]`)
- ignore_cache (`bool`)

**Returns**

- `Set[str]`
