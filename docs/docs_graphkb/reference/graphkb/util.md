# graphkb/util

## QUERY_CACHE

```python
QUERY_CACHE: Dict[Any, Any] = {}
```

## logger

```python
logger = logging.getLogger('graphkb')
```


## class GraphKBConnection



### GraphKBConnection.request()

Request wrapper to handle adding common headers and logging.

```python
def request(self, endpoint: str, method: str = 'GET', **kwargs) -> Dict:
```

**Args**

- endpoint (`str`): api endpoint, excluding the base uri
- method (`str`): the http method. Defaults to 'GET'.

**Returns**

- `Dict`: the json response as a python dict

### GraphKBConnection.post()

Convenience method for making post requests.

```python
def post(self, uri: str, data: Dict = {}, **kwargs) -> Dict:
```

**Args**

- uri (`str`)
- data (`Dict`)

**Returns**

- `Dict`



### GraphKBConnection.set\_cache\_data()

Explicitly add a query to the cache.

```python
def set_cache_data(self, request_body: Dict, result: List[Record]) -> None:
```

**Args**

- request_body (`Dict`)
- result (List\[[Record](../types/#record)\])

### GraphKBConnection.query()

Query GraphKB

```python
def query(
    self,
    request_body: Dict = {},
    paginate: bool = True,
    ignore_cache: bool = False,
    force_refresh: bool = False,
    limit: int = DEFAULT_LIMIT,
) -> List[Record]:
```

**Args**

- request_body (`Dict`)
- paginate (`bool`)
- ignore_cache (`bool`)
- force_refresh (`bool`)
- limit (`int`)

**Returns**

- List\[[Record](../types/#record)\]






## convert\_to\_rid\_list()

Given a list of records or record id strings, return their record IDs.

```python
def convert_to_rid_list(records: Iterable[Record]) -> List[str]:
```

**Args**

- records (Iterable\[[Record](../types/#record)\])

**Returns**

- `List[str]`

## looks\_like\_rid()

Check if an input string looks like a GraphKB ID.

```python
def looks_like_rid(rid: str) -> bool:
```

**Args**

- rid (`str`)

**Returns**

- `bool`

## convert\_aa\_3to1()

Convert an Input string from 3 letter AA notation to 1 letter AA notation.

```python
def convert_aa_3to1(three_letter_notation: str) -> str:
```

**Args**

- three_letter_notation (`str`)

**Returns**

- `str`

## join\_url()

Join parts of a URL into a full URL.

```python
def join_url(base_url: str, *parts) -> str:
```

**Args**

- base_url (`str`)

**Returns**

- `str`

## millis\_interval()

Millisecond time from start and end datetime instances.

```python
def millis_interval(start: datetime, end: datetime) -> int:
```

**Args**

- start (`datetime`)
- end (`datetime`)

**Returns**

- `int`

## cache\_key()

Create a cache key for a query request to GraphKB.

```python
def cache_key(request_body) -> str:
```

**Args**

- request_body

**Returns**

- `str`

## get\_rid()

Retrieve a record by name and target

```python
def get_rid(conn: GraphKBConnection, target: str, name: str) -> str:
```

**Args**

- conn ([GraphKBConnection](#class-graphkbconnection)): GraphKBConnection
- target (`str`): record type to query
- name (`str`): the name of the record to retrieve

**Returns**

- `str`: @rid of the record

**Raises**

- `AssertionError`: if the term was not found or more than 1 match was found (expected to be unique)
