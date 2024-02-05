# graphkb/types

Type annotations used for static type checking in this module

## Record

```python
Record: TypedDict = TypedDict('Record', {'@rid': str, '@class': str})
```

**Attributes**

- @rid (`str`)
- @class (`str`)

## EmbeddedRecord

```python
EmbeddedRecord: TypedDict = TypedDict('EmbeddedRecord', {'@class': str})
```

**Attributes**

- @class (`str`)

## RecordLink

```python
RecordLink = Union[str, Record]
```

## OntologyLink

```python
OntologyLink = Union[str, Ontology]
```

## Position

```python
Position = Union[BasicPosition, CytobandPosition]
```

## CategoryBaseTermMapping

```python
CategoryBaseTermMapping = List[Tuple[str, List[str]]]
```

## class Ontology

**inherits** [Record](#record)

**Attributes**

- sourceId (`str`)
- name (`str`)
- source ([RecordLink](#recordlink))
- displayName (`str`)

## class BasicPosition

**inherits** [EmbeddedRecord](#embeddedrecord)

**Attributes**

- pos (`int`)

## class CytobandPosition

**inherits** [EmbeddedRecord](#embeddedrecord)

**Attributes**

- arm (`str`)
- majorBand (`str`)
- minorBand (`str`)

## class Variant

**inherits** [Record](#record)

**Attributes**

- reference1 ([OntologyLink](#ontologylink))
- reference2 (Optional\[[OntologyLink](#ontologylink)\])
- type ([OntologyLink](#ontologylink))
- zygosity (`str`)
- germline (`bool`)
- displayName (`str`)

## class PositionalVariant

**inherits** [Variant](#class-variant)

**Attributes**

- break1Start (Union\[[Position](#position), [CytobandPosition](#class-cytobandposition)\])
- break1End (`Optional[Union]`)
- break2Start (`Optional[Union]`)
- break2End (`Optional[Union]`)
- refSeq (`Optional[str]`)
- untemplatedSeq (`Optional[str]`)
- untemplatedSeqSize (`Optional[int]`)

## class ParsedVariant

**inherits** `TypedDict`

**Attributes**

- reference1 (`str`)
- reference2 (`Optional[str]`)
- type (`str`)
- zygosity (`str`)
- germline (`bool`)
- break1Start (Union\[[Position](#position), [CytobandPosition](#class-cytobandposition)\])
- break1End (`Optional[Union]`)
- break2Start (`Optional[Union]`)
- break2End (`Optional[Union]`)
- refSeq (`Optional[str]`)
- untemplatedSeq (`Optional[str]`)
- untemplatedSeqSize (`Optional[int]`)

## class Statement

**inherits** [Record](#record)

**Attributes**

- relevance ([OntologyLink](#ontologylink))
- subject ([OntologyLink](#ontologylink))
- conditions (List\[[OntologyLink](#ontologylink)\])
- evidence (List\[[OntologyLink](#ontologylink)\])
- evidenceLevel (List\[[OntologyLink](#ontologylink)\])
- source ([RecordLink](#recordlink))
- sourceId (`str`)
