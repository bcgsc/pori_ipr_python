import hashlib
import json
import logging
import pandas as pd
from graphkb import GraphKBConnection
from graphkb.types import Ontology, Record
from graphkb.vocab import get_term_tree
from numpy import nan
from typing import Any, Dict, List, Sequence, Set, Tuple, cast

from .types import IprVariant

GENE_NEIGHBORS_MAX = 3

# name the logger after the package to make it simple to disable for packages using this one as a dependency
logger = logging.getLogger('ipr')
LOG_LEVELS = {
    'info': logging.INFO,
    'debug': logging.DEBUG,
    'warn': logging.WARN,
    'error': logging.ERROR,
}


class Hashabledict(dict):
    def __hash__(self):
        return hash(frozenset(self))


def get_terms_set(graphkb_conn: GraphKBConnection, base_terms: List[str]) -> Set[str]:
    terms = set()
    for base_term in base_terms:
        terms.update(
            convert_to_rid_set(get_term_tree(graphkb_conn, base_term, include_superclasses=False))
        )
    return terms


def hash_key(key: Tuple[str]) -> str:
    body = json.dumps({'key': key}, sort_keys=True)
    hash_code = hashlib.md5(body.encode('utf-8')).hexdigest()
    return hash_code


def convert_to_rid_set(records: Sequence[Record]) -> Set[str]:
    return {r['@rid'] for r in records}


def trim_empty_values(obj: IprVariant, empty_values: Sequence = ('', None, nan)):
    blacklist = ('gene1', 'gene2')  # allow null for sv genes
    keys = list(obj.keys())

    for key in keys:
        if obj[key] in empty_values and key not in blacklist:  # type: ignore
            del obj[key]  # type: ignore
    return obj


def create_variant_name_tuple(variant: IprVariant) -> Tuple[str, str]:
    """
    Given an IPR variant row, create the variant representation to be used as the name
    of the variant
    """
    variant_type = variant['variantType']
    gene = str(variant.get('gene', variant.get('gene1', '')))
    if variant_type == 'exp':
        return (gene, str(variant.get('expressionState', '')))
    elif variant_type == 'cnv':
        return (gene, str(variant.get('cnvState', '')))
    variant_split = (
        variant['variant'].split(':', 1)[1] if ':' in variant['variant'] else variant['variant']
    )

    gene2 = str(variant.get('gene2', ''))
    if gene and gene2:
        gene = f'{gene}, {gene2}'
    elif gene2:
        gene = gene2

    return (gene, variant_split)


def find_variant(
    all_variants: Sequence[IprVariant], variant_type: str, variant_key: str
) -> IprVariant:
    """
    Find a variant in a list of variants by its key and type
    """
    for variant in all_variants:
        if variant['key'] == variant_key and variant['variantType'] == variant_type:
            return variant
    raise KeyError(f'expected variant ({variant_key}, {variant_type}) does not exist')


def generate_ontology_preference_key(record: Ontology, sources_sort: Dict[str, int] = {}) -> Tuple:
    """Generate a tuple key for comparing preferred ontology terms."""
    return (
        record.get('name') == record.get('sourceId'),
        record.get('deprecated', False),
        record.get('alias', False),
        bool(record.get('dependency', '')),
        sources_sort.get(record['source'], 99999),
        record['sourceId'],
        record.get('sourceIdVersion', ''),
        record['name'],
    )


def get_alternatives(graphkb_conn: GraphKBConnection, record_id: str) -> List[Ontology]:
    rec_list = graphkb_conn.query(
        {'target': [record_id], 'queryType': 'similarTo', 'treeEdges': []}
    )
    return [cast(Ontology, rec) for rec in rec_list]


def get_preferred_drug_representation(
    graphkb_conn: GraphKBConnection, drug_record_id: str
) -> Ontology:
    """Given a Drug record, follow its linked records to find the preferred
    representation by following alias, deprecating, and cross reference links.
    """

    source_preference = {
        r['@rid']: r['sort']
        for r in graphkb_conn.query({'target': 'Source', 'returnProperties': ['sort', '@rid']})
    }
    drugs = sorted(
        get_alternatives(graphkb_conn, drug_record_id),
        key=lambda rec: generate_ontology_preference_key(rec, source_preference),
    )
    return cast(Ontology, drugs[0])


def get_preferred_gene_name(
    graphkb_conn: GraphKBConnection, record_id: str, neighbors: int = GENE_NEIGHBORS_MAX
) -> str:
    """Given some Feature record ID return the preferred gene name."""
    record = graphkb_conn.get_record_by_id(record_id)
    biotype = record.get('biotype', '')
    genes = []
    expanded_gene_names = graphkb_conn.query({'target': [record_id], 'neighbors': neighbors})
    assert len(expanded_gene_names) == 1, "get_preferred_gene_name should have single result"
    expanded: Dict[str, List] = expanded_gene_names[0]  # type: ignore
    if biotype != 'gene':
        for edge in expanded.get('out_ElementOf', []):
            target = edge['in']
            if target.get('biotype') == 'gene':
                genes.append(target)

    for edge_type in [
        'out_AliasOf',
        'in_AliasOf',
        'in_DeprecatedBy',
        'out_CrossReferenceOf',
        'in_CrossReferenceOf',
    ]:
        target_name = 'out' if edge_type.startswith('in') else 'in'
        for edge in expanded.get(edge_type, []):
            target = edge[target_name]
            if target.get('biotype') == 'gene':
                genes.append(target)
    genes = sorted(
        genes,
        key=lambda gene: (
            gene['deprecated'],
            bool(gene['dependency']),
            '_' in gene['name'],
            gene['name'].startswith('ens'),
        ),
    )
    if genes:
        return genes[0]['displayName']
    # fallback to the input displayName
    return str(record.get('displayName', ''))


def pandas_falsy(field: Any) -> bool:
    """Check if a field is python falsy or pandas null."""
    return bool(pd.isnull(field) or not field)
