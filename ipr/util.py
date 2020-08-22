import hashlib
import json
import logging
from typing import Dict, List, Set, Tuple

from graphkb.types import Record
from graphkb.vocab import get_term_tree
from graphkb import GraphKBConnection

from .types import IprVariant


# name the logger after the package to make it simple to disable for packages using this one as a dependency
# https://stackoverflow.com/questions/11029717/how-do-i-disable-log-messages-from-the-requests-library
VERBOSE_ERROR_CODE = (logging.INFO + logging.DEBUG) // 2
logging.addLevelName(VERBOSE_ERROR_CODE, 'VERBOSE')
logger = logging.getLogger('ipr')
# add shortbut for verbose logging
setattr(logger, 'verbose', lambda *pos, **kw: logger.log(VERBOSE_ERROR_CODE, *pos, **kw))
LOG_LEVELS = {
    'info': logging.INFO,
    'debug': logging.DEBUG,
    'warn': logging.WARN,
    'error': logging.ERROR,
    'verbose': VERBOSE_ERROR_CODE,
}


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


def convert_to_rid_set(records: List[Record]) -> Set[str]:
    return {r['@rid'] for r in records}


def trim_empty_values(obj: Dict, empty_values: List = ['', None]) -> Dict:
    blacklist = ['gene1', 'gene2']  # allow null for sv genes
    keys = list(obj.keys())

    for key in keys:
        if obj[key] in empty_values and key not in blacklist:
            del obj[key]
    return obj


def create_variant_name(variant: IprVariant) -> str:
    """
    Given an IPR variant row, create the variant representation to be used as the name
    of the variant
    """
    variant_type = variant['variantType']
    if variant_type == 'exp':
        gene = variant['gene']
        return f'{gene} ({variant["expressionState"]})'
    elif variant_type == 'cnv':
        gene = variant['gene']
        return f'{gene} ({variant["cnvState"]})'
    return variant['variant']


def create_variant_name_tuple(variant: IprVariant) -> Tuple[str, str]:
    """
    Given an IPR variant row, create the variant representation to be used as the name
    of the variant
    """
    variant_type = variant['variantType']
    gene = variant['gene'] if 'gene' in variant else variant['gene1']
    if variant_type == 'exp':
        gene = variant['gene']
        return (gene, variant['expressionState'])
    elif variant_type == 'cnv':
        gene = variant['gene']
        return (gene, variant['cnvState'])
    variant_split = variant['variant'].split(':', 1)[1]
    gene2 = variant.get('gene2')
    if gene and gene2:
        gene = f'{gene}, {gene2}'
    elif gene2:
        gene = gene2

    return (gene, variant_split)


def find_variant(all_variants: List[IprVariant], variant_type: str, variant_key: str) -> IprVariant:
    """
    Find a variant in a list of variants by its key and type
    """
    for variant in all_variants:
        if variant['key'] == variant_key and variant['variantType'] == variant_type:
            return variant
    raise KeyError(f'expected variant ({variant_key}, {variant_type}) does not exist')


def generate_ontology_preference_key(record: Dict, sources_sort: Dict[str, int] = {}) -> Tuple:
    """
    Generate a tuple key for comparing preferred ontology terms.
    """
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


def get_alternatives(graphkb_conn: GraphKBConnection, record_id: str) -> List[Dict]:
    return graphkb_conn.query({'target': [record_id], 'queryType': 'similarTo', 'treeEdges': []})


def get_preferred_drug_representation(graphkb_conn: GraphKBConnection, drug_record_id: str) -> Dict:
    """
    Given a Drug record, follow its linked records to find the preferred
    representation by following alias, deprecating, and cross reference links
    """
    source_preference = {
        r['@rid']: r['sort']
        for r in graphkb_conn.query({'target': 'Source', 'returnProperties': ['sort', '@rid']})
    }
    drugs = sorted(
        get_alternatives(graphkb_conn, drug_record_id),
        key=lambda rec: generate_ontology_preference_key(rec, source_preference),
    )
    return drugs[0]


def get_preferred_gene_name(graphkb_conn: GraphKBConnection, record_id: str) -> str:
    """
    Given some Feature record ID return the preferred gene name
    """
    record = graphkb_conn.get_record_by_id(record_id)
    biotype = record.get('biotype', '')
    genes = []
    expanded = graphkb_conn.query({'target': [record_id], 'neighbors': 3})[0]

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
    return record['displayName']
