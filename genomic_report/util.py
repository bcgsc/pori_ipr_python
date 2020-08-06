import hashlib
import json
import logging
from typing import Dict, List, Set, Tuple

from graphkb.types import Record, Ontology
from graphkb.vocab import get_term_tree
from graphkb import GraphKBConnection
from graphkb.match import get_equivalent_features


from .types import IprGene

# name the logger after the package to make it simple to disable for packages using this one as a dependency
# https://stackoverflow.com/questions/11029717/how-do-i-disable-log-messages-from-the-requests-library
VERBOSE_ERROR_CODE = (logging.INFO + logging.DEBUG) // 2
logging.addLevelName(VERBOSE_ERROR_CODE, 'VERBOSE')
logger = logging.getLogger('genomic_report')
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


def get_genes_from_defn(conn: GraphKBConnection, gene_defn: IprGene) -> List[Ontology]:
    """
    Get the GraphKB record for a gene from some input gene defintion
    """
    if gene_defn.get('source') and gene_defn.get('sourceId'):
        if gene_defn.get('sourceIdVersion'):
            features = conn.query(
                {
                    'target': 'Feature',
                    'filters': [
                        {
                            'source': {
                                'target': 'Source',
                                'filters': {'name': gene_defn.get('source')},
                            }
                        },
                        {'sourceId': gene_defn.get('sourceId')},
                        {
                            'OR': [
                                {'sourceIdVersion': gene_defn.get('sourceIdVersion')},
                                {'sourceIdVersion': None},
                            ]
                        },
                    ],
                }
            )
        else:
            features = conn.query(
                {
                    'target': 'Feature',
                    'filters': [
                        {
                            'source': {
                                'target': 'Source',
                                'filters': {'name': gene_defn.get('source')},
                            }
                        },
                        {'sourceId': gene_defn.get('sourceId')},
                    ],
                }
            )
        if gene_defn.get('name') and len(
            [f for f in features if f['name'] == gene_defn.get('name')]
        ):
            return sorted(
                [f for f in features if f['name'] == gene_defn.get('name')],
                key=generate_ontology_preference_key,
            )
        return sorted(features, key=generate_ontology_preference_key)
    else:
        return sorted(
            get_equivalent_features(conn, gene_defn['name']), key=generate_ontology_preference_key
        )
