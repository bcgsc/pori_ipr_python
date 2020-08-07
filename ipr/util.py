import hashlib
import json
import logging
from typing import Dict, List, Set, Tuple

from graphkb.types import Record
from graphkb.vocab import get_term_tree
from graphkb import GraphKBConnection

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
