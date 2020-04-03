"""
upload variant and report information to IPR
"""
from typing import List, Dict, Tuple
import requests
import json
import zlib

from graphkb import GraphKBConnection
from graphkb.util import IterableNamespace
from graphkb.vocab import get_term_tree

from .util import convert_to_rid_set

BASE_THERAPEUTIC_TERM = 'therapeutic efficacy'
BASE_DIAGNOSTIC_TERM = 'diagnostic indicator'
BASE_PROGNOSTIC_TERM = 'prognostic indicator'
BASE_BIOLOGICAL_TERM = 'functional effect'

VARIANT_CLASSES = {'Variant', 'CategoryVariant', 'PositionalVariant', 'CatalogueVariant'}
REPORT_KB_SECTIONS = IterableNamespace(
    therapeutic='therapeutic',
    prognostic='prognostic',
    biological='biological',
    unknown='unknown',
    novel='novel',
    diagnostic='diagnostic',
)

APPROVED_EVIDENCE_LEVELS = {
    # sourceIds of levels by source name
    'oncokb': ['1', 'r1'],
    'profyle': ['t1'],
    'cancer genome interpreter': [
        'cpic guidelines',
        'european leukemianet guidelines',
        'fda guidelines',
        'nccn guidelines',
        'nccn/cap guidelines',
    ],
}

DEFAULT_URL = 'https://iprdev-api.bcgsc.ca/api'
DEFAULT_LIMIT = 1000


def display_evidence_levels(statement: Dict) -> str:
    result = []

    for evidence_level in statement.get('evidenceLevel', []) or []:
        result.append(evidence_level['displayName'])

    return ';'.join(sorted(result))


def get_approved_evidence_levels(graphkb_conn: GraphKBConnection) -> List[Dict]:
    filters = []
    for source, names in APPROVED_EVIDENCE_LEVELS.items():
        filters.append(
            {
                'AND': [
                    {'source': {'target': 'Source', 'filters': {'name': source}}},
                    {'name': names, 'operator': 'IN'},
                ]
            }
        )
    return graphkb_conn.query({'target': 'EvidenceLevel', 'filters': {'OR': filters}})


def convert_statements_to_alterations(
    graphkb_conn: GraphKBConnection, statements: List[Dict], disease_name: str
) -> List[Dict]:
    """
    Given a set of statements matched from graphkb, convert these into their IPR equivalent representations

    Args:
        graphkb_conn: the graphkb connection object
        statements: list of statement records from graphkb
        disease_name: name of the cancer type for the patient being reported on

    Raises:
        ValueError: could not find the disease type in GraphKB

    Returns:
        IPR graphkb row representations
    """
    disease_matches = {
        r['@rid'] for r in get_term_tree(graphkb_conn, disease_name, ontology_class='Disease')
    }

    if not disease_matches:
        raise ValueError(f'failed to match disease ({disease_name}) to graphkb')

    rows = []

    approved = convert_to_rid_set(get_approved_evidence_levels(graphkb_conn))

    therapeutic_terms = convert_to_rid_set(
        get_term_tree(graphkb_conn, BASE_THERAPEUTIC_TERM, include_superclasses=False)
    )
    diagnostic_terms = convert_to_rid_set(
        get_term_tree(graphkb_conn, BASE_DIAGNOSTIC_TERM, include_superclasses=False)
    )
    prognostic_terms = convert_to_rid_set(
        get_term_tree(graphkb_conn, BASE_PROGNOSTIC_TERM, include_superclasses=False)
    )
    biological_terms = convert_to_rid_set(
        get_term_tree(graphkb_conn, BASE_BIOLOGICAL_TERM, include_superclasses=False)
    )

    for statement in statements:
        variants = [c for c in statement['conditions'] if c['@class'] in VARIANT_CLASSES]
        diseases = [c for c in statement['conditions'] if c['@class'] == 'Disease']
        pmid = ';'.join([e['displayName'] for e in statement['evidence']])

        ipr_section = REPORT_KB_SECTIONS.unknown  # table this goes to in IPR
        relevance_id = statement['relevance']['@rid']

        approved_therapy = False

        if relevance_id in therapeutic_terms:
            ipr_section = REPORT_KB_SECTIONS.therapeutic

            for level in statement['evidenceLevel'] or []:
                if level['@rid'] in approved:
                    approved_therapy = True
                    break

        elif relevance_id in diagnostic_terms:
            ipr_section = REPORT_KB_SECTIONS.diagnostic
        elif relevance_id in prognostic_terms:
            ipr_section = REPORT_KB_SECTIONS.prognostic
        elif relevance_id in biological_terms:
            ipr_section = REPORT_KB_SECTIONS.biological

        disease_match = len(diseases) == 1 and diseases[0]['@rid'] in disease_matches

        for variant in variants:
            row = {
                'approvedTherapy': approved_therapy,
                'category': ipr_section,
                'context': (statement['subject']['displayName'] if statement['subject'] else None),
                'disease': ';'.join(sorted(d['displayName'] for d in diseases)),
                'evidenceLevel': display_evidence_levels(statement),
                'kbStatementId': statement['@rid'],
                'kbVariant': variant['displayName'],
                'kbVariantId': variant['@rid'],
                'matchedCancer': disease_match,
                'reference': pmid,
                'relevance': statement['relevance']['displayName'],
                # TODO: remove, these columns are for debugging but not to output to IPR
                '_source': statement['source']['displayName'] if statement['source'] else None,
                '_sourceId': statement['sourceId'],
            }
            rows.append(row)
    return rows


def create_key_alterations(
    kb_matches: List[Dict],
    expression_variants: List[Dict],
    copy_variants: List[Dict],
    structural_variants: List[Dict],
    small_mutations: List[Dict],
) -> Tuple[List[Dict], Dict]:
    """
    Creates the list of genomic key alterations which summarizes all the variants matched by the KB
    This list of matches is also used to create the variant counts
    """

    def find_variant(kb_match: Dict, variant_type: str, variant_key: str) -> List[Dict]:
        variant_list = None

        if variant_type == 'exp':
            variant_list = expression_variants
        elif variant_type == 'cnv':
            variant_list = copy_variants
        elif variant_type == 'sv':
            variant_list = structural_variants
        else:
            variant_list = small_mutations

        return [v for v in variant_list if v['key'] == variant_key][0]

    alterations = []
    type_mapping = {
        'mut': 'smallMutations',
        'cnv': 'CNVs',
        'sv': "SVs",
        'exp': 'expressionOutliers',
    }
    counts = {v: set() for v in type_mapping.values()}

    for kb_match in kb_matches:
        variant_type = kb_match['variantType']
        variant_key = kb_match['variant']
        variant = find_variant(kb_match, variant_type, variant_key)
        counts[type_mapping[variant_type]].add(variant_key)

        if variant_type in ['exp', 'cnv']:
            gene = variant['gene']
            alterations.append(f'{gene} ({variant["variant"]})')
        else:
            alterations.append(variant['variant'])

    counted_variants = set.union(*counts.values())
    counts['variantsUnknown'] = set()

    # count the un-matched variants
    for variant in expression_variants + structural_variants + copy_variants + small_mutations:
        variant_key = variant['key']

        if variant['variant'] and variant_key not in counted_variants:
            counts['variantsUnknown'].add(variant_key)

    return alterations, {k: len(v) for k, v in counts.items()}


class IprConnection:
    def __init__(self, username: str, password: str, url: str = DEFAULT_URL):
        self.token = None
        self.url = url
        self.username = username
        self.password = password
        self.headers = {
            'Accept': 'application/json',
            'Content-Type': 'application/json',
            'Content-Encoding': 'deflate',
        }
        self.cache = {}
        self.request_count = 0

    def request(self, endpoint: str, method: str = 'GET', **kwargs) -> Dict:
        """Request wrapper to handle adding common headers and logging

        Args:
            endpoint (string): api endpoint, excluding the base uri
            method (str, optional): the http method. Defaults to 'GET'.

        Returns:
            dict: the json response as a python dict
        """
        url = f'{self.url}/{endpoint}'
        self.request_count += 1
        resp = requests.request(
            method, url, headers=self.headers, auth=(self.username, self.password), **kwargs
        )
        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as err:
            # try to get more error details
            message = str(err)
            try:
                message += ' ' + resp.json()['message']
            except Exception:
                pass

            raise requests.exceptions.HTTPError(message)
        return resp.json()

    def post(self, uri: str, data: Dict = {}, **kwargs) -> Dict:
        """Convenience method for making post requests"""
        return self.request(
            uri, method='POST', data=zlib.compress(json.dumps(data).encode('utf-8')), **kwargs
        )

    def upload_report(self, content):
        return self.post('/reports', content)
