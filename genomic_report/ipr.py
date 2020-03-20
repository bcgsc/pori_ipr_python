"""
upload variant and report information to IPR
"""
import requests
import json

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

DEFAULT_URL = 'https://graphkb-api.bcgsc.ca/api'
DEFAULT_LIMIT = 1000


def display_evidence_levels(statement):
    result = []

    for evidence_level in statement.get('evidenceLevel', []) or []:
        result.append(evidence_level['displayName'])

    return ';'.join(sorted(result))


def get_approved_evidence_levels(graphkb_conn):
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


def convert_statements_to_alterations(graphkb_conn, statements, disease_name):
    """
    Given a set of statements matched from graphkb, convert these into their IPR equivalent representations

    Args:
        graphkb_conn (GraphKBConnection): the graphkb connection object
        statements (list.<dict>): list of statement records from graphkb
        disease_name (str): name of the cancer type for the patient being reported on

    Raises:
        ValueError: could not find the disease type in GraphKB

    Returns:
        list.<dict>: IPR graphkb row representations
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
                'kb_entry_key': statement['@rid'],
                'kb_event_key': variant['@rid'],
                'kb_entry_type': ipr_section,
                'alterationType': ipr_section,
                'approvedTherapy': approved_therapy,
                'association': statement['relevance']['displayName'],
                'therapeuticContext': (
                    statement['subject']['displayName'] if statement['subject'] else None
                ),  # may as well include for all statement types. TODO: rename in IPR to be just context
                'kbVariant': variant['displayName'],
                'reference': pmid,
                'evidence': display_evidence_levels(statement),
                'disease': ';'.join(sorted(d['displayName'] for d in diseases)),
                'matched_cancer': disease_match,
                # TODO: remove, these columns are for debugging but not to output to IPR
                '_source': statement['source']['displayName'] if statement['source'] else None,
                '_sourceId': statement['sourceId'],
            }
            rows.append(row)
    return rows


class IprConnection:
    def __init__(self, username, password, url=DEFAULT_URL):
        self.token = None
        self.url = url
        self.username = username
        self.password = password
        self.headers = {'Accept': 'application/json', 'Content-Type': 'application/json'}
        self.cache = {}
        self.request_count = 0

    def request(self, endpoint, method='GET', **kwargs):
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
        resp.raise_for_status()
        return resp.json()

    def post(self, uri, data={}, **kwargs):
        """Convenience method for making post requests"""
        return self.request(uri, method='POST', data=json.dumps(data), **kwargs)

    def upload_report(self, content):
        return self.post('/reports', content)
