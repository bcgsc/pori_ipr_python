"""
upload variant and report information to IPR
"""

from graphkb.util import IterableNamespace, convert_to_rid_list
from graphkb.vocab import get_term_tree


BASE_THERAPEUTIC_TERM = 'therapeutic efficacy'
BASE_DIAGNOSITC_TERM = 'diagnostic indicator'
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


def display_evidence_levels(statement):
    result = []

    for evidence_level in statement.get('evidenceLevel', []) or []:
        result.append(evidence_level['displayName'])

    return ';'.join(sorted(result))


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

    therapeutic_terms = set(
        convert_to_rid_list(
            get_term_tree(graphkb_conn, BASE_THERAPEUTIC_TERM, include_superclasses=False)
        )
    )
    diagnostic_terms = set(
        convert_to_rid_list(
            get_term_tree(graphkb_conn, BASE_DIAGNOSITC_TERM, include_superclasses=False)
        )
    )
    prognostic_terms = set(
        convert_to_rid_list(
            get_term_tree(graphkb_conn, BASE_PROGNOSTIC_TERM, include_superclasses=False)
        )
    )
    biological_terms = set(
        convert_to_rid_list(
            get_term_tree(graphkb_conn, BASE_BIOLOGICAL_TERM, include_superclasses=False)
        )
    )

    for statement in statements:
        variants = [c for c in statement['conditions'] if c['@class'] in VARIANT_CLASSES]
        diseases = [c for c in statement['conditions'] if c['@class'] == 'Disease']
        pmid = ';'.join([e['displayName'] for e in statement['evidence']])

        ipr_section = REPORT_KB_SECTIONS.unknown  # table this goes to in IPR
        relevance_id = statement['relevance']['@rid']

        if relevance_id in therapeutic_terms:
            ipr_section = REPORT_KB_SECTIONS.therapeutic
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
