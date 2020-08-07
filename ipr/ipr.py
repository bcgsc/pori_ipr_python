"""
upload variant and report information to IPR
"""
from typing import Dict, Iterable, List, Set, Tuple

from graphkb import GraphKBConnection
from graphkb.types import Ontology, Statement
from graphkb.vocab import get_term_tree

from .types import ImageDefinition, IprGene, IprStructuralVariant, IprVariant, KbMatch
from .util import (
    convert_to_rid_set,
    find_variant,
    get_terms_set,
    create_variant_name,
)
from .constants import (
    BASE_BIOLOGICAL_TERMS,
    BASE_DIAGNOSTIC_TERM,
    BASE_PROGNOSTIC_TERM,
    BASE_THERAPEUTIC_TERMS,
    VARIANT_CLASSES,
    REPORT_KB_SECTIONS,
    APPROVED_EVIDENCE_LEVELS,
)


def display_evidence_levels(statement: Statement) -> str:
    result = []

    for evidence_level in statement.get('evidenceLevel', []) or []:
        result.append(evidence_level['displayName'])

    return ';'.join(sorted(result))


def get_approved_evidence_levels(graphkb_conn: GraphKBConnection) -> List[Ontology]:
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


def filter_structural_variants(
    structural_variants: List[IprStructuralVariant],
    kb_matches: List[KbMatch],
    gene_annotations: List[IprGene],
) -> List[IprStructuralVariant]:
    """
    Filter structural variants to remove non-high quality events unless they are matched/annotated or
    they involve a gene that is a known fusion partner
    """
    matched_svs = {match['variant'] for match in kb_matches if match['variantType'] == 'sv'}
    fusion_genes = {
        gene['name'] for gene in gene_annotations if gene.get('knownFusionPartner', False)
    }

    result = []

    for structural_variant in structural_variants:
        if any(
            [
                structural_variant['highQuality'],
                structural_variant['key'] in matched_svs,
                structural_variant['gene1'] in fusion_genes,
                structural_variant['gene2'] in fusion_genes,
            ]
        ):
            result.append(structural_variant)
    return result


def convert_statements_to_alterations(
    graphkb_conn: GraphKBConnection,
    statements: List[Statement],
    disease_name: str,
    variant_matches: Iterable[str],
) -> List[KbMatch]:
    """
    Given a set of statements matched from graphkb, convert these into their IPR equivalent representations

    Args:
        graphkb_conn: the graphkb connection object
        statements: list of statement records from graphkb
        disease_name: name of the cancer type for the patient being reported on
        variant_matches: the list of RIDs the variant matched for these statements

    Raises:
        ValueError: could not find the disease type in GraphKB

    Returns:
        IPR graphkb row representations

    Notes:
        - only report disease matched prognostic markers https://www.bcgsc.ca/jira/browse/GERO-72
    """
    disease_matches = {
        r['@rid'] for r in get_term_tree(graphkb_conn, disease_name, ontology_class='Disease')
    }

    if not disease_matches:
        raise ValueError(f'failed to match disease ({disease_name}) to graphkb')

    rows = []

    approved = convert_to_rid_set(get_approved_evidence_levels(graphkb_conn))

    therapeutic_terms = get_terms_set(graphkb_conn, BASE_THERAPEUTIC_TERMS)
    diagnostic_terms = get_terms_set(graphkb_conn, [BASE_DIAGNOSTIC_TERM])
    prognostic_terms = get_terms_set(graphkb_conn, [BASE_PROGNOSTIC_TERM])
    biological_terms = get_terms_set(graphkb_conn, BASE_BIOLOGICAL_TERMS)

    for statement in statements:
        variants = [c for c in statement['conditions'] if c['@class'] in VARIANT_CLASSES]
        diseases = [c for c in statement['conditions'] if c['@class'] == 'Disease']
        pmid = ';'.join([e['displayName'] for e in statement['evidence']])

        ipr_section = REPORT_KB_SECTIONS.unknown  # table this goes to in IPR
        relevance_id = statement['relevance']['@rid']

        approved_therapy = False

        disease_match = len(diseases) == 1 and diseases[0]['@rid'] in disease_matches

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
            if not disease_match:
                continue  # GERO-72
        elif relevance_id in biological_terms:
            ipr_section = REPORT_KB_SECTIONS.biological

        for variant in variants:
            if variant['@rid'] not in variant_matches:
                continue
            row = KbMatch(
                {
                    'approvedTherapy': approved_therapy,
                    'category': ipr_section,
                    'context': (
                        statement['subject']['displayName'] if statement['subject'] else None
                    ),
                    'kbContextId': (statement['subject']['@rid'] if statement['subject'] else None),
                    'disease': ';'.join(sorted(d['displayName'] for d in diseases)),
                    'evidenceLevel': display_evidence_levels(statement),
                    'kbStatementId': statement['@rid'],
                    'kbVariant': variant['displayName'],
                    'kbVariantId': variant['@rid'],
                    'matchedCancer': disease_match,
                    'reference': pmid,
                    'relevance': statement['relevance']['displayName'],
                    'kbRelevanceId': statement['relevance']['@rid'],
                }
            )
            rows.append(row)
    return rows


def select_expression_plots(
    kb_matches: List[KbMatch], all_variants: List[IprVariant]
) -> List[Dict[str, ImageDefinition]]:
    """
    Given the list of expression variants, determine which expression
    historgram plots should be included in the IPR upload. This filters them
    based on the graphkb annotations to avoid loading more images than are required

    Args:
        kb_matches: the IPR graphkb annoations for all variants
        expression_variants: the list of expression variants loaded

    Returns:
        list of expression images to be loaded by IPR
    """

    selected_variants = {
        (match['variantType'], match['variant'])
        for match in kb_matches
        if match['category'] == 'therapeutic'
    }
    images_by_gene: Dict[str, ImageDefinition] = {}
    selected_genes = set()
    for variant in all_variants:
        if (variant['variantType'], variant['key']) in selected_variants:
            for key in ['gene', 'gene1', 'gene2']:
                if key in variant and variant[key]:
                    selected_genes.add(variant[key])
        if variant.get('histogramImage', ''):
            gene = variant['gene']
            images_by_gene[gene] = ImageDefinition(
                {'key': f'expDensity.{gene}', 'path': variant['histogramImage']}
            )
    return [images_by_gene[gene] for gene in selected_genes if gene in images_by_gene]


def create_key_alterations(
    kb_matches: List[KbMatch], all_variants: List[IprVariant],
) -> Tuple[List[Dict], Dict]:
    """
    Creates the list of genomic key alterations which summarizes all the variants matched by the KB
    This list of matches is also used to create the variant counts
    """

    alterations = []
    type_mapping = {
        'mut': 'smallMutations',
        'cnv': 'CNVs',
        'sv': "SVs",
        'exp': 'expressionOutliers',
    }
    counts: Dict[str, Set] = {v: set() for v in type_mapping.values()}

    for kb_match in kb_matches:
        variant_type = kb_match['variantType']
        variant_key = kb_match['variant']
        variant = find_variant(all_variants, variant_type, variant_key)
        counts[type_mapping[variant_type]].add(variant_key)

        if kb_match['category'] == 'unknown':
            continue

        alterations.append(create_variant_name(variant))

    counted_variants = set.union(*counts.values())
    counts['variantsUnknown'] = set()

    # count the un-matched variants
    for variant in all_variants:
        variant_key = variant['key']

        if variant['variant'] and variant_key not in counted_variants:
            counts['variantsUnknown'].add(variant_key)

    return (
        [{'geneVariant': alt} for alt in set(alterations)],
        {k: len(v) for k, v in counts.items()},
    )