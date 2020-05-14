"""
handles annotating variants with annotation information from graphkb
"""
from typing import Dict, List, Set

from graphkb import GraphKBConnection
from graphkb.constants import BASE_RETURN_PROPERTIES, GENERIC_RETURN_PROPERTIES
from graphkb.genes import get_oncokb_oncogenes, get_oncokb_tumour_supressors
from graphkb.match import (
    get_equivalent_features,
    match_copy_variant,
    match_expression_variant,
    match_positional_variant,
)
from graphkb.util import FeatureNotFoundError, convert_to_rid_list
from progressbar import progressbar
from requests.exceptions import HTTPError

from .ipr import BASE_THERAPEUTIC_TERMS, convert_statements_to_alterations, get_terms_set
from .util import convert_to_rid_set, logger


def get_therapeutic_associated_genes(graphkb_conn: GraphKBConnection) -> Set[str]:
    therapeutic_relevance = get_terms_set(graphkb_conn, BASE_THERAPEUTIC_TERMS)
    statements = graphkb_conn.query(
        {
            'target': 'Statement',
            'filters': {'relevance': list(therapeutic_relevance)},
            'returnProperties': [
                'conditions.@rid',
                'conditions.@class',
                'conditions.reference1.@class',
                'conditions.reference1.@rid',
                'conditions.reference2.@class',
                'conditions.reference2.@rid',
            ],
        }
    )
    genes = set()

    for statement in statements:
        for condition in statement['conditions']:
            if condition['@class'] == 'Feature':
                genes.add(condition['@rid'])
            elif condition['@class'].endswith('Variant'):
                if condition['reference1'] and condition['reference1']['@class'] == 'Feature':
                    genes.add(condition['reference1']['@rid'])
                if condition['reference2'] and condition['reference2']['@class'] == 'Feature':
                    genes.add(condition['reference2']['@rid'])
    return genes


def get_gene_information(graphkb_conn: GraphKBConnection, gene_names: List[str]) -> List[Dict]:
    """
    Create the Gene Info object for upload to IPR with the other report information

    Args:
        graphkb_conn ([type]): [description]
        gene_names ([type]): [description]
    """
    logger.info('fetching variant related genes list')
    variants = graphkb_conn.query(
        {'target': 'Variant', 'returnProperties': ['@class', 'reference1', 'reference2']}
    )

    gene_flags = {
        'cancerRelated': set(),
        'knownFusionPartner': set(),
        'knownSmallMutation': set(),
    }

    for variant in variants:
        gene_flags['cancerRelated'].add(variant['reference1'])
        if variant['reference2']:
            gene_flags['cancerRelated'].add(variant['reference2'])
            gene_flags['knownFusionPartner'].add(variant['reference1'])
            gene_flags['knownFusionPartner'].add(variant['reference2'])
        elif variant['@class'] == 'PositionalVariant':
            gene_flags['knownSmallMutation'].add(variant['reference1'])

    logger.info('fetching oncogenes list')
    gene_flags['oncogene'] = convert_to_rid_set(get_oncokb_oncogenes(graphkb_conn))
    logger.info('fetching tumour supressors list')
    gene_flags['tumourSuppressor'] = convert_to_rid_set(get_oncokb_tumour_supressors(graphkb_conn))
    logger.info('fetching therapeutic associated genes lists')
    gene_flags['therapeuticAssociated'] = get_therapeutic_associated_genes(graphkb_conn)

    result = []

    for gene_name in gene_names:
        equivalent = convert_to_rid_set(get_equivalent_features(graphkb_conn, gene_name))

        row = {'name': gene_name}

        for flag in gene_flags:
            row[flag] = bool(equivalent & gene_flags[flag])

        flags = [c for c in row.keys() if c != 'name']

        if any(row[c] for c in flags):
            result.append(row)

            # make smaller JSON to upload since all default to false already
            for flag in flags:
                if not row[flag]:
                    del row[flag]

    return result


def get_statements_from_variants(
    graphkb_conn: GraphKBConnection, variants: List[Dict]
) -> List[Dict]:
    """
    Given a list of variant records from GraphKB, return all the related statements

    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object
        variants (list.<dict>): list of variant records

    Returns:
        list.<dict>: list of Statement records from graphkb
    """
    return_props = (
        BASE_RETURN_PROPERTIES
        + ['sourceId', 'source.name', 'source.displayName']
        + [f'conditions.{p}' for p in GENERIC_RETURN_PROPERTIES]
        + [f'subject.{p}' for p in GENERIC_RETURN_PROPERTIES]
        + [f'evidence.{p}' for p in GENERIC_RETURN_PROPERTIES]
        + [f'relevance.{p}' for p in GENERIC_RETURN_PROPERTIES]
        + [f'evidenceLevel.{p}' for p in GENERIC_RETURN_PROPERTIES]
    )

    statements = graphkb_conn.query(
        {
            'target': 'Statement',
            'filters': {'conditions': convert_to_rid_list(variants), 'operator': 'CONTAINSANY'},
            'returnProperties': return_props,
        }
    )
    return statements


def annotate_category_variants(
    graphkb_conn: GraphKBConnection,
    variants: List[Dict],
    disease_name: str,
    copy_variant: bool = True,
) -> List[Dict]:
    """
    Annotate variant calls with information from GraphKB and return these annotations in the IPR
    alterations format

    Args:
        graphkb_conn: the graphkb api connection object
        variants: list of variants

    Returns:
        list of kbMatches records for IPR
    """
    skipped = 0
    alterations = []
    problem_genes = set()

    logger.info(f"Starting annotation of {len(variants)} category_variants")
    for row in progressbar(variants):
        gene = row['gene']
        variant = row['variant']

        if not variant:
            skipped += 1
            continue

        try:
            if copy_variant:
                matches = match_copy_variant(graphkb_conn, gene, variant)
            else:
                matches = match_expression_variant(graphkb_conn, gene, variant)

            if matches:
                statements = get_statements_from_variants(graphkb_conn, matches)
                for ipr_row in convert_statements_to_alterations(
                    graphkb_conn, statements, disease_name, convert_to_rid_set(matches)
                ):
                    new_row = {'variant': row['key'], 'variantType': row['variantType']}
                    new_row.update(ipr_row)
                    alterations.append(new_row)
        except FeatureNotFoundError as err:
            problem_genes.add(gene)
            logger.debug(f'Unrecognized gene ({gene} {variant}): {err}')
        except ValueError as err:
            logger.error(f'failed to match variants ({gene} {variant}): {err}')

    if skipped:
        logger.info(f'skipped matching {skipped} non variant information rows')
    if problem_genes:
        logger.debug(f'gene finding failures for {sorted(problem_genes)}')
        logger.error(f'gene finding falure for {len(problem_genes)} genes')
    logger.info(
        f'matched {len(variants)} category variants to {len(alterations)} graphkb annotations'
    )
    return alterations


def annotate_positional_variants(
    graphkb_conn: GraphKBConnection, variants: List[Dict], disease_name: str
) -> List[Dict]:
    """
    Annotate variant calls with information from GraphKB and return these annotations in the IPR
    alterations format

    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object
        variants (list.<dict>): list of variants. Defaults to [].

    Returns:
        list of kbMatches records for IPR
    """
    errors = 0
    alterations = []
    problem_genes = set()

    for row in progressbar(variants):
        variant = row['variant']

        if not row.get('gene', '') and (not row.get('gene1', '') or not row.get('gene2', '')):
            # https://www.bcgsc.ca/jira/browse/GERO-56?focusedCommentId=1234791&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-1234791
            # should not match single gene SVs
            continue

        try:
            matches = match_positional_variant(graphkb_conn, variant)

            if matches:
                statements = get_statements_from_variants(graphkb_conn, matches)

                for ipr_row in convert_statements_to_alterations(
                    graphkb_conn, statements, disease_name, convert_to_rid_set(matches)
                ):
                    new_row = {'variant': row['key'], 'variantType': row['variantType']}
                    new_row.update(ipr_row)
                    alterations.append(new_row)
        except FeatureNotFoundError as err:
            logger.debug(f'failed to match positional variants ({variant}): {err}')
            errors += 1
            if 'gene' in row:
                problem_genes.add(row['gene'])
            elif 'gene1' in row and f"({row['gene1']})" in str(err):
                problem_genes.add(row['gene1'])
            elif 'gene2' in row and f"({row['gene2']})" in str(err):
                problem_genes.add(row['gene2'])
            elif 'gene1' in row and 'gene2' in row:
                problem_genes.add(row['gene1'])
                problem_genes.add(row['gene2'])
            else:
                raise err
        except HTTPError as err:
            errors += 1
            logger.error(f'failed to match positional variants ({variant}): {err}')

    if problem_genes:
        logger.error(f'gene finding failures for {sorted(problem_genes)}')
        logger.error(f'{len(problem_genes)} gene finding failures for positional variants')
    if errors:
        logger.error(f'skipped {errors} positional variants due to errors')
    logger.info(
        f'matched {len(variants)} positional variants to {len(alterations)} graphkb annotations'
    )

    return alterations
