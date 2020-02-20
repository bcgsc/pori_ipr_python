"""
handles annotating variants with annotation information from graphkb
"""
from graphkb.match import (
    match_copy_variant,
    match_positional_variant,
)
from graphkb.util import convert_to_rid_list
from graphkb.constants import BASE_RETURN_PROPERTIES, GENERIC_RETURN_PROPERTIES

from .ipr import convert_statements_to_alterations
from .util import logger


def get_statements_from_variants(graphkb_conn, variants):
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


def annotate_copy_variants(graphkb_conn, variants):
    """
    Annotate variant calls with information from GraphKB and return these annotations in the IPR
    alterations format

    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object
        variants (list.<dict>): list of copy number variants

    Returns:
        [type]: [description]
    """
    skipped = 0
    errors = 0
    alterations = []
    for row in variants:
        gene = row['gene']
        variant = row['variant']

        if not variant:
            skipped += 1
            continue

        try:
            matches = match_copy_variant(graphkb_conn, gene, variant)

            if matches:
                statements = get_statements_from_variants(graphkb_conn, matches)
                for ipr_row in convert_statements_to_alterations(
                    graphkb_conn, statements, 'colorectal cancer'
                ):
                    new_row = {
                        'gene': gene,
                        'variant': variant,
                        '_variant_key': row['key'],
                    }
                    new_row.update(ipr_row)
                    alterations.append(new_row)
        except ValueError as err:
            logger.warning(f'failed to match copy variants ({gene} {variant}): {err}')
            errors += 1

    logger.info(f'skipped matching {skipped} neutral copy variants')
    logger.info(f'skipped {errors} copy number variants due to errors')
    return alterations


def annotate_small_mutations(graphkb_conn, variants):
    """
    Annotate variant calls with information from GraphKB and return these annotations in the IPR
    alterations format

    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object
        small_mutations (list.<dict>): list of small mutations. Defaults to [].

    Returns:
        [type]: [description]
    """
    errors = 0
    alterations = []

    for row in variants:
        hgvsp = '{}:{}'.format(row['gene'], row['proteinChange'])
        try:
            matches = match_positional_variant(graphkb_conn, hgvsp)

            if matches:
                statements = get_statements_from_variants(graphkb_conn, matches)

                for ipr_row in convert_statements_to_alterations(
                    graphkb_conn, statements, 'colorectal cancer'
                ):
                    new_row = {
                        'gene': row['gene'],
                        'variant': row['proteinChange'],
                        '_variant_key': row['key'],
                    }
                    new_row.update(ipr_row)
                    alterations.append(new_row)
        except ValueError as err:
            errors += 1
            logger.warning(f'failed to match small mutation ({hgvsp}): {err}')
        except Exception as err:
            errors += 1
            logger.error(f'failed to match small mutation ({hgvsp}): {err}')

    logger.info(f'skipped {errors} small mutations due to errors')

    return alterations


def annotate_structural_variants():
    raise NotImplementedError('TODO')


def annotate_expression_variants():
    raise NotImplementedError('TODO')
