"""
handles annotating variants with annotation information from graphkb
"""
from graphkb.match import (
    match_copy_variant,
    match_positional_variant,
    match_expression_variant,
    get_equivalent_features,
)
from graphkb.genes import get_oncokb_oncogenes, get_oncokb_tumour_supressors
from graphkb.util import convert_to_rid_list
from graphkb.constants import BASE_RETURN_PROPERTIES, GENERIC_RETURN_PROPERTIES

from .ipr import convert_statements_to_alterations
from .util import logger, convert_to_rid_set


def get_cancer_related_genes(graphkb_conn):
    """
    Get the list of genes on any variant in GKB

    Args:
        graphkb_conn (GraphKBConnection): graphkb connection object

    Returns:
        set.<str>: set of records IDs
    """
    # cancer related means any gene that has a variant in GKB
    variants = graphkb_conn.query(
        {'target': 'Variant', 'returnProperties': ['reference1', 'reference2']}
    )

    genes = set()
    for variant in variants:
        genes.add(variant['reference1'])
        if variant['reference2']:
            genes.add(variant['reference2'])
    return genes


def get_gene_information(graphkb_conn, gene_names):
    """
    Create the Gene Info object for upload to IPR with the other report information

    Args:
        graphkb_conn ([type]): [description]
        gene_names ([type]): [description]
    """
    logger.verbose('fetching oncogenes list')
    oncogenes = convert_to_rid_set(get_oncokb_oncogenes(graphkb_conn))
    logger.verbose('fetching tumour supressors list')
    tumour_suppressors = convert_to_rid_set(get_oncokb_tumour_supressors(graphkb_conn))
    logger.verbose('fetching cancer related genes list')
    cancer_related = get_cancer_related_genes(graphkb_conn)

    result = []

    for gene_name in gene_names:
        equivalent = convert_to_rid_set(get_equivalent_features(graphkb_conn, gene_name))

        row = {
            'name': gene_name,
            'oncogene': bool(equivalent & oncogenes),
            'tumourSuppressor': bool(equivalent & tumour_suppressors),
            'drugTargetable': False,
            'cancerGene': bool(equivalent & cancer_related),
        }
        flags = [c for c in row.keys() if c != 'name']

        if any(row[c] for c in flags):
            result.append(row)

            # make smaller JSON to upload since all default to false already
            for flag in flags:
                if not row[flag]:
                    del row[flag]

    return result


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


def annotate_category_variants(graphkb_conn, variants, disease_name, copy_variant=True):
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
            if copy_variant:
                matches = match_copy_variant(graphkb_conn, gene, variant)
            else:
                matches = match_expression_variant(graphkb_conn, gene, variant)

            if matches:
                statements = get_statements_from_variants(graphkb_conn, matches)
                for ipr_row in convert_statements_to_alterations(
                    graphkb_conn, statements, disease_name
                ):
                    new_row = {
                        'gene': gene,
                        'variant': variant,
                        '_variant_key': row['key'],
                    }
                    new_row.update(ipr_row)
                    alterations.append(new_row)
        except ValueError as err:
            logger.warning(f'failed to match variants ({gene} {variant}): {err}')
            errors += 1

    logger.info(f'skipped matching {skipped} non variant information rows')
    logger.info(f'skipped {errors} variants due to errors')
    logger.info(f'matched {len(variants)} variants to {len(alterations)} graphkb annotations')
    return alterations


def annotate_positional_variants(graphkb_conn, variants, disease_name):
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
        variant = row['variant']
        try:
            matches = match_positional_variant(graphkb_conn, variant)

            if matches:
                statements = get_statements_from_variants(graphkb_conn, matches)

                for ipr_row in convert_statements_to_alterations(
                    graphkb_conn, statements, disease_name
                ):
                    new_row = {
                        'variant': variant,
                        '_variant_key': row['key'],
                    }
                    new_row.update(ipr_row)
                    alterations.append(new_row)
        except ValueError as err:
            errors += 1
            logger.warning(f'failed to match positional variants ({variant}): {err}')
        except Exception as err:
            errors += 1
            logger.error(f'failed to match positional variants ({variant}): {err}')

    logger.info(f'skipped {errors} positional variants due to errors')
    logger.info(f'matched {len(variants)} variants to {len(alterations)} graphkb annotations')

    return alterations
