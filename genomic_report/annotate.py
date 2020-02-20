"""
handles annotating variants with annotation information from graphkb
"""
from Bio.Data.IUPACData import protein_letters_3to1
from graphkb.match import (
    match_copy_variant,
    # match_expression_variant,
    match_positional_variant,
)
from graphkb.util import convert_to_rid_list
from graphkb.constants import BASE_RETURN_PROPERTIES, GENERIC_RETURN_PROPERTIES

from .ipr import convert_statements_to_alterations


protein_letters_3to1.setdefault('Ter', '*')


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


def annotate_variants(
    graphkb_conn,
    copy_variants=[],
    expression_variants=[],
    structural_variants=[],
    small_mutations=[],
):
    """
    Annotate variant calls with information from GraphKB and return these annotations in the IPR
    alterations format

    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object
        copy_variants (list, optional): list of copy number variants. Defaults to [].
        expression_variants (list, optional): list of expression variants. Defaults to [].
        structural_variants (list, optional): list of structural variants. Defaults to [].
        small_mutations (list, optional): list of small mutations. Defaults to [].

    Returns:
        [type]: [description]
    """
    skipped = 0
    errors = 0
    alterations = []

    if copy_variants:
        for row in copy_variants:

            if not row['variant']:
                skipped += 1
                continue

            try:
                matches = match_copy_variant(graphkb_conn, row['gene'], row['variant'])

                if matches:
                    statements = get_statements_from_variants(graphkb_conn, matches)
                    for ipr_row in convert_statements_to_alterations(
                        graphkb_conn, statements, 'colorectal cancer'
                    ):
                        new_row = {
                            'gene': row['gene'],
                            'variant': row['variant'],
                            '_variant_key': row['key'],
                        }
                        new_row.update(ipr_row)
                        alterations.append(new_row)
            except ValueError:
                errors += 1

        print(f'skipped matching {skipped} neutral copy variants')
        print(f'skipped {errors} copy number variants due to errors')

    if small_mutations:

        skipped = 0
        errors = 0

        for row in small_mutations:
            for longAA, shortAA in protein_letters_3to1.items():
                row['proteinChange'] = row['proteinChange'].replace(longAA, shortAA)

        for row in small_mutations:
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
            except ValueError:
                errors += 1
            except Exception as err:
                errors += 1
                print(hgvsp, 'ERROR')
                print(err)

        print(f'skipped {errors} small mutations due to errors')

    return alterations
