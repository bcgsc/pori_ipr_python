"""
handles annotating variants with annotation information from graphkb
"""
from requests.exceptions import HTTPError

from graphkb import GraphKBConnection
from graphkb import genes as gkb_genes
from graphkb import match as gkb_match
from graphkb.constants import STATEMENT_RETURN_PROPERTIES
from graphkb.genes import get_therapeutic_associated_genes
from graphkb.match import INPUT_COPY_CATEGORIES
from graphkb.types import Variant
from graphkb.util import FeatureNotFoundError, convert_to_rid_list
from progressbar import progressbar
from typing import Any, Dict, List, Sequence, Set, cast

from .constants import FAILED_REVIEW_STATUS
from .ipr import convert_statements_to_alterations
from .types import (
    GkbStatement,
    IprCopyVariant,
    IprExprVariant,
    IprGene,
    IprStructuralVariant,
    KbMatch,
)
from .util import Hashabledict, convert_to_rid_set, logger

REPORTED_COPY_VARIANTS = (INPUT_COPY_CATEGORIES.AMP, INPUT_COPY_CATEGORIES.DEEP)


def get_gene_information(
    graphkb_conn: GraphKBConnection, gene_names: Sequence[str]
) -> List[IprGene]:
    """Create the Gene Info object for upload to IPR with the other report information.

    Args:
        graphkb_conn ([type]): [description]
        gene_names ([type]): [description]
    """
    logger.info('fetching variant related genes list')
    body: Dict[str, Any] = {
        'target': 'Variant',
        'returnProperties': ['@class', 'reference1', 'reference2'],
    }
    if len(gene_names) < 100:
        # SDEV-3148 - Filter by gene_ids to improve speed
        gene_ids = set()
        for gene_name in gene_names:
            gene_ids.update(
                convert_to_rid_set(gkb_match.get_equivalent_features(graphkb_conn, gene_name))
            )
        genes = sorted(gene_ids)
        filters = [{'reference1': genes}, {'reference2': genes}]
        variants = []
        for ref_filter in filters:
            body['filters'] = ref_filter
            variants.extend(graphkb_conn.query(body))
    else:
        variants = graphkb_conn.query(body)

    gene_flags: Dict[str, Set[str]] = {
        'cancerRelated': set(),
        'knownFusionPartner': set(),
        'knownSmallMutation': set(),
    }

    for variant in [cast(Variant, v) for v in variants]:
        if 'reference1' not in variant:
            continue
        gene_flags['cancerRelated'].add(variant['reference1'])
        if variant['reference2']:
            gene_flags['cancerRelated'].add(variant['reference2'])
            gene_flags['knownFusionPartner'].add(variant['reference1'])
            gene_flags['knownFusionPartner'].add(variant['reference2'])
        elif variant['@class'] == 'PositionalVariant':
            gene_flags['knownSmallMutation'].add(variant['reference1'])

    logger.info('fetching oncogenes list')
    gene_flags['oncogene'] = convert_to_rid_set(gkb_genes.get_oncokb_oncogenes(graphkb_conn))
    logger.info('fetching tumour supressors list')
    gene_flags['tumourSuppressor'] = convert_to_rid_set(
        gkb_genes.get_oncokb_tumour_supressors(graphkb_conn)
    )
    logger.info('fetching therapeutic associated genes lists')
    gene_flags['therapeuticAssociated'] = convert_to_rid_set(
        get_therapeutic_associated_genes(graphkb_conn)
    )

    result = []
    for gene_name in gene_names:
        equivalent = convert_to_rid_set(gkb_match.get_equivalent_features(graphkb_conn, gene_name))
        row = IprGene({'name': gene_name})
        flagged = False
        for flag in gene_flags:
            # make smaller JSON to upload since all default to false already
            if equivalent & gene_flags[flag]:
                row[flag] = flagged = True
        if flagged:
            result.append(row)

    return result


def get_statements_from_variants(
    graphkb_conn: GraphKBConnection, variants: List[Variant]
) -> List[GkbStatement]:
    """Given a list of variant records from GraphKB, return all the related statements.

    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object
        variants (list.<dict>): list of variant records

    Returns:
        list.<dict>: list of Statement records from graphkb
    """
    statements = graphkb_conn.query(
        {
            'target': 'Statement',
            'filters': {'conditions': convert_to_rid_list(variants), 'operator': 'CONTAINSANY'},
            'returnProperties': STATEMENT_RETURN_PROPERTIES,
        }
    )
    return [
        cast(GkbStatement, s) for s in statements if s.get('reviewStatus') != FAILED_REVIEW_STATUS
    ]


def get_second_pass_variants(
    graphkb_conn: GraphKBConnection, statements: List[GkbStatement]
) -> List[Variant]:
    """Given a list of statements that have been matched, convert these to
    new category variants to be used in a second-pass matching.
    """
    # second-pass matching
    all_inferred_matches: Dict[str, Variant] = {}
    inferred_variants = {
        (s['subject']['@rid'], s['relevance']['name'])
        for s in statements
        if s['subject'] and s['subject']['@class'] in ('Feature', 'Signature')
    }

    for reference1, variant_type in inferred_variants:
        variants = gkb_match.match_category_variant(graphkb_conn, reference1, variant_type)

        for variant in variants:
            all_inferred_matches[variant['@rid']] = variant
    inferred_matches: List[Variant] = list(all_inferred_matches.values())
    return inferred_matches


def get_ipr_statements_from_variants(
    graphkb_conn: GraphKBConnection, matches: List[Variant], disease_name: str
) -> List[KbMatch]:
    """IPR upload formatted GraphKB statements from the list of variants.

    Matches to GraphKB statements from the list of input variants. From these results matches
    again with the inferred variants. Then returns the results formatted for upload to IPR
    """
    if not matches:
        return []
    rows = []

    statements = get_statements_from_variants(graphkb_conn, matches)
    existing_statements = {s['@rid'] for s in statements}

    for ipr_row in convert_statements_to_alterations(
        graphkb_conn, statements, disease_name, convert_to_rid_set(matches)
    ):
        rows.append(ipr_row)

    # second-pass matching
    inferred_matches = get_second_pass_variants(graphkb_conn, statements)

    inferred_statements = [
        s
        for s in get_statements_from_variants(graphkb_conn, inferred_matches)
        if s['@rid'] not in existing_statements  # do not duplicate if non-inferred match
    ]

    for ipr_row in convert_statements_to_alterations(
        graphkb_conn, inferred_statements, disease_name, convert_to_rid_set(inferred_matches)
    ):
        new_row = KbMatch({'kbData': {'inferred': True}})
        new_row.update(ipr_row)
        rows.append(new_row)

    return rows


def annotate_expression_variants(
    graphkb_conn: GraphKBConnection,
    variants: List[IprExprVariant],
    disease_name: str,
    show_progress: bool = False,
) -> List[KbMatch]:
    """Annotate expression variants with GraphKB in the IPR alterations format.

    Args:
        graphkb_conn: the graphkb api connection object
        variants: list of variants

    Returns:
        list of kbMatches records for IPR
    """
    skipped = 0
    alterations = []
    problem_genes = set()

    logger.info(f"Starting annotation of {len(variants)} expression category_variants")
    iterfunc = progressbar if show_progress else iter
    for row in iterfunc(variants):
        gene = row['gene']
        variant = row['variant']

        if not variant:
            skipped += 1
            logger.debug(f"Skipping malformed Expression {gene}: {row}")
            continue

        try:
            matches = gkb_match.match_expression_variant(graphkb_conn, gene, variant)
            for ipr_row in get_ipr_statements_from_variants(graphkb_conn, matches, disease_name):
                ipr_row['variant'] = row['key']
                ipr_row['variantType'] = row.get('variantType', 'exp')
                alterations.append(ipr_row)
        except FeatureNotFoundError as err:
            problem_genes.add(gene)
            logger.debug(f'Unrecognized gene ({gene} {variant}): {err}')
        except ValueError as err:
            logger.error(f'failed to match variants ({gene} {variant}): {err}')

    if skipped:
        logger.info(f'skipped matching {skipped} expression information rows')
    if problem_genes:
        logger.error(f'gene finding failures for expression {sorted(problem_genes)}')
        logger.error(f'gene finding falure for {len(problem_genes)} expression genes')
    logger.info(
        f'matched {len(variants)} expression variants to {len(alterations)} graphkb annotations'
    )
    return alterations


def annotate_copy_variants(
    graphkb_conn: GraphKBConnection,
    variants: List[IprCopyVariant],
    disease_name: str,
    show_progress: bool = False,
) -> List[KbMatch]:
    """Annotate allowed copy variants with GraphKB in the IPR alterations format.

    Args:
        graphkb_conn: the graphkb api connection object
        variants: list of variants

    Returns:
        list of kbMatches records for IPR
    """
    skipped = 0
    alterations = []
    problem_genes = set()

    logger.info(f"Starting annotation of {len(variants)} copy category_variants")
    iterfunc = progressbar if show_progress else iter
    for row in iterfunc(variants):
        gene = row['gene']
        variant = row['variant']

        if variant not in REPORTED_COPY_VARIANTS:
            # https://www.bcgsc.ca/jira/browse/GERO-77
            skipped += 1
            logger.debug(f"Dropping {gene} copy change '{variant}' - not in REPORTED_COPY_VARIANTS")
            continue

        try:
            matches = gkb_match.match_copy_variant(graphkb_conn, gene, variant)
            for ipr_row in get_ipr_statements_from_variants(graphkb_conn, matches, disease_name):
                ipr_row['variant'] = row['key']
                ipr_row['variantType'] = row.get('variantType', 'cnv')
                alterations.append(ipr_row)
        except FeatureNotFoundError as err:
            problem_genes.add(gene)
            logger.debug(f'Unrecognized gene ({gene} {variant}): {err}')
        except ValueError as err:
            logger.error(f'failed to match variants ({gene} {variant}): {err}')

    if skipped:
        logger.info(
            f'skipped matching {skipped} copy number variants not in {REPORTED_COPY_VARIANTS}'
        )
    if problem_genes:
        logger.error(f'gene finding failures for copy variants {sorted(problem_genes)}')
        logger.error(f'gene finding failure for {len(problem_genes)} copy variant genes')
    logger.info(
        f'matched {len(variants)} copy category variants to {len(alterations)} graphkb annotations'
    )
    return alterations


def annotate_positional_variants(
    graphkb_conn: GraphKBConnection,
    variants: Sequence[IprStructuralVariant],
    disease_name: str,
    show_progress: bool = False,
) -> List[KbMatch]:
    """Annotate SNP, INDEL or fusion variant calls with GraphKB and return in IPR match format.

    Args:
        graphkb_conn (GraphKBConnection): the graphkb api connection object
        variants (list.<dict>): list of variants. Defaults to [].
        disease_name (str): GraphKB disease name for statement matching.  'cancer' is most general
        show_progress (bool): Progressbar displayed for long runs.

    Returns:
        list of kbMatches records for IPR
    """
    VARIANT_KEYS = ('variant', 'hgvsProtein', 'hgvsCds', 'hgvsGenomic')
    errors = 0
    alterations = []
    problem_genes = set()

    iterfunc = progressbar if show_progress else iter
    for row in iterfunc(variants):
        if not row.get('gene') and (not row.get('gene1') or not row.get('gene2')):
            # https://www.bcgsc.ca/jira/browse/GERO-56?focusedCommentId=1234791&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-1234791
            # should not match single gene SVs
            continue

        for var_key in VARIANT_KEYS:
            variant = row.get(var_key)
            if not variant:
                continue
            try:
                matches = gkb_match.match_positional_variant(graphkb_conn, variant)

                # GERO-299 - check for conflicting nonsense and missense categories
                missense = [
                    m for m in matches if 'missense' in m.get('type', m).get('displayName', '')
                ]
                nonsense = [
                    m for m in matches if 'nonsense' in m.get('type', m).get('displayName', '')
                ]
                missense_cat = [m for m in missense if m.get('@class', '') == 'CategoryVariant']
                nonsense_cat = [m for m in nonsense if m.get('@class', '') == 'CategoryVariant']
                if missense_cat and nonsense_cat:
                    conflict_names = sorted(
                        set([m.get('displayName', '') for m in nonsense + missense])
                    )
                    logger.error(
                        f'Conflicting nonsense and misense categories for {variant}: {conflict_names}'
                    )
                    # Check if cross referenced positional variants resolve the category conflict.
                    if nonsense_cat == nonsense and missense_cat != missense:
                        cross_match = sorted(
                            set(
                                [
                                    m.get('displayName', '')
                                    for m in missense
                                    if m not in missense_cat
                                ]
                            )
                        )
                        logger.error(f"GERO-299 - dropping nonsense category due to: {cross_match}")
                        matches = [m for m in matches if m not in nonsense_cat]
                    elif nonsense_cat != nonsense and missense_cat == missense:
                        cross_match = sorted(
                            set(
                                [
                                    m.get('displayName', '')
                                    for m in nonsense
                                    if m not in nonsense_cat
                                ]
                            )
                        )
                        logger.error(f"GERO-299 - dropping missense category due to: {cross_match}")
                        matches = [m for m in matches if m not in missense_cat]
                    else:
                        conflict_names = sorted(
                            set([m.get('displayName', '') for m in nonsense_cat + missense_cat])
                        )
                        logger.error(
                            f"GERO-299 Dropping all conflicting missense/nonsense categories: {conflict_names}"
                        )
                        matches = [
                            m for m in matches if m not in missense_cat and m not in nonsense_cat
                        ]

                for ipr_row in get_ipr_statements_from_variants(
                    graphkb_conn, matches, disease_name
                ):
                    ipr_row['variant'] = row['key']
                    ipr_row['variantType'] = row.get(
                        'variantType', 'mut' if row.get('gene') else 'sv'
                    )
                    alterations.append(Hashabledict(ipr_row))

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

    # drop duplicates
    alterations: List[KbMatch] = list(set(alterations))
    logger.info(
        f'matched {len(variants)} positional variants to {len(alterations)} graphkb annotations'
    )

    return alterations


def annotate_msi(
    graphkb_conn: GraphKBConnection,
    msi_category: str,
    disease_name: str,
) -> List[KbMatch]:
    """Annotate microsatellite instablity from GraphKB in the IPR alterations format.

    Match to GraphKb Category variants with similar names
    Args:
        graphkb_conn: the graphkb api connection object
        msi_category: such as 'microsatellite instability'

    Returns:
        list of kbMatches records for IPR
    """
    gkb_matches = []
    msi_categories = graphkb_conn.query(
        {
            'target': {
                'target': 'CategoryVariant',
                'filters': {
                    'reference1': {'target': 'Signature', 'filters': {'name': msi_category}}
                },
            },
            'queryType': 'similarTo',
            'returnProperties': ['@rid', 'displayName'],
        },
    )
    if msi_categories:
        for ipr_row in get_ipr_statements_from_variants(graphkb_conn, msi_categories, disease_name):
            ipr_row['variant'] = msi_category
            ipr_row['variantType'] = 'msi'
            gkb_matches.append(ipr_row)
    return gkb_matches
