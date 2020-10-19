"""
Read/Validate the variant input files
"""
import os
import re
from Bio.Data.IUPACData import protein_letters_3to1
from csv import DictReader
from graphkb.match import INPUT_COPY_CATEGORIES, INPUT_EXPRESSION_CATEGORIES
from typing import Callable, Dict, Iterable, List, Set, Tuple, cast

from .types import IprGeneVariant, IprStructuralVariant, IprVariant
from .util import hash_key, logger

protein_letters_3to1.setdefault('Ter', '*')

NULLABLE_FLOAT_REGEX = r'^-?((inf)|(\d+(\.\d+)?)|)$'
# 'cnvState' is for display
COPY_REQ = ['gene', 'kbCategory']
COPY_KEY = ['gene']
COPY_OPTIONAL = [
    'cnvState',
    'copyChange',
    'lohState',  # Loss of Heterzygosity state - informative detail to analyst
    'chromosomeBand',
    'start',
    'end',
    'size',
    'log2Cna',
    'cna',
]

SMALL_MUT_REQ = ['gene', 'proteinChange']
# alternate details in the key, can distinguish / subtype events.
SMALL_MUT_KEY = SMALL_MUT_REQ + [
    'altSeq',
    'chromosome',
    'endPosition',
    'refSeq',
    'startPosition',
    'transcript',
]
SMALL_MUT_OPTIONAL = [
    'altSeq',
    'chromosome',
    'endPosition',
    'hgvsCds',
    'hgvsGenomic',
    'hgvsProtein',
    'ncbiBuild',
    'normalAltCount',
    'normalDepth',
    'normalRefCount',
    'refSeq',
    'rnaAltCount',
    'rnaDepth',
    'rnaRefCount',
    'startPosition',
    'transcript',
    'tumourAltCount',
    'tumourDepth',
    'tumourRefCount',
    'zygosity',
]

EXP_REQ = ['gene', 'kbCategory']
EXP_KEY = ['gene']
EXP_OPTIONAL = [
    'biopsySiteFoldChange',
    'biopsySitePercentile',
    'biopsySiteQC',
    'biopsySiteZScore',
    'biopsySitekIQR',
    'diseaseFoldChange',
    'diseasekIQR',
    'diseasePercentile',
    'diseaseQC',
    'diseaseZScore',
    'expressionState',
    'histogramImage',
    'primarySiteFoldChange',
    'primarySitekIQR',
    'primarySitePercentile',
    'primarySiteQC',
    'primarySiteZScore',
    'rnaReads',
    'rpkm',
    'tpm',
]

SV_REQ = [
    'eventType',
    'breakpoint',
    'gene1',  # prev: nterm_hugo
    'gene2',  # prev: cterm_hugo
    'exon1',  # n-terminal
    'exon2',  # c-terminal
]
SV_KEY = SV_REQ[:]
SV_OPTIONAL = [
    'ctermTranscript',
    'ntermTranscript',
    'ctermGene',  # combined hugo ensembl form
    'ntermGene',  # combined hugo ensembl form
    'detectedIn',
    'conventionalName',
    'svg',
    'svgTitle',
    'name',
    'frame',
    'omicSupport',
    'highQuality',
]


def validate_variant_rows(
    rows: Iterable[Dict], required: List[str], optional: List[str], row_to_key: Callable
) -> List[IprVariant]:
    """
    - check that the required columns are present
    - check that a unique key can be formed for each row
    - drop any non-defined columns

    Args:
        rows: the input files rows
        required list of required column names
        optional: list of optional column names
        row_to_key: function to generate a key for a given row

    Raises:
        ValueError: row keys are not unique
        ValueError: A required column is missing

    Returns:
        the rows from the tab file as dictionaries
    """
    header = required + optional + ['key']

    result = []
    keys = set()

    header_validated = False

    for row in rows:
        if not header_validated:
            for req_col in required:
                if req_col not in row:
                    raise ValueError(f'header missing required column ({req_col})')
            header_validated = True
        row_key = hash_key(row_to_key(row))
        if row_key in keys:
            raise ValueError(f'duplicate row key ({row_key}) from ({row_to_key(row)})')
        row['key'] = row_key
        keys.add(row_key)

        result.append(cast(IprVariant, {col: row.get(col, '') for col in header}))

    return result


def read_tabbed_file(filename: str) -> List[Dict]:
    """
    Load a tab delimited file as dictionary rows
    """
    with open(filename, 'r') as fh:
        reader = DictReader(fh, delimiter='\t')
        return [row for row in reader]


def validate_row_patterns(
    rows: List[IprVariant], patterns: Dict, row_key_columns: List[str]
) -> None:
    """
    Validate rows against a regex for some set of columns

    Args:
        rows (List.<dict>): input rows read from a delimited file
        patterns (dict.<str,str>): mapping of column names to regex patterns the column are expected to match

    Raises:
        ValueError: A row does not match the expected pattern for a given column
    """
    for row in rows:
        for col, pattern in patterns.items():
            if not re.match(pattern, '' if row.get(col, None) is None else str(row[col])):
                row_repr_dict = dict((key, row.get(key, '')) for key in row_key_columns)
                raise ValueError(
                    f'{row_repr_dict} column {col}: "{row[col]}" re pattern failure: "{pattern}"'
                )


def preprocess_copy_variants(rows: Iterable[Dict]) -> List[IprVariant]:
    """
    Validate the input rows contain the minimum required fields and
    generate any default values where possible
    """
    # default map for display - concise names
    display_name_mapping = {
        INPUT_COPY_CATEGORIES.DEEP: "deep deletion",
        INPUT_COPY_CATEGORIES.AMP: "amplification",
        INPUT_COPY_CATEGORIES.GAIN: "copy gain",
        INPUT_COPY_CATEGORIES.LOSS: "copy loss",
    }

    def row_key(row: Dict) -> Tuple[str, ...]:
        return tuple(['cnv'] + [row[key] for key in COPY_KEY])

    result = validate_variant_rows(rows, COPY_REQ, COPY_OPTIONAL, row_key)

    patterns = {
        'chromosomeBand': r'^(\S+:\S+?)?$',
    }
    validate_row_patterns(result, patterns, COPY_KEY)

    for row in result:
        if row['kbCategory']:
            if row['kbCategory'] not in INPUT_COPY_CATEGORIES.values():
                raise ValueError(f'invalid copy variant kbCategory value ({row["kbCategory"]})')
            if not row['cnvState']:  # apply default short display name
                row['cnvState'] = display_name_mapping[row['kbCategory']]
        row['variant'] = row['kbCategory']
        row['variantType'] = 'cnv'

    return result


def preprocess_small_mutations(rows: Iterable[Dict]) -> List[IprGeneVariant]:
    """
    Validate the input rows contain the minimum required fields and
    generate any default values where possible
    """

    def row_key(row: Dict) -> Tuple[str, ...]:
        return tuple(['small mutation'] + [row.get(key, '') for key in SMALL_MUT_KEY])

    result = validate_variant_rows(rows, SMALL_MUT_REQ, SMALL_MUT_OPTIONAL, row_key)
    if not result:
        return result

    # 'location' and 'refAlt' are not currently used for matching; still optional and allowed blank
    patterns = {
        'hgvsProtein': r'^(\S+:p\.\S+)?$',
        'hgvsCds': r'^(\S+:[crn]\.\S+)?$',
        'hgvsGenomic': r'^(\S+:g\.\S+)?$',
    }
    validate_row_patterns(result, patterns, SMALL_MUT_KEY)

    # change 3 letter AA to 1 letter AA notation
    for row in result:
        for longAA, shortAA in protein_letters_3to1.items():
            row['proteinChange'] = row['proteinChange'].replace(longAA, shortAA)
        hgvsp = '{}:{}'.format(row['gene'], row['proteinChange'])
        row['variant'] = hgvsp
        row['variantType'] = 'mut'

        if row.get('startPosition') and not row.get('endPosition'):
            row['endPosition'] = row['startPosition']

        # check integer columns
        for col in [
            'endPosition',
            'normalAltCount',
            'normalDepth',
            'normalRefCount',
            'rnaAltCount',
            'rnaDepth',
            'rnaRefCount',
            'startPosition',
            'tumourAltCount',
            'tumourDepth',
            'tumourRefCount',
        ]:
            if row.get(col, ''):
                row[col] = int(row[col])

        # default depth to alt + ref if not given
        for sample_type in ['normal', 'rna', 'tumour']:
            if (
                row.get(f'{sample_type}RefCount', '')
                and row.get(f'{sample_type}AltCount', '')
                and not row.get(f'{sample_type}Depth', '')
            ):
                row[f'{sample_type}Depth'] = (
                    row[f'{sample_type}RefCount'] + row[f'{sample_type}AltCount']
                )

    return result


def preprocess_expression_variants(rows: Iterable[Dict]) -> List[IprGeneVariant]:
    """
    Validate the input rows contain the minimum required fields and
    generate any default values where possible
    """

    def row_key(row: Dict) -> Tuple[str, ...]:
        return tuple(['expression'] + [row[key] for key in EXP_KEY])

    result = validate_variant_rows(rows, EXP_REQ, EXP_OPTIONAL, row_key)

    patterns = {}

    float_columns = [
        col
        for col in EXP_REQ + EXP_OPTIONAL
        if col.endswith('kIQR')
        or col.endswith('Percentile')
        or col.endswith('FoldChange')
        or col.endswith('QC')
        or col.endswith('ZScore')
        or col in ['tpm', 'rpkm']
    ]
    for col in float_columns:
        if col not in patterns:
            patterns[col] = NULLABLE_FLOAT_REGEX
    validate_row_patterns(result, patterns, EXP_KEY)

    errors = []
    for row in result:
        row['variant'] = row['kbCategory']
        if not row['expressionState'] and row['kbCategory']:
            row['expressionState'] = row['kbCategory']

        if row['variant']:
            if row['variant'] not in INPUT_EXPRESSION_CATEGORIES.values():
                err_msg = f"{row['gene']} variant '{row['variant']}' not in {INPUT_EXPRESSION_CATEGORIES.values()}"
                errors.append(err_msg)
                logger.error(err_msg)
        row['variantType'] = 'exp'

        for col in float_columns:
            if row[col] in ['inf', '+inf', '-inf']:
                row[col] = row[col].replace('inf', 'Infinity')

        # check images exist
        if row['histogramImage'] and not os.path.exists(row['histogramImage']):
            raise FileNotFoundError(f'missing image ({row["histogramImage"]})')

    if errors:
        raise ValueError(f'{len(errors)} Invalid expression variants in file')

    return result


def create_graphkb_sv_notation(row: IprStructuralVariant) -> str:
    """
    Generate GKB style structural variant notation from a structural variant input row
    """
    gene1 = row['gene1'] if row['gene1'] else '?'
    gene2 = row['gene2'] if row['gene2'] else '?'
    exon1 = row['exon1'] if row['exon1'] else '?'
    exon2 = row['exon2'] if row['exon2'] else '?'
    if not row['gene1']:
        gene1, gene2 = gene2, gene1
        exon1, exon2 = exon2, exon1
    if gene1 == '?':
        raise ValueError(
            f'both genes cannot be blank for a structural variant {row["key"]}. At least 1 gene must be entered'
        )
    return f'({gene1},{gene2}):fusion(e.{exon1},e.{exon2})'


def preprocess_structural_variants(rows: Iterable[Dict]) -> List[IprVariant]:
    """
    Validate the input rows contain the minimum required fields and
    generate any default values where possible
    """

    def row_key(row: Dict) -> Tuple[str, ...]:
        return tuple(['sv'] + [row[key] for key in SV_KEY])

    result = validate_variant_rows(rows, SV_REQ, SV_OPTIONAL, row_key)
    # genes are optional for structural variants
    EXON_PATTERN = r'^(\d+)?$'
    patterns = {
        'gene1': r'^((\w|-)+)?$',
        'gene2': r'^((\w|-)+)?$',
        'breakpoint': r'^(\w+:\d+\|\w+:\d+)?$',
        'exon1': EXON_PATTERN,
        'exon2': EXON_PATTERN,
    }
    validate_row_patterns(result, patterns, SV_KEY)

    for row in result:
        row['variant'] = create_graphkb_sv_notation(row)
        row['variantType'] = 'sv'

        # check and load the svg file where applicable
        if row['svg']:
            if not os.path.exists(row['svg']):
                raise FileNotFoundError(row['svg'])
            with open(row['svg'], 'r') as fh:
                row['svg'] = fh.read()

        for bool_col in ['highQuality', 'omicSupport']:
            if row[bool_col]:
                if row[bool_col].lower() not in ['true', 'false']:
                    raise ValueError(f'{bool_col} flag must be true or false if given')
                row[bool_col] = bool(row[bool_col].lower() == 'true')

    return result


def check_variant_links(
    small_mutations: List[IprGeneVariant],
    expression_variants: List[IprGeneVariant],
    copy_variants: List[IprGeneVariant],
    structural_variants: List[IprStructuralVariant],
) -> Set[str]:
    """
    Check matching information for any genes with variants.
    Warn about genes with only one experimental measure.

    Args:
        small_mutations: list of small mutations
        expression_variants: list of expression variants
        copy_variants: list of copy variants
        structural_variants: list of structural variants

    Returns:
        set of gene names with variants (used for filtering before upload to IPR)
    """
    # filter excess variants not required for extra gene information
    missing_information_genes = set()
    missing_information_errors = set()

    copy_variant_genes = {variant['gene'] for variant in copy_variants}
    expression_variant_genes = {variant['gene'] for variant in expression_variants}
    genes_with_variants = set()  # filter excess copy variants

    for variant in copy_variants:
        gene = variant['gene']
        if not gene:
            logger.error("copy_variant data cannot be applied to an empty genename")
        elif variant['variant']:
            genes_with_variants.add(gene)

            if expression_variant_genes and gene not in expression_variant_genes:
                missing_information_genes.add(gene)
                missing_information_errors.add(
                    f'gene ({gene}) has a copy variant but is missing expression information'
                )

    for variant in expression_variants:
        gene = variant['gene']
        if not gene:
            logger.error("expression_variant data cannot be applied to an empty genename")
        elif variant['variant']:
            genes_with_variants.add(gene)

            if copy_variant_genes and gene not in copy_variant_genes:
                missing_information_genes.add(gene)
                missing_information_errors.add(
                    f'gene ({gene}) has an expression variant but is missing copy number information'
                )

    for variant in small_mutations:
        gene = variant['gene']
        if not gene:
            logger.error("small_mutation data cannot be applied to an empty genename")
            continue

        if copy_variant_genes and gene not in copy_variant_genes:
            missing_information_genes.add(gene)
            missing_information_errors.add(
                f'gene ({gene}) has a small mutation but is missing copy number information'
            )
        if expression_variant_genes and gene not in expression_variant_genes:
            missing_information_genes.add(gene)
            missing_information_errors.add(
                f'gene ({gene}) has a small mutation but is missing expression information'
            )
        genes_with_variants.add(gene)

    for variant in structural_variants:
        for gene in [variant['gene1'], variant['gene2']]:
            if gene:  # genes are optional for structural variants
                if gene not in copy_variant_genes:
                    missing_information_genes.add(gene)
                    missing_information_errors.add(
                        f'gene ({gene}) has a structural variant but is missing copy number information'
                    )
                if gene not in expression_variant_genes:
                    missing_information_genes.add(gene)
                    missing_information_errors.add(
                        f'gene ({gene}) has a structural variant but is missing expression information'
                    )
                genes_with_variants.add(gene)

    if missing_information_genes:
        for err_msg in sorted(missing_information_errors):
            logger.verbose(err_msg)  # type: ignore
        link_err_msg = (
            f'Missing information variant links on {len(missing_information_genes)} genes'
        )
        logger.warning(link_err_msg)
    return genes_with_variants


def check_comparators(content: Dict, expresssionVariants: Iterable[Dict] = []) -> None:
    """
    Given the optional content dictionary, check that based on the analyses present the
    correct/sufficient comparators have also been specified
    """
    mutation_burden = 'mutationBurden'
    comparator_roles = {c['analysisRole'] for c in content.get('comparators', [])}

    for image in content.get('images', []):
        key = image['key']
        if key.startswith(mutation_burden):
            comp_type = key.split('.')[-1]
            role = f'mutation burden ({comp_type})'
            if role in comparator_roles:
                continue
            if '_sv.' in key:
                sv_role = f'mutation burden SV ({comp_type})'
                if sv_role in comparator_roles:
                    continue
            raise ValueError(f'missing required comparator definition ({role})')

    if expresssionVariants:
        required_comparators = {'expression (disease)'}

        def all_none(row: Dict, columns: List[str]) -> bool:
            return all([row.get(col) is None for col in columns])

        for exp in expresssionVariants:
            if not all_none(
                exp,
                [
                    'primarySitekIQR',
                    'primarySitePercentile',
                    'primarySiteZScore',
                    'primarySiteFoldChange',
                ],
            ):
                required_comparators.add('expression (primary site)')

            if not all_none(
                exp,
                [
                    'biopsySitekIQR',
                    'biopsySitePercentile',
                    'biopsySiteZScore',
                    'biopsySiteFoldChange',
                ],
            ):
                required_comparators.add('expression (biopsy site)')

        if required_comparators - comparator_roles:
            missing = '; '.join(sorted(list(required_comparators - comparator_roles)))
            raise ValueError(f'missing required comparator definitions ({missing})')
