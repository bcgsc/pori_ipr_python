"""
Read/Validate the variant input files
"""
from typing import List, Dict, Tuple, Set
from csv import DictReader
import logging
import re
import os

from graphkb.match import INPUT_COPY_CATEGORIES, INPUT_EXPRESSION_CATEGORIES
from Bio.Data.IUPACData import protein_letters_3to1

from .util import hash_key

protein_letters_3to1.setdefault('Ter', '*')

NULLABLE_FLOAT_REGEX = r'^-?((inf)|(\d+(\.\d+)?)|)$'
COPY_REQ = ['gene', 'variant']  # 'variant' in INPUT_COPY_CATEGORIES
COPY_OPTIONAL = [
    'ploidyCorrCpChange',
    'lohState',  # Loss of Heterzygosity state - informative detail to analyst
    'chromosomeBand',
    'start',
    'end',
]

SMALL_MUT_REQ = ['location', 'refAlt', 'gene', 'proteinChange', 'transcript']
SMALL_MUT_OPTIONAL = ['zygosity', 'tumourReads', 'rnaReads', 'detectedIn']

# 'expression_class', is for display
EXP_REQ = ['gene', 'variant', 'expression_class']
EXP_OPTIONAL = [
    'rnaReads',
    'rpkm',
    'foldChange',
    'tcgaPerc',
    'tcgaPercCol',
    'tcgakIQR',
    'tcgaQC',
    'tcgaQCCol',
    'tcgaAvgPerc',
    'tcgaAvgkIQR',
    'tcgaAvgQC',
    'tcgaAvgQCCol',
    'tcgaNormPerc',
    'tcgaNormkIQR',
    'ptxPerc',
    'ptxkIQR',
    'ptxQC',
    'ptxPercCol',
    'ptxTotSampObs',
    'ptxPogPerc',
    'gtexComp',
    'gtexFC',
    'gtexkIQR',
    'gtexAvgPerc',
    'gtexAvgFC',
    'gtexAvgkIQR',
]

SV_KEY = ['eventType', 'breakpoint', 'gene1', 'gene2', 'exon1', 'exon2']
SV_REQ = [
    'eventType',
    'breakpoint',
    'gene1',  # prev: nterm_hugo
    'gene2',  # prev: cterm_hugo
    'exon1',  # n-terminal
    'exon2',  # c-terminal
]
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
]


def load_variant_file(
    filename: str, required: List[str], optional: List[str], row_to_key
) -> List[Dict]:
    """
    Load a tab delimited file and
    - check that the required columns are present
    - check that a unique key can be formed for each row
    - drop any non-defined columns

    Args:
        filename (str): the file to be read
        required (List.<str>): list of required column names
        optional (List.<str>): list of optional column names
        row_to_key (Function): function to generate a key for a given row

    Raises:
        ValueError: row keys are not unique
        ValueError: A required column is missing

    Returns:
        List.<dict>: the rows from the tab file as dictionaries
    """
    header = required + optional + ['key']

    result = []
    keys = set()

    with open(filename, 'r') as fh:
        reader = DictReader(fh, delimiter='\t')
        header_validated = False

        for row in reader:
            if not header_validated:
                for req_col in required:
                    if req_col not in row:
                        raise ValueError(
                            f'header missing required column ({req_col}) in {filename}'
                        )
                header_validated = True
            row_key = hash_key(row_to_key(row))
            if row_key in keys:
                raise ValueError(
                    f'duplicate row key ({row_key}) from ({row_to_key(row)}) in {filename}'
                )
            row['key'] = row_key
            keys.add(row_key)

            result.append({col: row.get(col, '') for col in header})

    return result


def validate_row_patterns(rows: List[Dict], patterns: Dict):
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
            if not re.match(pattern, '' if row[col] is None else row[col]):
                raise ValueError(
                    f'row value ({row[col]}) does not match expected column ({col}) pattern of "{pattern}"'
                )


def load_copy_variants(filename: str) -> List[Dict]:
    # default map for display - concise names
    COPY_VARIANT2CNVSTATE = {
        INPUT_COPY_CATEGORIES.DEEP: "Deep Loss",
        INPUT_COPY_CATEGORIES.AMP: "Amplification",
        INPUT_COPY_CATEGORIES.GAIN: "Gain",
        INPUT_COPY_CATEGORIES.LOSS: "Loss",
    }

    def row_key(row):
        return ('cnv', row['gene'])

    result = load_variant_file(filename, COPY_REQ, COPY_OPTIONAL, row_key)

    # verify the copy number category is valid or blank
    patterns = {'variant': f'({"|".join(INPUT_COPY_CATEGORIES.values())}|)'}
    validate_row_patterns(result, patterns)

    # Create a 'cnvState' displayed variant label
    for row in result:
        if row['variant'] in COPY_VARIANT2CNVSTATE:
            row['cnvState'] = COPY_VARIANT2CNVSTATE[row['variant']]
        # any non-blank measurement, without another category, is Neutral.
        else:
            row['cnvState'] = ''  # no measurement

        row['variantType'] = 'cnv'

    return result


def load_small_mutations(filename: str) -> List[Dict]:
    def row_key(row: Dict) -> Tuple[str]:
        return (
            'small mutation',
            row['location'],
            row['refAlt'],
            row['gene'],
            row['proteinChange'],
            row['transcript'],
        )

    result = load_variant_file(filename, SMALL_MUT_REQ, SMALL_MUT_OPTIONAL, row_key,)

    patterns = {'location': r'^\w+:\d+$', 'refAlt': r'^[A-Z]+>[A-Z]+$'}

    validate_row_patterns(result, patterns)

    # change 3 letter AA to 1 letter AA notation
    for row in result:
        for longAA, shortAA in protein_letters_3to1.items():
            row['proteinChange'] = row['proteinChange'].replace(longAA, shortAA)
        hgvsp = '{}:{}'.format(row['gene'], row['proteinChange'])
        row['variant'] = hgvsp
        row['variantType'] = 'mut'

    return result


def load_expression_variants(filename):
    def row_key(row):
        return ('expression', row['gene'])

    result = load_variant_file(filename, EXP_REQ, EXP_OPTIONAL, row_key)

    patterns = {}

    float_columns = [
        col
        for col in EXP_REQ + EXP_OPTIONAL
        if col.endswith('kIQR') or col.endswith('Perc') or col.endswith('FC')
    ]
    for col in float_columns:
        if col not in patterns:
            patterns[col] = NULLABLE_FLOAT_REGEX
    validate_row_patterns(result, patterns)

    errors = []
    for row in result:
        if row['variant']:
            if row['variant'] not in INPUT_EXPRESSION_CATEGORIES.values():
                err_msg = (
                    f"{row['gene']} variant '{row['variant']}' not in {INPUT_EXPRESSION_CATEGORIES}"
                )
                errors.append(err_msg)
                logging.error(err_msg)
        else:
            row['expression_class'] = ''
        row['variantType'] = 'exp'

        for col in float_columns:
            if row[col] in ['inf', '+inf', '-inf']:
                row[col] = row[col].replace('inf', 'Infinity')

    if errors:
        raise ValueError(f"{len(errors)} Invalid expression variants in file - {filename}")

    return result


def load_structural_variants(filename: str) -> List[Dict]:
    def row_key(row):
        return ('sv', row['eventType'], row['breakpoint'])

    result = load_variant_file(filename, SV_REQ, SV_OPTIONAL, row_key)
    exon_pattern = r'^(\d+)?$'
    patterns = {
        'gene1': r'^(\w|-)+$',
        'gene2': r'^(\w|-)+$',
        'breakpoint': r'^\w+:\d+\|\w+:\d+$',
        'exon1': exon_pattern,
        'exon2': exon_pattern,
    }
    validate_row_patterns(result, patterns)

    for row in result:
        row[
            'variant'
        ] = f'({row["gene1"]},{row["gene2"]}):fusion(e.{row["exon1"]},e.{row["exon2"]})'
        row['variantType'] = 'sv'

        # check and load the svg file where applicable
        if row['svg']:
            if not os.path.exists(row['svg']):
                raise FileNotFoundError(row['svg'])
            with open(row['svg'], 'r') as fh:
                row['svg'] = fh.read()

    return result


def check_variant_links(
    small_mutations: List[Dict],
    expression_variants: List[Dict],
    copy_variants: List[Dict],
    structural_variants: List[Dict],
) -> Set[str]:
    """
    Check that there is matching expression and copy variant information for any genes with variants

    Args:
        small_mutations: list of small mutations
        expression_variants: list of expression variants
        copy_variants: list of copy variants
        structural_variants: list of structural variants

    Raises:
        KeyError: A variant is called on a gene without expression or without copy number information

    Returns:
        set of gene names with variants (used for filtering before upload to IPR)
    """
    # filter excess variants not required for extra gene information
    copy_variant_genes = {variant['gene'] for variant in copy_variants}
    expression_variant_genes = {variant['gene'] for variant in expression_variants}
    genes_with_variants = set()  # filter excess copy variants

    for variant in copy_variants:
        gene = variant['gene']
        if variant['variant']:
            genes_with_variants.add(gene)

            if expression_variant_genes and gene not in expression_variant_genes:
                raise KeyError(
                    f'gene ({gene}) has a copy variant but is missing expression information'
                )

    for variant in expression_variants:
        gene = variant['gene']
        if variant['variant']:
            genes_with_variants.add(gene)

            if copy_variant_genes and gene not in copy_variant_genes:
                raise KeyError(
                    f'gene ({gene}) has an expression variant but is missing copy number information'
                )

    for variant in small_mutations:
        gene = variant['gene']
        if copy_variant_genes and gene not in copy_variant_genes:
            raise KeyError(
                f'gene ({gene}) has a small mutation but is missing copy number information'
            )
        if expression_variant_genes and gene not in expression_variant_genes:
            raise KeyError(
                f'gene ({gene}) has a small mutation but is missing expression information'
            )
        genes_with_variants.add(gene)

    for variant in structural_variants:
        for gene in [variant['gene1'], variant['gene2']]:
            if gene:  # genes are optional for structural variants
                if gene not in copy_variant_genes:
                    raise KeyError(
                        f'gene ({gene}) has a structural variant but is missing copy number information'
                    )
                if gene not in expression_variant_genes:
                    raise KeyError(
                        f'gene ({gene}) has a structural variant but is missing expression information'
                    )
                genes_with_variants.add(gene)

    return genes_with_variants
