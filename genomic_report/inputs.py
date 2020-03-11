"""
Read/Validate the variant input files
"""
from csv import DictReader
import re

from graphkb.match import INPUT_COPY_CATEGORIES, INPUT_EXPRESSION_CATEGORIES
from Bio.Data.IUPACData import protein_letters_3to1

from .util import hash_key

protein_letters_3to1.setdefault('Ter', '*')


COPY_REQ = ['gene', 'variant']  # 'variant' in INPUT_COPY_CATEGORIES
COPY_OPTIONAL = [
    'ploidyCorrCpChange',
    'lohState',  # Loss of Heterzygosity state - informative detail to analyst
    'chromosomeBand',
    'start',
    'end',
]

SMALL_MUT_REQ = ['location', 'refAlt', 'gene', 'proteinChange', 'transcript']
SMALL_MUT_OPTIONAL = ['zygosity', 'tumourReads', 'RNAReads']

EXP_REQ = ['gene', 'expression_class']
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
    'ptsPogPerc',
    'gtexComp',
    'gtexFC',
    'gtexkIQR',
    'gtexAvgPerc',
    'gtexAvgFC',
    'gtexAvgkIQR',
]

SV_REQ = [
    'ctermGene',
    'ntermGene',
    'ctermTranscript',
    'ntermTranscript',
    'exons',
    'eventType',
    'genes',
    'breakpoint',
]
SV_OPTIONAL = ['detectedIn', 'conventionalName', 'svg', 'svgTitle', 'name', 'frame']


def load_variant_file(filename, required, optional, row_to_key):
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
                        raise ValueError(f'header missing required column ({req_col})')
                header_validated = True
            row_key = hash_key(row_to_key(row))
            if row_key in keys:
                raise ValueError(f'duplicate row key ({row_key})')
            row['key'] = row_key
            keys.add(row_key)

            result.append({col: row.get(col, '') for col in header})

    return result


def validate_row_patterns(rows, patterns):
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
            if not re.match(pattern, row[col]):
                raise ValueError(
                    f'row value ({row[col]}) does not match expected column ({col}) pattern of "{pattern}"'
                )


def load_copy_variants(filename):
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
        elif 'ploidyCorrCpChange' in row.keys() and row['ploidyCorrCpChange'] not in ('', 'na'):
            row['cnvState'] = 'Neutral'
        else:
            row['cnvState'] = ''  # no measurement

    return result


def load_small_mutations(filename):
    def row_key(row):
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

    return result


def load_expression_variants(filename):
    def row_key(row):
        return ('expression', row['gene'])

    result = load_variant_file(filename, EXP_REQ, EXP_OPTIONAL, row_key,)
    patterns = {
        'expression_class': r'^(overexpressed|outlier_high|high_percentile|outlier_low|low_percentile|underexpressed|no_category|na|)$'
    }

    validate_row_patterns(result, patterns)

    # transform expression class to 'variant' column of GraphKB vocabulary
    for row in result:
        variant = ''

        if row['expression_class'] in {'outlier_low', 'underexpressed', 'low_percentile'}:
            variant = INPUT_EXPRESSION_CATEGORIES.DOWN
        elif row['expression_class'] in {'outlier_high', 'overexpressed', 'high_percentile'}:
            variant = INPUT_EXPRESSION_CATEGORIES.UP

        row['variant'] = variant

    return result


def load_structural_variants(filename):
    def row_key(row):
        return ('sv', row['eventType'], row['breakpoint'])

    result = load_variant_file(filename, SV_REQ, SV_OPTIONAL, row_key)
    patterns = {
        'genes': r'^(\w|-)+::(\w|-)+$',
        'breakpoint': r'^\w+:\d+\|\w+:\d+$',
        # "e:e" just means no exon data.
        'exons': r'^e(\d+)?:e(\d+)?$',
    }
    validate_row_patterns(result, patterns)

    for row in result:
        exon1, exon2 = re.match(r'^e(\d+)?:e(\d+)?$', row['exons']).group(1, 2)
        gene1, gene2 = row['genes'].split('::')

        row[
            'variant'
        ] = f'({gene1},{gene2}):fusion(e.{exon1 if exon1 else "?"},e.{exon2 if exon2 else "?"})'
        del row['genes']
        del row['exons']
        row['gene1'] = gene1
        row['gene2'] = gene2
        row['exon1'] = exon1
        row['exon2'] = exon2

    return result


def check_variant_links(small_mutations, expression_variants, copy_variants, structural_variants):
    """
    Check that there is matching expression and copy variant information for any genes with variants

    Args:
        small_mutations (list.<dict>): list of small mutations
        expression_variants (list.<dict>): list of expression variants
        copy_variants (list.<dict>): list of copy variants
        structural_variants (list.<dict>): list of structural variants

    Raises:
        KeyError: A variant is called on a gene without expression or without copy number information

    Returns:
        set.<str>: set of gene names with variants (used for filtering before upload to IPR)
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
