"""
Read/Validate the variant input files
"""
from csv import DictReader
import re

from graphkb.match import INPUT_COPY_CATEGORIES

from .util import hash_key


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
                    f'row value ({row[col]}) does not match expected column ({col}) pattern'
                )


def load_copy_variants(filename):
    def row_key(row):
        return ('cnv', row['gene'])

    result = load_variant_file(
        filename,
        ['gene', 'cnvState'],
        ['chromosomeBand', 'ploidyCorrCpChange', 'start', 'end', 'lohState'],
        row_key,
    )

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

    result = load_variant_file(
        filename,
        ['location', 'refAlt', 'gene', 'proteinChange', 'transcript'],
        ['zygosity', 'tumourReads', 'RNAReads'],
        row_key,
    )

    patterns = {'location': r'^\w+:\d+$', 'refAlt': r'^[A-Z]+>[A-Z]+$'}

    validate_row_patterns(result, patterns)

    return result


def load_expression_variants(filename):
    def row_key(row):
        return ('expression', row['gene'])

    result = load_variant_file(
        filename,
        ['gene', 'expression_class'],
        [
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
        ],
        row_key,
    )

    return result


def load_structural_variants(filename):
    def row_key(row):
        return ('sv', row['eventType'], row['breakpoint'])

    result = load_variant_file(
        filename,
        [
            'ctermGene',
            'ntermGene',
            'ctermTranscript',
            'ntermTranscript',
            'exons',
            'eventType',
            'genes',
            'breakpoint',
        ],
        ['detectedIn', 'conventionalName', 'svg', 'svgTitle', 'name', 'frame'],
        row_key,
    )
    patterns = {
        'genes': r'^\w+::\w+$',
        'breakpoint': r'^\w+:\d+\|\w+:\d+$',
        'exons': r'^e\d+:e\d+$',
    }
    validate_row_patterns(result, patterns)

    return result
