"""
Read/Validate the variant input files
"""
import json
import jsonschema
import os
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
from graphkb.match import INPUT_COPY_CATEGORIES, INPUT_EXPRESSION_CATEGORIES
from typing import Callable, Dict, Iterable, List, Set, Tuple, cast

from .types import (
    IprCopyVariant,
    IprExprVariant,
    IprFusionVariant,
    IprSmallMutationVariant,
    IprVariant,
)
from .util import hash_key, logger, pandas_falsy

protein_letters_3to1.setdefault('Ter', '*')

SPECIFICATION = os.path.join(os.path.dirname(__file__), 'content.spec.json')
# content in the local specification should match the values in IPR_API_SPEC_JSON_URL
IPR_API_SPEC_JSON_URL = 'https://ipr-api.bcgsc.ca/api/spec.json'


# TODO: GERO-307 - use SPECIFICATION json to derive the variant required and optional details defined below

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
    'comments',
    'library',
    'germline',
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
    'comments',
    'chromosome',
    'endPosition',
    'germline',
    'hgvsCds',
    'hgvsGenomic',
    'hgvsProtein',
    'library',
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
    'tumourAltCopies',
    'tumourDepth',
    'tumourRefCount',
    'tumourRefCopies',
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
    'comments',
    'diseaseFoldChange',
    'diseasekIQR',
    'diseasePercentile',
    'diseaseQC',
    'diseaseZScore',
    'expressionState',
    'histogramImage',
    'library',
    'primarySiteFoldChange',
    'primarySitekIQR',
    'primarySitePercentile',
    'primarySiteQC',
    'primarySiteZScore',
    'internalPancancerFoldChange',
    'internalPancancerkIQR',
    'internalPancancerPercentile',
    'internalPancancerQC',
    'internalPancancerZScore',
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
    'comments',
    'library',
    # GERO-307 - tumourAltCount and tumourDepth are available but not rnaAltCount and rnaDepth
    'tumourAltCount',
    'tumourDepth',
    'germline',
    'mavis_product_id',
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
        for k, v in row.items():
            if v is pd.NA:
                row[k] = ''

        result.append(cast(IprVariant, {col: row.get(col, '') for col in header}))

    return result


def preprocess_copy_variants(rows: Iterable[Dict]) -> List[IprCopyVariant]:
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
    display_name_mapping.update(dict([(v, v) for v in display_name_mapping.values()]))

    def row_key(row: Dict) -> Tuple[str, ...]:
        return tuple(['cnv'] + [row[key] for key in COPY_KEY])

    result = validate_variant_rows(rows, COPY_REQ, COPY_OPTIONAL, row_key)
    ret_list = [cast(IprCopyVariant, var) for var in result]
    for row in ret_list:
        kb_cat = row.get('kbCategory')
        kb_cat = '' if pd.isnull(kb_cat) else str(kb_cat)
        if kb_cat:
            if kb_cat not in INPUT_COPY_CATEGORIES.values():
                raise ValueError(f'invalid copy variant kbCategory value ({kb_cat})')
            if not row.get('cnvState'):  # apply default short display name
                row['cnvState'] = display_name_mapping[kb_cat]
        row['variant'] = kb_cat
        row['variantType'] = 'cnv'

    return ret_list


def preprocess_small_mutations(rows: Iterable[Dict]) -> List[IprSmallMutationVariant]:
    """
    Validate the input rows contain the minimum required fields and
    generate any default values where possible
    """

    def row_key(row: IprSmallMutationVariant) -> Tuple[str, ...]:
        key_vals = []
        for kval in [row.get(key, '') for key in SMALL_MUT_KEY]:
            key_vals.append(str(kval) if pd.notnull(kval) else '')
        return tuple(['small mutation'] + key_vals)

    result = validate_variant_rows(rows, SMALL_MUT_REQ, SMALL_MUT_OPTIONAL, row_key)
    if not result:
        return []

    def pick_variant(row: IprSmallMutationVariant) -> str:
        protein_change = row.get('proteinChange')
        if not pandas_falsy(protein_change):
            for longAA, shortAA in protein_letters_3to1.items():
                protein_change = str(protein_change).replace(longAA, shortAA)
            hgvsp = '{}:{}'.format(row['gene'], protein_change)
            return hgvsp

        for field in ['hgvsProtein', 'hgvsCds', 'hgvsGenomic']:
            if not pandas_falsy(row.get(field)):
                return str(row.get(field))

        raise ValueError(
            'Variant field cannot be empty. Must include proteinChange or one of the hgvs fields (hgvsProtein, hgvsCds, hgvsGenomic) to build the variant string'
        )

    # 'location' and 'refAlt' are not currently used for matching; still optional and allowed blank

    # change 3 letter AA to 1 letter AA notation
    # for row in result:
    def convert_sm(row: IprVariant) -> IprSmallMutationVariant:
        ret = cast(IprSmallMutationVariant, row)
        ret['variant'] = pick_variant(ret)
        ret['variantType'] = 'mut'

        if ret.get('startPosition') and not ret.get('endPosition'):
            ret['endPosition'] = ret['startPosition']

        # default depth to alt + ref if not given
        for sample_type in ('normal', 'rna', 'tumour'):
            if (
                ret.get(f'{sample_type}RefCount')
                and ret.get(f'{sample_type}AltCount')
                and not ret.get(f'{sample_type}Depth')
            ):
                ret[f'{sample_type}Depth'] = (  # type: ignore
                    ret[f'{sample_type}RefCount'] + ret[f'{sample_type}AltCount']  # type: ignore
                )
        return ret

    res_list = [convert_sm(var) for var in result]

    return res_list


def preprocess_expression_variants(rows: Iterable[Dict]) -> List[IprExprVariant]:
    """
    Validate the input rows contain the minimum required fields and
    generate any default values where possible
    """

    def row_key(row: Dict) -> Tuple[str, ...]:
        return tuple(['expression'] + [row[key] for key in EXP_KEY])

    variants = validate_variant_rows(rows, EXP_REQ, EXP_OPTIONAL, row_key)
    result = [cast(IprExprVariant, var) for var in variants]
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

    errors = []
    for row in result:
        row['variant'] = row['kbCategory']
        if not row['expressionState'] and row['kbCategory']:
            row['expressionState'] = row['kbCategory']

        if row['variant'] and not pd.isnull(row['variant']):
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


def create_graphkb_sv_notation(row: IprFusionVariant) -> str:
    """Generate GKB/IPR fusion style notation from a structural variant."""
    gene1 = row['gene1'] or '?'
    gene2 = row['gene2'] or '?'
    exon1 = str(row['exon1']) if row['exon1'] else '?'
    exon2 = str(row['exon2']) if row['exon2'] else '?'
    if not row['gene1']:
        gene1, gene2 = gene2, gene1
        exon1, exon2 = exon2, exon1
    if gene1 == '?':
        raise ValueError(
            f'both genes cannot be blank for a structural variant {row["key"]}. At least 1 gene must be entered'
        )
    # force exons to integer repr string
    exon1 = exon1[:-2] if exon1.endswith(".0") else exon1
    exon2 = exon2[:-2] if exon2.endswith(".0") else exon2
    return f'({gene1},{gene2}):fusion(e.{exon1},e.{exon2})'


def preprocess_structural_variants(rows: Iterable[Dict]) -> List[IprFusionVariant]:
    """
    Validate the input rows contain the minimum required fields and
    generate any default values where possible
    """

    def row_key(row: Dict) -> Tuple[str, ...]:
        return tuple(['sv'] + [row[key] for key in SV_KEY])

    variants = validate_variant_rows(rows, SV_REQ, SV_OPTIONAL, row_key)
    result = [cast(IprFusionVariant, var) for var in variants]
    # genes are optional for structural variants
    for row in result:
        row['variant'] = create_graphkb_sv_notation(row)
        row['variantType'] = 'sv'

        # check and load the svg file where applicable
        if row['svg'] and not pd.isnull(row['svg']):
            if not os.path.exists(row['svg']):
                raise FileNotFoundError(row['svg'])
            with open(row['svg'], 'r') as fh:
                row['svg'] = fh.read()

    return result


def check_variant_links(
    small_mutations: List[IprSmallMutationVariant],
    expression_variants: List[IprExprVariant],
    copy_variants: List[IprCopyVariant],
    structural_variants: List[IprFusionVariant],
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

    variant = IprVariant  # to silence type errors
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
            logger.debug(err_msg)
        link_err_msg = (
            f'Missing information variant links on {len(missing_information_genes)} genes'
        )
        logger.warning(link_err_msg)
    return genes_with_variants


def check_comparators(content: Dict, expresssionVariants: List[IprExprVariant] = []) -> None:
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

        def all_none(row: IprExprVariant, columns: List[str]) -> bool:
            return all([row.get(col) is None or row.get(col) == '' for col in columns])

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

            if not all_none(
                exp,
                [
                    'internalPancancerkIQR',
                    'internalPancancerPercentile',
                    'internalPancancerZScore',
                    'internalPancancerFoldChange',
                ],
            ):
                required_comparators.add('expression (internal pancancer cohort)')

        if required_comparators - comparator_roles:
            missing = '; '.join(sorted(list(required_comparators - comparator_roles)))
            raise ValueError(f'missing required comparator definitions ({missing})')


def extend_with_default(validator_class):
    # https://python-jsonschema.readthedocs.io/en/latest/faq/#why-doesn-t-my-schema-s-default-property-set-the-default-on-my-instance
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_defaults(validator, properties, instance, schema):
        for property, subschema in properties.items():
            if "default" in subschema:
                instance.setdefault(property, subschema["default"])

        for error in validate_properties(validator, properties, instance, schema):
            yield error

    def check_null(checker, instance):
        return (
            validator_class.TYPE_CHECKER.is_type(instance, "null")
            or pd.isnull(instance)
            or instance == ""
        )

    type_checker = validator_class.TYPE_CHECKER.redefine("null", check_null)

    return jsonschema.validators.extend(
        validator_class, validators={"properties": set_defaults}, type_checker=type_checker
    )


# Customize the default jsonschema behaviour to add default values and treat np.nan as null
DefaultValidatingDraft7Validator = extend_with_default(jsonschema.Draft7Validator)


def validate_report_content(content: Dict, schema_file: str = SPECIFICATION) -> None:
    """
    Validate a report content input JSON object against the schema specification

    Adds defaults as reccommended by: https://python-jsonschema.readthedocs.io/en/latest/faq/#why-doesn-t-my-schema-s-default-property-set-the-default-on-my-instance
    """
    with open(schema_file, 'r') as fh:
        schema = json.load(fh)

    return DefaultValidatingDraft7Validator(schema).validate(content)
