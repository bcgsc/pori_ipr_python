import json
import numpy as np
import os
import pandas as pd
import pytest
from graphkb.match import INPUT_COPY_CATEGORIES
from unittest import mock

from ipr.inputs import (
    COPY_OPTIONAL,
    check_comparators,
    check_variant_links,
    create_graphkb_sv_notation,
    preprocess_copy_variants,
    preprocess_expression_variants,
    preprocess_small_mutations,
    preprocess_structural_variants,
    validate_report_content,
)
from ipr.types import IprFusionVariant, IprGeneVariant
from ipr.util import logger

DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')
NON_EMPTY_STRING_NULLS = ['', None, np.nan, pd.NA]


def read_data_file(filename):
    pass


class TestPreProcessSmallMutations:
    def test_load_test_file(self) -> None:
        records = preprocess_small_mutations(
            pd.read_csv(os.path.join(DATA_DIR, 'small_mutations.tab'), sep='\t').to_dict('records')
        )
        assert records
        assert len(records) == 2614

    def test_maintains_optional_fields(self):
        original = {
            'gene': 'A1BG',
            'proteinChange': 'p.V460M',
            'zygosity': 'het',
            'tumourAltCount': 42,
            'tumourRefCount': 48,
            'hgvsProtein': '',
            'transcript': 'ENST1000',
            'hgvsCds': '',
            'hgvsGenomic': '',
            'key': '02fe85a3477784b5ac0f8ecffb300d10',
            'variant': 'blargh',
            'chromosome': '2',
            'startPosition': 1234,
        }
        records = preprocess_small_mutations([original])
        record = records[0]
        assert record['variantType'] == 'mut'
        for col in original:
            assert col in record
        assert record['variant'] == 'A1BG:p.V460M'
        assert 'endPosition' in record
        assert record['endPosition'] == record['startPosition']
        assert 'tumourDepth' in record
        assert record['tumourDepth'] == 90

    def test_null(self):
        original = {
            'gene': 'A1BG',
            'proteinChange': 'p.V460M',
            'tumourAltCount': 42,
            'tumourRefCount': 48,
            'startPosition': 1234,
        }
        # Make sure TEST_KEYS are appropriate.
        # For some fields, like 'ref' and 'alt', NA is _not_ equivalent to a null string.
        TEST_KEYS = ['startPosition', 'endPosition', 'tumourAltCount', 'tumourRefCount']
        for key in TEST_KEYS:
            for null in NON_EMPTY_STRING_NULLS:
                small_mut = original.copy()
                small_mut[key] = null
                records = preprocess_small_mutations([small_mut])
                record = records[0]
                assert record['variantType'] == 'mut'
                for col in original:
                    assert col in record
                assert record['variant'] == 'A1BG:p.V460M'
                assert 'endPosition' in record

    def test_load_small_mutations_probe(self) -> None:
        records = preprocess_small_mutations(
            pd.read_csv(os.path.join(DATA_DIR, 'small_mutations_probe.tab'), sep='\t').to_dict(
                'records'
            )
        )
        assert records
        assert len(records) == 4
        assert records[0]['variantType'] == 'mut'
        assert 'variant' in records[0]


class TestPreProcessCopyVariants:
    def test_load_copy_variants(self) -> None:
        records = preprocess_copy_variants(
            pd.read_csv(os.path.join(DATA_DIR, 'copy_variants.tab'), sep='\t').to_dict('records')
        )
        assert records
        assert len(records) == 4603
        assert records[0]['variantType'] == 'cnv'
        assert 'variant' in records[0]

    def test_null(self):
        for kb_cat in list(INPUT_COPY_CATEGORIES.values()) + NON_EMPTY_STRING_NULLS:
            original = {'gene': 'ERBB2', 'kbCategory': kb_cat}
            for key in COPY_OPTIONAL:
                for null in NON_EMPTY_STRING_NULLS:
                    copy_var = original.copy()
                    copy_var[key] = null
                    records = preprocess_copy_variants([copy_var])
                    record = records[0]
                    assert record['variantType'] == 'cnv'


def test_load_structural_variants() -> None:
    records = preprocess_structural_variants(
        pd.read_csv(os.path.join(DATA_DIR, 'fusions.tab'), sep='\t').to_dict('records')
    )
    assert records
    assert len(records) == 5
    assert records[0]['variantType'] == 'sv'
    assert 'variant' in records[0]


def test_load_expression_variants() -> None:
    records = preprocess_expression_variants(
        pd.read_csv(os.path.join(DATA_DIR, 'expression.tab'), sep='\t').to_dict('records')
    )
    assert records
    assert len(records) == 4603
    assert records[0]['variantType'] == 'exp'
    assert 'variant' in records[0]


class TestCheckVariantLinks:
    def test_sm_missing_copy_empty_ok(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({'gene': 'KRAS'})],  # type: ignore
            copy_variants=[],
            expression_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],  # type: ignore
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_sm_missing_exp_empty_ok(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({'gene': 'KRAS'})],  # type: ignore
            copy_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],  # type: ignore
            expression_variants=[],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_sm_missing_copy(self) -> None:
        with mock.patch.object(logger, 'debug') as mock_debug:
            check_variant_links(
                small_mutations=[IprGeneVariant({'gene': 'KRAS'})],  # type: ignore
                copy_variants=[IprGeneVariant({'gene': 'CDK', 'variant': ''})],  # type: ignore
                expression_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],  # type: ignore
                structural_variants=[],
            )
            assert mock_debug.called

    def test_sm_missing_exp(self) -> None:
        with mock.patch.object(logger, 'debug') as mock_debug:
            check_variant_links(
                small_mutations=[IprGeneVariant({'gene': 'KRAS'})],  # type: ignore
                copy_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],  # type: ignore
                expression_variants=[IprGeneVariant({'gene': 'CDK', 'variant': ''})],  # type: ignore
                structural_variants=[],
            )
            assert mock_debug.called

    def test_with_valid_inputs(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({'gene': 'KRAS'})],  # type: ignore
            copy_variants=[
                IprGeneVariant({'gene': 'KRAS', 'variant': ''}),  # type: ignore
                IprGeneVariant({'gene': 'CDK', 'variant': ''}),  # type: ignore
            ],
            expression_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],  # type: ignore
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_copy_missing_exp(self) -> None:
        with mock.patch.object(logger, 'debug') as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[
                    IprGeneVariant({'gene': 'BRAF', 'variant': 'copy gain'}),  # type: ignore
                    IprGeneVariant({'gene': 'KRAS', 'variant': ''}),  # type: ignore
                ],
                expression_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],  # type: ignore
                structural_variants=[],
            )
            assert mock_debug.called

    def test_exp_missing_copy(self) -> None:
        with mock.patch.object(logger, 'debug') as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],  # type: ignore
                expression_variants=[
                    IprGeneVariant({'gene': 'BRAF', 'variant': 'increased expression'})  # type: ignore
                ],
                structural_variants=[],
            )
            assert mock_debug.called


class TestCreateGraphkbSvNotation:
    def test_both_genes_and_exons(self) -> None:
        notation = create_graphkb_sv_notation(
            IprFusionVariant({'gene1': 'A', 'gene2': 'B', 'exon1': 1, 'exon2': 2})  # type: ignore
        )
        assert notation == '(A,B):fusion(e.1,e.2)'

    def test_one_exon_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprFusionVariant({'gene1': 'A', 'gene2': 'B', 'exon1': '', 'exon2': 2})  # type: ignore
        )
        assert notation == '(A,B):fusion(e.?,e.2)'

    def test_one_gene_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprFusionVariant({'gene1': 'A', 'gene2': '', 'exon1': 1, 'exon2': 2})  # type: ignore
        )
        assert notation == '(A,?):fusion(e.1,e.2)'

    def test_first_gene_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprFusionVariant({'gene1': '', 'gene2': 'B', 'exon1': 1, 'exon2': 2})  # type: ignore
        )
        assert notation == '(B,?):fusion(e.2,e.1)'

    def test_no_genes_error(self) -> None:
        with pytest.raises(ValueError):
            create_graphkb_sv_notation(
                IprFusionVariant({'gene1': '', 'gene2': '', 'exon1': 1, 'exon2': 2, 'key': 'x'})  # type: ignore
            )


class TestCheckComparators:
    def test_missing_disease_expression_error(self):
        content = {'comparators': [{'analysisRole': 'expression (primary site)'}]}
        variants = [{}]

        with pytest.raises(ValueError):
            check_comparators(content, variants)

    def test_missing_primary_expression_error(self):
        content = {'comparators': [{'analysisRole': 'expression (disease)'}]}
        variants = [{'primarySiteFoldChange': 1}]

        with pytest.raises(ValueError):
            check_comparators(content, variants)

    def test_missing_biopsy_expression_error(self):
        content = {'comparators': [{'analysisRole': 'expression (disease)'}]}
        variants = [{'biopsySitePercentile': 1}]

        with pytest.raises(ValueError):
            check_comparators(content, variants)

    def test_expression_not_required_without_variants(self):
        content = {'comparators': []}
        variants = []

        assert check_comparators(content, variants) is None

    def test_missing_mutation_burden(self):
        content = {
            'comparators': [{'analysisRole': 'mutation burden (secondary)'}],
            'images': [{'key': 'mutationBurden.density_snv.primary'}],
        }
        variants = []

        with pytest.raises(ValueError):
            check_comparators(content, variants)


@pytest.mark.parametrize("example_name", ['no_variants', 'sm_and_exp', 'sm_only'])
def test_valid_json_inputs(example_name: str):
    with open(os.path.join(DATA_DIR, 'json_examples', f'{example_name}.json'), 'r') as fh:
        content = json.load(fh)
    validate_report_content(content)
