import os
from unittest import mock

import pytest

from genomic_report.inputs import (
    load_small_mutations,
    load_copy_variants,
    check_variant_links,
    load_expression_variants,
    load_structural_variants,
)
from genomic_report.util import logger


DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')


def test_load_small_mutations():
    records = load_small_mutations(os.path.join(DATA_DIR, 'small_mutations.tab'))
    assert records
    assert len(records) == 2614


def test_load_copy_variants():
    records = load_copy_variants(os.path.join(DATA_DIR, 'copy_variants.tab'))
    assert records
    assert len(records) == 4603


def test_load_structural_variants():
    records = load_structural_variants(os.path.join(DATA_DIR, 'fusions.tab'))
    assert records
    assert len(records) == 3


def test_load_expression_variants():
    records = load_expression_variants(os.path.join(DATA_DIR, 'expression.tab'))
    assert records
    assert len(records) == 4603


class TestCheckVariantLinks:
    def test_sm_missing_copy_empty_ok(self):
        genes = check_variant_links(
            small_mutations=[{'gene': 'KRAS'}],
            copy_variants=[],
            expression_variants=[{'gene': 'KRAS', 'variant': ''}],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_sm_missing_exp_empty_ok(self):
        genes = check_variant_links(
            small_mutations=[{'gene': 'KRAS'}],
            copy_variants=[{'gene': 'KRAS', 'variant': ''}],
            expression_variants=[],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_sm_missing_copy(self):
        with mock.patch.object(logger, 'warning') as mock_debug:
            check_variant_links(
                small_mutations=[{'gene': 'KRAS'}],
                copy_variants=[{'gene': 'CDK', 'variant': ''}],
                expression_variants=[{'gene': 'KRAS', 'variant': ''}],
                structural_variants=[],
            )
            assert mock_debug.called

    def test_sm_missing_exp(self):
        with mock.patch.object(logger, 'warning') as mock_debug:
            check_variant_links(
                small_mutations=[{'gene': 'KRAS'}],
                copy_variants=[{'gene': 'KRAS', 'variant': ''}],
                expression_variants=[{'gene': 'CDK', 'variant': ''}],
                structural_variants=[],
            )
            assert mock_debug.called

    @pytest.mark.skip('TODO')
    def test_sv_missing_copy(self):
        with mock.patch.object(logger, 'warning') as mock_debug:
            pass

    @pytest.mark.skip('TODO')
    def test_sv_missing_exp(self):
        with mock.patch.object(logger, 'warning') as mock_debug:
            pass

    def test_with_valid_inputs(self):
        genes = check_variant_links(
            small_mutations=[{'gene': 'KRAS'}],
            copy_variants=[{'gene': 'KRAS', 'variant': ''}, {'gene': 'CDK', 'variant': ''}],
            expression_variants=[{'gene': 'KRAS', 'variant': ''}],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_copy_missing_exp(self):
        with mock.patch.object(logger, 'warning') as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[
                    {'gene': 'BRAF', 'variant': 'copy gain'},
                    {'gene': 'KRAS', 'variant': ''},
                ],
                expression_variants=[{'gene': 'KRAS', 'variant': ''}],
                structural_variants=[],
            )
            assert mock_debug.called

    def test_exp_missing_copy(self):
        with mock.patch.object(logger, 'warning') as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[{'gene': 'KRAS', 'variant': ''}],
                expression_variants=[{'gene': 'BRAF', 'variant': 'increased expression'}],
                structural_variants=[],
            )
            assert mock_debug.called
