import os
from unittest import mock

import pytest

from genomic_report.inputs import (
    check_variant_links,
    create_graphkb_sv_notation,
    load_copy_variants,
    load_expression_variants,
    load_small_mutations,
    load_structural_variants,
)
from genomic_report.util import logger

DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')


def test_load_small_mutations() -> None:
    records = load_small_mutations(os.path.join(DATA_DIR, 'small_mutations.tab'))
    assert records
    assert len(records) == 2614


def test_load_copy_variants() -> None:
    records = load_copy_variants(os.path.join(DATA_DIR, 'copy_variants.tab'))
    assert records
    assert len(records) == 4603


def test_load_structural_variants() -> None:
    records = load_structural_variants(os.path.join(DATA_DIR, 'fusions.tab'))
    assert records
    assert len(records) == 3


def test_load_expression_variants() -> None:
    records = load_expression_variants(os.path.join(DATA_DIR, 'expression.tab'))
    assert records
    assert len(records) == 4603


class TestCheckVariantLinks:
    def test_sm_missing_copy_empty_ok(self) -> None:
        genes = check_variant_links(
            small_mutations=[{'gene': 'KRAS'}],
            copy_variants=[],
            expression_variants=[{'gene': 'KRAS', 'variant': ''}],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_sm_missing_exp_empty_ok(self) -> None:
        genes = check_variant_links(
            small_mutations=[{'gene': 'KRAS'}],
            copy_variants=[{'gene': 'KRAS', 'variant': ''}],
            expression_variants=[],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_sm_missing_copy(self) -> None:
        with mock.patch.object(logger, 'warning') as mock_debug:
            check_variant_links(
                small_mutations=[{'gene': 'KRAS'}],
                copy_variants=[{'gene': 'CDK', 'variant': ''}],
                expression_variants=[{'gene': 'KRAS', 'variant': ''}],
                structural_variants=[],
            )
            assert mock_debug.called

    def test_sm_missing_exp(self) -> None:
        with mock.patch.object(logger, 'warning') as mock_debug:
            check_variant_links(
                small_mutations=[{'gene': 'KRAS'}],
                copy_variants=[{'gene': 'KRAS', 'variant': ''}],
                expression_variants=[{'gene': 'CDK', 'variant': ''}],
                structural_variants=[],
            )
            assert mock_debug.called

    @pytest.mark.skip('TODO')
    def test_sv_missing_copy(self) -> None:
        with mock.patch.object(logger, 'warning') as mock_debug:
            pass

    @pytest.mark.skip('TODO')
    def test_sv_missing_exp(self) -> None:
        with mock.patch.object(logger, 'warning') as mock_debug:
            pass

    def test_with_valid_inputs(self) -> None:
        genes = check_variant_links(
            small_mutations=[{'gene': 'KRAS'}],
            copy_variants=[{'gene': 'KRAS', 'variant': ''}, {'gene': 'CDK', 'variant': ''}],
            expression_variants=[{'gene': 'KRAS', 'variant': ''}],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_copy_missing_exp(self) -> None:
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

    def test_exp_missing_copy(self) -> None:
        with mock.patch.object(logger, 'warning') as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[{'gene': 'KRAS', 'variant': ''}],
                expression_variants=[{'gene': 'BRAF', 'variant': 'increased expression'}],
                structural_variants=[],
            )
            assert mock_debug.called


class TestCreateGraphkbSvNotation:
    def test_both_genes_and_exons(self) -> None:
        notation = create_graphkb_sv_notation({'gene1': 'A', 'gene2': 'B', 'exon1': 1, 'exon2': 2})
        assert notation == '(A,B):fusion(e.1,e.2)'

    def test_one_exon_missing(self) -> None:
        notation = create_graphkb_sv_notation({'gene1': 'A', 'gene2': 'B', 'exon1': '', 'exon2': 2})
        assert notation == '(A,B):fusion(e.?,e.2)'

    def test_one_gene_missing(self) -> None:
        notation = create_graphkb_sv_notation({'gene1': 'A', 'gene2': '', 'exon1': 1, 'exon2': 2})
        assert notation == '(A,?):fusion(e.1,e.2)'

    def test_first_gene_missing(self) -> None:
        notation = create_graphkb_sv_notation({'gene1': '', 'gene2': 'B', 'exon1': 1, 'exon2': 2})
        assert notation == '(B,?):fusion(e.2,e.1)'

    def test_no_genes_error(self) -> None:
        with pytest.raises(ValueError):
            create_graphkb_sv_notation(
                {'gene1': '', 'gene2': '', 'exon1': 1, 'exon2': 2, 'key': 'x'}
            )
