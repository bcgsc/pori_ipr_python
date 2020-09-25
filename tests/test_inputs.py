import os
import pytest
from unittest import mock

from ipr.inputs import (
    check_comparators,
    check_variant_links,
    create_graphkb_sv_notation,
    preprocess_copy_variants,
    preprocess_expression_variants,
    preprocess_small_mutations,
    preprocess_structural_variants,
    read_tabbed_file,
)
from ipr.types import IprGeneVariant, IprStructuralVariant
from ipr.util import logger

DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')


class TestPreProcessSmallMutations:
    def test_load_test_file(self) -> None:
        records = preprocess_small_mutations(
            read_tabbed_file(os.path.join(DATA_DIR, 'small_mutations.tab'))
        )
        assert records
        assert len(records) == 2614

    def test_error_on_missing_gene(self):
        original = {
            'proteinChange': 'p.V460M',
            'zygosity': 'het',
            'tumourReads': '48/42',
            'rnaReads': '26/0',
            'hgvsProtein': '',
            'transcript': 'ENST1000',
            'hgvsCds': '',
            'hgvsGenomic': '',
            'key': '02fe85a3477784b5ac0f8ecffb300d10',
            'variant': 'A1BG:p.V460M',
            'location': '2:1234',
        }
        with pytest.raises(ValueError):
            preprocess_small_mutations([original])

    def test_error_on_missing_change(self):
        original = {
            'zygosity': 'het',
            'tumourReads': '48/42',
            'rnaReads': '26/0',
            'hgvsProtein': '',
            'transcript': 'ENST1000',
            'hgvsCds': '',
            'hgvsGenomic': '',
            'key': '02fe85a3477784b5ac0f8ecffb300d10',
            'location': '2:1234',
            'gene': 'KRAS',
        }
        with pytest.raises(ValueError):
            preprocess_small_mutations([original])

    def test_maintains_optional_fields(self):
        original = {
            'gene': 'A1BG',
            'proteinChange': 'p.V460M',
            'zygosity': 'het',
            'tumourAltCount': '42',
            'tumourRefCount': '48',
            'hgvsProtein': '',
            'transcript': 'ENST1000',
            'hgvsCds': '',
            'hgvsGenomic': '',
            'key': '02fe85a3477784b5ac0f8ecffb300d10',
            'variant': 'blargh',
            'chromosome': '2',
            'startPosition': '1234',
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


def test_load_small_mutations_probe() -> None:
    records = preprocess_small_mutations(
        read_tabbed_file(os.path.join(DATA_DIR, 'small_mutations_probe.tab'))
    )
    assert records
    assert len(records) == 4
    assert records[0]['variantType'] == 'mut'
    assert 'variant' in records[0]


def test_load_copy_variants() -> None:
    records = preprocess_copy_variants(
        read_tabbed_file(os.path.join(DATA_DIR, 'copy_variants.tab'))
    )
    assert records
    assert len(records) == 4603
    assert records[0]['variantType'] == 'cnv'
    assert 'variant' in records[0]


def test_load_structural_variants() -> None:
    records = preprocess_structural_variants(
        read_tabbed_file(os.path.join(DATA_DIR, 'fusions.tab'))
    )
    assert records
    assert len(records) == 5
    assert records[0]['variantType'] == 'sv'
    assert 'variant' in records[0]


def test_load_expression_variants() -> None:
    records = preprocess_expression_variants(
        read_tabbed_file(os.path.join(DATA_DIR, 'expression.tab'))
    )
    assert records
    assert len(records) == 4603
    assert records[0]['variantType'] == 'exp'
    assert 'variant' in records[0]


class TestCheckVariantLinks:
    def test_sm_missing_copy_empty_ok(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({'gene': 'KRAS'})],
            copy_variants=[],
            expression_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_sm_missing_exp_empty_ok(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({'gene': 'KRAS'})],
            copy_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],
            expression_variants=[],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_sm_missing_copy(self) -> None:
        with mock.patch.object(logger, 'verbose') as mock_debug:
            check_variant_links(
                small_mutations=[IprGeneVariant({'gene': 'KRAS'})],
                copy_variants=[IprGeneVariant({'gene': 'CDK', 'variant': ''})],
                expression_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],
                structural_variants=[],
            )
            assert mock_debug.called

    def test_sm_missing_exp(self) -> None:
        with mock.patch.object(logger, 'verbose') as mock_debug:
            check_variant_links(
                small_mutations=[IprGeneVariant({'gene': 'KRAS'})],
                copy_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],
                expression_variants=[IprGeneVariant({'gene': 'CDK', 'variant': ''})],
                structural_variants=[],
            )
            assert mock_debug.called

    def test_with_valid_inputs(self) -> None:
        genes = check_variant_links(
            small_mutations=[IprGeneVariant({'gene': 'KRAS'})],
            copy_variants=[
                IprGeneVariant({'gene': 'KRAS', 'variant': ''}),
                IprGeneVariant({'gene': 'CDK', 'variant': ''}),
            ],
            expression_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],
            structural_variants=[],
        )
        assert genes == {'KRAS'}

    def test_copy_missing_exp(self) -> None:
        with mock.patch.object(logger, 'verbose') as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[
                    IprGeneVariant({'gene': 'BRAF', 'variant': 'copy gain'}),
                    IprGeneVariant({'gene': 'KRAS', 'variant': ''}),
                ],
                expression_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],
                structural_variants=[],
            )
            assert mock_debug.called

    def test_exp_missing_copy(self) -> None:
        with mock.patch.object(logger, 'verbose') as mock_debug:
            check_variant_links(
                small_mutations=[],
                copy_variants=[IprGeneVariant({'gene': 'KRAS', 'variant': ''})],
                expression_variants=[
                    IprGeneVariant({'gene': 'BRAF', 'variant': 'increased expression'})
                ],
                structural_variants=[],
            )
            assert mock_debug.called


class TestCreateGraphkbSvNotation:
    def test_both_genes_and_exons(self) -> None:
        notation = create_graphkb_sv_notation(
            IprStructuralVariant({'gene1': 'A', 'gene2': 'B', 'exon1': 1, 'exon2': 2})
        )
        assert notation == '(A,B):fusion(e.1,e.2)'

    def test_one_exon_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprStructuralVariant({'gene1': 'A', 'gene2': 'B', 'exon1': '', 'exon2': 2})
        )
        assert notation == '(A,B):fusion(e.?,e.2)'

    def test_one_gene_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprStructuralVariant({'gene1': 'A', 'gene2': '', 'exon1': 1, 'exon2': 2})
        )
        assert notation == '(A,?):fusion(e.1,e.2)'

    def test_first_gene_missing(self) -> None:
        notation = create_graphkb_sv_notation(
            IprStructuralVariant({'gene1': '', 'gene2': 'B', 'exon1': 1, 'exon2': 2})
        )
        assert notation == '(B,?):fusion(e.2,e.1)'

    def test_no_genes_error(self) -> None:
        with pytest.raises(ValueError):
            create_graphkb_sv_notation(
                IprStructuralVariant({'gene1': '', 'gene2': '', 'exon1': 1, 'exon2': 2, 'key': 'x'})
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
