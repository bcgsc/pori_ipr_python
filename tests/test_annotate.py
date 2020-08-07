import os
from typing import Dict, List

import pytest
from graphkb import GraphKBConnection

from ipr.annotate import get_gene_information

from .constants import EXCLUDE_INTEGRATION_TESTS


@pytest.fixture(scope='class')
def genes() -> List[Dict]:
    graphkb_conn = GraphKBConnection()
    graphkb_conn.login(os.environ['IPR_USER'], os.environ['IPR_PASS'])

    return get_gene_information(graphkb_conn, ['kras', 'cdkn2a', 'blargh-monkeys', 'ewsr1'])


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
class TestGetGeneInformation:
    def test_fetches_tumour_suppressors(self, genes: List[Dict]) -> None:
        assert genes
        gene = [g for g in genes if g['name'] == 'cdkn2a']
        assert gene
        assert 'tumourSuppressor' in gene[0]
        assert gene[0]['tumourSuppressor']

    def test_fetches_oncogenes(self, genes: List[Dict]) -> None:
        assert genes
        gene = [g for g in genes if g['name'] == 'kras']
        assert gene
        assert 'oncogene' in gene[0]
        assert gene[0]['oncogene']

    def test_fetches_cancer_genes(self, genes: List[Dict]) -> None:
        assert genes
        cancer_genes = [g for g in genes if g.get('cancerRelated', False)]
        assert cancer_genes

    def test_ignores_noninteresting_genes(self, genes: List[Dict]) -> None:
        assert genes
        names = [g['name'] for g in genes]
        assert 'blargh-monkeys' not in names

    def test_fetches_fusion_partner_genes(self, genes: List[Dict]) -> None:
        assert genes
        names = [g['name'] for g in genes if g.get('knownFusionPartner')]
        assert 'ewsr1' in names

    def test_fetches_small_mutation_genes(self, genes: List[Dict]) -> None:
        assert genes
        names = [g['name'] for g in genes if g.get('knownSmallMutation')]
        assert 'kras' in names

    def test_fetches_therapeutic_genes(self, genes: List[Dict]) -> None:
        assert genes
        names = [g['name'] for g in genes if g.get('therapeuticAssociated')]
        assert 'cdkn2a' in names
