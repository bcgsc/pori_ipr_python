import os

import pytest
from graphkb import GraphKBConnection

from genomic_report.annotate import get_gene_information


@pytest.fixture(scope='class')
def genes():
    graphkb_conn = GraphKBConnection()
    graphkb_conn.login(os.environ['USERNAME'], os.environ['PASSWORD'])

    return get_gene_information(graphkb_conn, ['kras', 'cdkn2a', 'blargh-monkeys', 'ewsr1'])


class TestGetGeneInformation:
    def test_fetches_tumour_suppressors(self, genes):
        assert genes
        gene = [g for g in genes if g['name'] == 'cdkn2a']
        assert gene
        assert 'tumourSuppressor' in gene[0]
        assert gene[0]['tumourSuppressor']

    def test_fetches_oncogenes(self, genes):
        assert genes
        gene = [g for g in genes if g['name'] == 'kras']
        assert gene
        assert 'oncogene' in gene[0]
        assert gene[0]['oncogene']

    def test_fetches_cancer_genes(self, genes):
        assert genes
        cancer_genes = [g for g in genes if g.get('cancerRelated', False)]
        assert cancer_genes

    def test_ignores_noninteresting_genes(self, genes):
        assert genes
        names = [g['name'] for g in genes]
        assert 'blargh-monkeys' not in names

    def test_fetches_fusion_partner_genes(self, genes):
        assert genes
        names = [g['name'] for g in genes if g.get('knownFusionPartner')]
        assert 'ewsr1' in names
