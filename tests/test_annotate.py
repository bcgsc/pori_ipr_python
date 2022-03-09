import pytest
from graphkb import genes as gkb_genes
from graphkb import match as gkb_match
from graphkb import vocab as gkb_vocab
from unittest.mock import Mock

from ipr.annotate import get_gene_information
from ipr.constants import FAILED_REVIEW_STATUS

from .util import QueryMock


@pytest.mark.parametrize(
    'gene,flags',
    [
        ['fusionPartnerGene', ['knownFusionPartner', 'cancerRelated']],
        ['smallMutationGene', ['knownSmallMutation', 'cancerRelated']],
        ['cancerRelatedGene', ['cancerRelated']],
        ['therapyAssociatedGene1', ['therapeuticAssociated']],
        ['therapyAssociatedGene2', ['therapeuticAssociated']],
        ['therapyAssociatedGene3', ['therapeuticAssociated']],
        ['oncoGene', ['oncogene']],
        ['tumourSuppressorGene', ['tumourSuppressor']],
    ],
)
def test_get_gene_information(gene, flags, monkeypatch):

    # mock the API connection class
    graphkb_conn = Mock(
        query=QueryMock(
            [
                # variants to return
                [
                    {'reference1': 'fusionGene1', 'reference2': 'fusionPartnerGene'},
                    {
                        'reference1': 'smallMutationGene',
                        '@class': 'PositionalVariant',
                        'reference2': None,
                    },
                    {
                        'reference1': 'cancerRelatedGene',
                        '@class': 'CategoryVariant',
                        'reference2': None,
                    },
                ],
                # mock the return statements
                [
                    {'reviewStatus': FAILED_REVIEW_STATUS},
                    {
                        'reviewStatus': '',
                        'conditions': [
                            {'@class': 'Feature', '@rid': 'therapyAssociatedGene1'},
                            {
                                '@class': 'Variant',
                                'reference1': {
                                    '@rid': 'therapyAssociatedGene2',
                                    '@class': 'Feature',
                                },
                                'reference2': None,
                            },
                            {
                                '@class': 'PositionalVariant',
                                'reference1': {
                                    '@rid': 'therapyAssociatedGene2',
                                    '@class': 'Feature',
                                },
                                'reference2': None,
                            },
                            {
                                '@class': 'CategoryVariant',
                                'reference1': {'@rid': 'someGene', '@class': 'Feature'},
                                'reference2': {
                                    '@rid': 'therapyAssociatedGene3',
                                    '@class': 'Feature',
                                },
                            },
                        ],
                    },
                ],
            ]
        ),
        cache={},
    )

    monkeypatch.setattr(gkb_vocab, 'get_terms_set', lambda conn, term: ['anything'])
    monkeypatch.setattr(gkb_genes, 'get_oncokb_oncogenes', lambda conn: [{'@rid': 'oncoGene'}])
    monkeypatch.setattr(
        gkb_genes, 'get_oncokb_tumour_supressors', lambda conn: [{'@rid': 'tumourSuppressorGene'}]
    )
    monkeypatch.setattr(gkb_match, 'get_equivalent_features', lambda conn, term: [{'@rid': term}])

    info = get_gene_information(
        graphkb_conn,
        [gene],
    )

    assert info
    gene_info = [i for i in info if i['name'] == gene]
    assert len(gene_info) == 1
    gene_info = gene_info[0]
    for flag in flags:
        assert flag in gene_info
        assert gene_info[flag]

    for attr in gene_info:
        if attr not in {'name'} | set(flags):
            assert not gene_info[attr], f'expected {attr} to be False'
