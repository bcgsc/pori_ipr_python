import os
import pytest
from graphkb import GraphKBConnection
from graphkb import genes as gkb_genes
from graphkb import match as gkb_match
from graphkb import vocab as gkb_vocab
from unittest.mock import Mock

from ipr.annotate import (
    annotate_positional_variants,
    get_gene_information,
    get_therapeutic_associated_genes,
)
from ipr.constants import FAILED_REVIEW_STATUS
from ipr.types import IprSmallMutationVariant

from .util import QueryMock

# TP53 examples from https://www.bcgsc.ca/jira/browse/SDEV-3122
# Mutations are actually identical but on alternate transcripts.
TP53_ALT = IprSmallMutationVariant(  # type: ignore
    {
        'key': '1',
        'gene': 'TP53',
        'hgvsGenomic': 'chr17:g.7674252C>T',
        'hgvsCds': 'ENST00000610292:c.594G>A',
        'hgvsProtein': 'TP53:p.M198I',
        'variantType': 'mut',
    }
)

TP53_PREF = IprSmallMutationVariant(  # type: ignore
    {
        'key': '2',
        'gene': 'TP53',
        'hgvsGenomic': 'chr17:g.7674252C>T',
        'hgvsCds': 'ENST00000269305:c.711G>A',
        'hgvsProtein': 'TP53:p.M237I',
        'variantType': 'mut',
    }
)


@pytest.fixture(scope='module')
def graphkb_conn():
    username = os.environ['IPR_USER']
    password = os.environ['IPR_PASS']
    graphkb_conn = GraphKBConnection()
    graphkb_conn.login(username, password)
    return graphkb_conn


def test_annotate_structural_variants(graphkb_conn):
    """Verify alternate TP53 variants match."""
    disease = 'cancer'
    tp53_pref = annotate_positional_variants(graphkb_conn, [TP53_PREF], disease)
    tp53_alt = annotate_positional_variants(graphkb_conn, [TP53_ALT], disease)
    # Verify the missed cross matches are only TP53:p.M237X to TP53:p.M198I
    alt_ids = set([m['kbVariantId'] for m in tp53_alt])
    pref_ids = set([m['kbVariantId'] for m in tp53_pref])
    assert pref_ids.issuperset(alt_ids)
    missed_matches = [m['kbVariant'] for m in tp53_pref if m['kbVariantId'] not in alt_ids]
    if missed_matches:
        assert (sorted(set(missed_matches))) == ['TP53:p.M237X']


def test_get_therapeutic_associated_genes(graphkb_conn):
    gene_list = get_therapeutic_associated_genes(graphkb_conn=graphkb_conn)
    assert gene_list, 'No get_therapeutic_associated_genes found'
    assert (
        len(gene_list) > 500
    ), f'Expected over 500 get_therapeutic_associated_genes but found {len(gene_list)}'


@pytest.mark.parametrize(
    'gene,flags',
    [
        ['fusionPartnerGene', ['knownFusionPartner', 'cancerRelated']],
        ['smallMutationGene', ['knownSmallMutation', 'cancerRelated']],
        ['cancerRelatedGene', ['cancerRelated']],
        # ['therapyAssociatedGene1', ['therapeuticAssociated']],
        # ['therapyAssociatedGene2', ['therapeuticAssociated']],
        # ['therapyAssociatedGene3', ['therapeuticAssociated']],
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

    info = get_gene_information(graphkb_conn, [gene])

    assert info, f"get_gene_information failed for {gene}"
    gene_info = [i for i in info if i['name'] == gene]
    assert len(gene_info) == 1
    gene_info = gene_info[0]
    for flag in flags:
        assert flag in gene_info
        assert gene_info[flag]

    for attr in gene_info:
        if attr not in {'name'} | set(flags):
            assert not gene_info[attr], f'expected {attr} to be False'
