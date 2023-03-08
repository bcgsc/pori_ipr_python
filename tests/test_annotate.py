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

from .constants import EXCLUDE_INTEGRATION_TESTS
from .util import QueryMock

# TP53 examples from https://www.bcgsc.ca/jira/browse/SDEV-3122
# Mutations are actually identical but on alternate transcripts.

TP53_MUT_DICT = {
    'pref': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'SDEV-3122_preferred',
            'gene': 'TP53',
            'hgvsGenomic': 'chr17:g.7674252C>T',
            'hgvsCds': 'ENST00000269305:c.711G>A',
            'hgvsProtein': 'TP53:p.M237I',
        }
    ),
    'intersect': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'SDEV-3122_alt',
            'gene': 'TP53',
            'hgvsGenomic': 'chr17:g.7674252C>T',
            'hgvsCds': 'ENST00000610292:c.594G>A',
            'hgvsProtein': 'TP53:p.M198I',
        }
    ),
    'prot_only': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'prot_only',
            'gene': 'TP53',
            'hgvsProtein': 'TP53:p.M237I',
        }
    ),
    'cds_only': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'cds_only',
            'gene': 'TP53',
            'hgvsCds': 'ENST00000269305:c.711G>A',
        }
    ),
    'genome_only': IprSmallMutationVariant(  # type: ignore
        {
            'key': 'genome_only',
            'gene': 'TP53',
            'hgvsGenomic': 'chr17:g.7674252C>T',
        }
    ),
}


@pytest.fixture(scope='module')
def graphkb_conn():
    username = os.environ['IPR_USER']
    password = os.environ['IPR_PASS']
    graphkb_conn = GraphKBConnection()
    graphkb_conn.login(username, password)
    return graphkb_conn


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="SDEV-3381 - github workflow failures.")
def test_annotate_nonsense_vs_missense(graphkb_conn):
    """Verify missense (point mutation) is not mistaken for a nonsense (stop codon) mutation."""
    disease = 'cancer'
    for key in ('prot_only', 'cds_only', 'genome_only', 'pref'):
        matched = annotate_positional_variants(graphkb_conn, [TP53_MUT_DICT[key]], disease)
        # nonsense - stop codon - should not match.  This is missense not nonsense (#164:933).
        nonsense = [a for a in matched if a['kbVariant'] == 'TP53 nonsense']
        assert not nonsense, f"nonsense matched to {key}: {TP53_MUT_DICT[key]}"
        assert matched, f"should have matched in {key}: {TP53_MUT_DICT[key]}"


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="SDEV-3381 - github workflow failures.")
def test_annotate_nonsense_vs_missense_protein(graphkb_conn):
    """Verify missense (point mutation) is not mistaken for a nonsense (stop codon) mutation."""
    disease = 'cancer'
    for key in ('prot_only', 'pref'):
        matched = annotate_positional_variants(graphkb_conn, [TP53_MUT_DICT[key]], disease)
        # nonsense - stop codon - should not match.  This is missense not nonsense (#164:933).
        nonsense = [a for a in matched if 'nonsense' in a['kbVariant']]
        assert not nonsense, f"nonsense matched to {key}: {TP53_MUT_DICT[key]}"
        assert matched, f"should have matched in {key}: {TP53_MUT_DICT[key]}"


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="SDEV-3381 - github workflow failures.")
def test_annotate_structural_variants_tp53(graphkb_conn):
    """Verify alternate TP53 variants match."""
    disease = 'cancer'
    ref_key = 'prot_only'
    pref = annotate_positional_variants(graphkb_conn, [TP53_MUT_DICT[ref_key]], disease)
    known_issues = set(['TP53:p.M237X'])  # SDEV-3122 -
    # GERO-299 - nonsense - stop codon - should not match.  This is missense not nonsense (#164:933).
    nonsense = [a for a in pref if a['kbVariant'] == 'TP53 nonsense']
    assert not nonsense
    pref_vars = set([m['kbVariant'] for m in pref])
    assert pref_vars, f"No matches to {TP53_MUT_DICT[pref]}"
    print(pref_vars)
    for key, alt_rep in TP53_MUT_DICT.items():
        if key == ref_key:
            continue
        alt = annotate_positional_variants(graphkb_conn, [alt_rep], disease)
        alt_vars = set([m['kbVariant'] for m in alt])
        diff = pref_vars.symmetric_difference(alt_vars)
        missing = pref_vars.difference(alt_vars)

        known_issues = set([])
        if 'hgvsCds' in alt_rep:
            known_issues.add('TP53 nonsense')  # GERO-299
        if 'p.M237' not in alt_rep:
            known_issues.add('TP53:p.M237X')  # SDEV-3122 - not matching imprecise mutations
        if key == 'genome_only':
            known_issues.add('TP53 mutation')

        # strangely genome_only matched to more precise type 'TP53 deleterious mutation' but not 'TP53 mutation'
        missing = pref_vars.difference(alt_vars).difference(known_issues)
        print(alt_vars)
        assert not missing, f"{key} missing{missing}: {diff}"


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="SDEV-3381 - github workflow failures.")
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
