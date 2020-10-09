from ipr.summary import (
    substitute_sentence_template,
    get_preferred_drug_representation,
    GRAPHKB_GUI,
)
from unittest.mock import MagicMock


class TestGetPreferredDrugRepresentation:
    def test_prefers_non_alias(self):
        api = MagicMock(
            query=MagicMock(
                side_effect=[
                    [],
                    [
                        {'sourceId': '1', 'alias': False, 'source': 'source', 'name': 'name'},
                        {'sourceId': '2', 'alias': True, 'source': 'source', 'name': 'name'},
                    ],
                ]
            )
        )
        rec = get_preferred_drug_representation(api, 'anything')
        assert rec['sourceId'] == '1'

    def test_prefers_non_deprecated(self):
        api = MagicMock(
            query=MagicMock(
                side_effect=[
                    [],
                    [
                        {'sourceId': '1', 'deprecated': True, 'source': 'source', 'name': 'name'},
                        {'sourceId': '2', 'deprecated': False, 'source': 'source', 'name': 'name'},
                    ],
                ]
            )
        )
        rec = get_preferred_drug_representation(api, 'anything')
        assert rec['sourceId'] == '2'

    def test_prefers_lower_sort_source(self):
        api = MagicMock(
            query=MagicMock(
                side_effect=[
                    [{'@rid': 'source2', 'sort': 0}, {'@rid': 'source1', 'sort': 1}],
                    [
                        {'sourceId': '1', 'deprecated': False, 'source': 'source1', 'name': 'name'},
                        {'sourceId': '2', 'deprecated': False, 'source': 'source2', 'name': 'name'},
                    ],
                ]
            )
        )
        rec = get_preferred_drug_representation(api, 'anything')
        assert rec['sourceId'] == '2'

    def test_prefers_newer_version(self):
        api = MagicMock(
            query=MagicMock(
                side_effect=[
                    [],
                    [
                        {
                            'sourceId': '2',
                            'deprecated': True,
                            'source': 'source',
                            'name': 'name',
                            'sourceIdVersion': '1',
                        },
                        {
                            'sourceId': '2',
                            'deprecated': True,
                            'source': 'source',
                            'name': 'name',
                            'sourceIdVersion': '2',
                        },
                    ],
                ]
            )
        )
        rec = get_preferred_drug_representation(api, 'anything')
        assert rec['sourceIdVersion'] == '1'


class TestSubstituteSentenceTemplate:
    def test_multiple_diseases_no_matches(self):
        template = "{conditions:variant} is associated with {relevance} to {subject} in {conditions:disease} ({evidence})"
        relevance = {'displayName': 'senitivity'}
        disease_matches = {'1'}
        diseases = [
            {'@class': 'Disease', '@rid': '2', 'displayName': f'disease 1'},
            {'@class': 'Disease', '@rid': '3', 'displayName': f'disease 2'},
        ]
        variants = [
            {
                '@class': 'CategoryVariant',
                'displayName': 'KRAS increased RNA expression',
                '@rid': '4',
            }
        ]
        subjects = [{'@class': 'Therapy', 'displayName': 'some drug', '@rid': '5'}]
        sentence = substitute_sentence_template(
            template,
            diseases + variants,
            subjects,
            relevance,
            [],
            ['6', '7'],
            disease_matches,
        )
        assert (
            sentence
            == f'KRAS increased RNA expression is associated with senitivity to some drug in other disease types (<a href="{GRAPHKB_GUI}/data/table?complex=eyJ0YXJnZXQiOiBbIjYiLCAiNyJdfQ%3D%3D&%40class=Statement" target="_blank" rel="noopener"></a>)'
        )

    def test_multiple_diseases_some_matches(self):
        template = "{conditions:variant} is associated with {relevance} to {subject} in {conditions:disease} ({evidence})"
        relevance = {'displayName': 'senitivity'}
        disease_matches = {'1'}
        diseases = [
            {'@class': 'Disease', '@rid': '2', 'displayName': f'disease 2'},
            {'@class': 'Disease', '@rid': '1', 'displayName': f'disease 1'},
            {'@class': 'Disease', '@rid': '3', 'displayName': f'disease 3'},
        ]
        variants = [
            {
                '@class': 'CategoryVariant',
                'displayName': 'KRAS increased RNA expression',
                '@rid': '4',
            }
        ]
        subjects = [{'@class': 'Therapy', 'displayName': 'some drug', '@rid': '5'}]
        sentence = substitute_sentence_template(
            template,
            diseases + variants,
            subjects,
            relevance,
            [],
            ['6', '7'],
            disease_matches,
        )
        assert (
            sentence
            == f'KRAS increased RNA expression is associated with senitivity to some drug in disease 1, and other disease types (<a href="{GRAPHKB_GUI}/data/table?complex=eyJ0YXJnZXQiOiBbIjYiLCAiNyJdfQ%3D%3D&%40class=Statement" target="_blank" rel="noopener"></a>)'
        )

    def test_multiple_diseases_only_matches(self):
        template = "{conditions:variant} is associated with {relevance} to {subject} in {conditions:disease} ({evidence})"
        relevance = {'displayName': 'senitivity'}
        disease_matches = {'1', '2', '3'}
        diseases = [
            {'@class': 'Disease', '@rid': '2', 'displayName': f'disease 2'},
            {'@class': 'Disease', '@rid': '1', 'displayName': f'disease 1'},
            {'@class': 'Disease', '@rid': '3', 'displayName': f'disease 3'},
        ]
        variants = [
            {
                '@class': 'CategoryVariant',
                'displayName': 'KRAS increased RNA expression',
                '@rid': '4',
            }
        ]
        subjects = [{'@class': 'Therapy', 'displayName': 'some drug', '@rid': '5'}]
        sentence = substitute_sentence_template(
            template,
            diseases + variants,
            subjects,
            relevance,
            [],
            ['6', '7'],
            disease_matches,
        )
        assert (
            sentence
            == f'KRAS increased RNA expression is associated with senitivity to some drug in disease 1, disease 2, and disease 3 (<a href="{GRAPHKB_GUI}/data/table?complex=eyJ0YXJnZXQiOiBbIjYiLCAiNyJdfQ%3D%3D&%40class=Statement" target="_blank" rel="noopener"></a>)'
        )
