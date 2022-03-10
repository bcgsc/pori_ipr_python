import pytest
from graphkb import statement as gkb_statement
from graphkb import vocab as gkb_vocab
from graphkb.types import Statement
from unittest.mock import Mock

from ipr.ipr import convert_statements_to_alterations

DISEASE_RIDS = ['#138:12', '#138:13']
APPROVED_EVIDENCE_RIDS = ['approved1', 'approved2']


@pytest.fixture
def graphkb_conn():
    class QueryMock:
        return_values = [
            # get approved evidence levels
            [{'@rid': v} for v in APPROVED_EVIDENCE_RIDS],
        ]
        index = -1

        def __call__(self, *args, **kwargs):
            self.index += 1
            ret_val = self.return_values[self.index] if self.index < len(self.return_values) else []
            return ret_val

    conn = Mock(query=QueryMock(), cache={})

    return conn


def base_graphkb_statement(disease_id: str = 'disease', relevance_rid: str = 'other') -> Statement:
    statement = Statement(
        {
            'conditions': [
                {'@class': 'Disease', '@rid': disease_id, 'displayName': 'disease_display_name'},
                {
                    '@class': 'CategoryVariant',
                    '@rid': 'variant_rid',
                    'displayName': 'KRAS increased expression',
                },
            ],
            'evidence': [],
            'subject': None,
            'source': None,
            'sourceId': None,
            'relevance': {'@rid': relevance_rid, 'displayName': 'relevance_display_name'},
            '@rid': 'statement_rid',
        }
    )
    return statement


@pytest.fixture(autouse=True)
def mock_get_term_tree(monkeypatch):
    def mock_func(*pos, **kwargs):
        return [{'@rid': d} for d in DISEASE_RIDS]

    monkeypatch.setattr(gkb_vocab, 'get_term_tree', mock_func)


@pytest.fixture(autouse=True)
def mock_categorize_relevance(monkeypatch):
    def mock_func(_, relevance_id):
        return relevance_id

    monkeypatch.setattr(gkb_statement, 'categorize_relevance', mock_func)


class TestConvertStatementsToAlterations:
    def test_disease_match(self, graphkb_conn, mock_get_term_tree) -> None:
        statement = base_graphkb_statement(DISEASE_RIDS[0])
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )

        assert len(result) == 1
        row = result[0]
        assert row['kbVariantId'] == 'variant_rid'
        assert row['kbStatementId'] == 'statement_rid'
        assert row['matchedCancer']
        assert row['kbVariant'] == 'KRAS increased expression'
        assert row['relevance'] == 'relevance_display_name'

    def test_no_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement('other')
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )

        assert len(result) == 1
        row = result[0]
        assert not row['matchedCancer']

    def test_multiple_disease_not_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement('disease')
        statement['conditions'].append(
            {'@class': 'Disease', '@rid': 'other', 'displayName': 'disease_display_name'}
        )
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )

        assert len(result) == 1
        row = result[0]
        assert not row['matchedCancer']

    def test_biological(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'biological'

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'biological'

    def test_prognostic_no_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'prognostic'

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 0

    def test_prognostic_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement(DISEASE_RIDS[0])
        statement['relevance']['@rid'] = 'prognostic'

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'prognostic'

    def test_diagnostic(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'diagnostic'

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'diagnostic'

    def test_unapproved_therapeutic(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'therapeutic'
        statement['evidenceLevel'] = [{'@rid': 'other', 'displayName': 'level'}]

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'therapeutic'

    def test_approved_therapeutic(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement['relevance']['@rid'] = 'therapeutic'
        statement['evidenceLevel'] = [{'@rid': APPROVED_EVIDENCE_RIDS[0], 'displayName': 'level'}]

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], 'disease', {'variant_rid'}
        )
        assert len(result) == 1
        row = result[0]
        assert row['category'] == 'therapeutic'
