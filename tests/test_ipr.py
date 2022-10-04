import pytest
from graphkb import statement as gkb_statement
from graphkb import vocab as gkb_vocab
from unittest.mock import Mock, patch

from ipr.ipr import convert_statements_to_alterations, germline_kb_matches
from ipr.types import GkbStatement

DISEASE_RIDS = ["#138:12", "#138:13"]
APPROVED_EVIDENCE_RIDS = ["approved1", "approved2"]
GERMLINE_VARIANTS = [
    {
        "key": "1",
        "germline": True,
        "hgvsCds": "SLC28A3:c.1381C>T",
        "hgvsGenomic": "chr9:g.84286011G>A",
        "hgvsProtein": "SLC28A3:p.L461L",
        "ncbiBuild": "GRCh38",
        "normalAltCount": 37,
        "normalDepth": 37,
        "normalRefCount": 0,
        "proteinChange": "p.L461L",
        "rnaAltCount": "",
        "rnaDepth": "",
        "rnaRefCount": "",
        "startPosition": 84286011,
        "transcript": "ENST00000376238",
        "tumourAltCount": "",
        "tumourDepth": "",
        "tumourRefCount": "",
        "variant": "SLC28A3:p.L461L",
        "variantType": "mut",
        "zygosity": "",
    },
    {
        "key": "2",
        "germline": True,
        "hgvsCds": "BRCA1:c.4837A>",
        "hgvsGenomic": "chr17:g.43071077T>C",
        "hgvsProtein": "BRCA1:p.S1613G",
        "normalAltCount": 33,
        "normalDepth": 33,
        "normalRefCount": 0,
        "tumourAltCount": 37,
        "tumourDepth": 37,
        "tumourRefCount": 0,
    },
]

SOMATIC_VARIANTS = [
    {
        "key": "1",
        "gene": "SLC28A3",
        "germline": False,
        "hgvsCds": "SLC28A3:c.1381C>T",
        "hgvsGenomic": "chr9:g.84286011G>A",
        "hgvsProtein": "SLC28A3:p.L461L",
        "ncbiBuild": "GRCh38",
        "normalAltCount": 0,
        "normalDepth": 37,
        "normalRefCount": 37,
        "tumourAltCount": 37,
        "tumourDepth": 37,
        "tumourRefCount": 0,
        "variant": "SLC28A3:p.L461L",
        "variantType": "mut",
        "zygosity": "",
    },
    {
        "key": "2",
        "germline": False,
        "hgvsCds": "BRCA1:c.4837A>",
        "hgvsGenomic": "chr17:g.43071077T>C",
        "hgvsProtein": "BRCA1:p.S1613G",
        "normalAltCount": 1,
        "normalDepth": 33,
        "normalRefCount": 32,
        "tumourAltCount": 37,
        "tumourDepth": 37,
        "tumourRefCount": 0,
    },
]

GERMLINE_KB_MATCHES = [
    {
        "variant": "1",
        "approvedTherapy": False,
        "category": "pharmacogenomic",
        "context": "anthracyclines",
        "kbContextId": "#122:20944",
        "kbRelevanceId": "#147:38",
        "kbStatementId": "#154:13387",
        "kbVariant": "SLC28A3:c.1381C>T",
        "kbVariantId": "#159:5426",
        "matchedCancer": False,
        "reference": "PMID: 27197003",
        "relevance": "decreased toxicity",
        "reviewStatus": "initial",
    },
    {
        "variant": "2",
        "approvedTherapy": True,
        "category": "cancer predisposition",
        "kbContextId": "#135:8764",
        "kbRelevanceId": "#147:32",
        "kbStatementId": "#155:13511",
        "kbVariant": "BRCA1 mutation",
        "kbVariantId": "#161:938",
        "matchedCancer": False,
        "reference": "MOAlmanac FDA-56",
        "relevance": "therapy",
        "reviewStatus": None,
    },
]

SOMATIC_KB_MATCHES = [
    {
        "variant": "1",
        "approvedTherapy": False,
        "category": "prognostic",
        "kbContextId": "somatic_test",
        "kbRelevanceId": "#147:38",
        "kbStatementId": "#154:13387",
        "kbVariant": "SLC28A3:c.1381C>T",
        "kbVariantId": "#159:5426",
        "relevance": "prognostic",
        "reviewStatus": "initial",
    },
    {
        "variant": "2",
        "approvedTherapy": True,
        "category": "therapy",
        "kbContextId": "#135:8764",
        "kbRelevanceId": "#147:32",
        "kbStatementId": "#155:13511",
        "kbVariant": "BRCA1 mutation",
        "kbVariantId": "#161:938",
        "matchedCancer": False,
        "reference": "MOAlmanac FDA-56",
        "relevance": "therapy",
        "reviewStatus": None,
    },
]


@pytest.fixture
def graphkb_conn():
    class QueryMock:
        return_values = [
            # get approved evidence levels
            [{"@rid": v} for v in APPROVED_EVIDENCE_RIDS],
        ]
        index = -1

        def __call__(self, *args, **kwargs):
            self.index += 1
            ret_val = self.return_values[self.index] if self.index < len(self.return_values) else []
            return ret_val

    def mock_get_source(source):
        return {"@rid": 0}

    conn = Mock(query=QueryMock(), cache={}, get_source=mock_get_source)

    return conn


def base_graphkb_statement(
    disease_id: str = "disease", relevance_rid: str = "other"
) -> GkbStatement:
    statement = GkbStatement(  # type: ignore
        {
            "conditions": [
                {
                    "@class": "Disease",
                    "@rid": disease_id,
                    "displayName": "disease_display_name",
                },
                {
                    "@class": "CategoryVariant",
                    "@rid": "variant_rid",
                    "displayName": "KRAS increased expression",
                },
            ],
            "evidence": [],
            "subject": None,
            "source": None,
            "sourceId": None,
            "relevance": {
                "@rid": relevance_rid,
                "displayName": "relevance_display_name",
            },
            "@rid": "statement_rid",
        }
    )
    return statement


@pytest.fixture(autouse=True)
def mock_get_term_tree(monkeypatch):
    def mock_func(*pos, **kwargs):
        return [{"@rid": d} for d in DISEASE_RIDS]

    monkeypatch.setattr(gkb_vocab, "get_term_tree", mock_func)


@pytest.fixture(autouse=True)
def mock_categorize_relevance(monkeypatch):
    def mock_func(_, relevance_id):
        return relevance_id

    monkeypatch.setattr(gkb_statement, "categorize_relevance", mock_func)


class TestConvertStatementsToAlterations:
    def test_disease_match(self, graphkb_conn, mock_get_term_tree) -> None:
        statement = base_graphkb_statement(DISEASE_RIDS[0])
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )

        assert len(result) == 1
        row = result[0]
        assert row["kbVariantId"] == "variant_rid"
        assert row["kbStatementId"] == "statement_rid"
        assert row["matchedCancer"]
        assert row["kbVariant"] == "KRAS increased expression"
        assert row["relevance"] == "relevance_display_name"

    def test_no_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement("other")
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )

        assert len(result) == 1
        row = result[0]
        assert not row["matchedCancer"]

    def test_multiple_disease_not_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement("disease")
        statement["conditions"].append(
            {"@class": "Disease", "@rid": "other", "displayName": "disease_display_name"}  # type: ignore
        )
        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )

        assert len(result) == 1
        row = result[0]
        assert not row["matchedCancer"]

    def test_biological(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "biological"

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "biological"

    def test_prognostic_no_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "prognostic"

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )
        assert len(result) == 0

    def test_prognostic_disease_match(self, graphkb_conn) -> None:
        statement = base_graphkb_statement(DISEASE_RIDS[0])
        statement["relevance"]["@rid"] = "prognostic"

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "prognostic"

    def test_diagnostic(self, graphkb_conn) -> None:
        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "diagnostic"

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "diagnostic"

    @patch("ipr.ipr.get_evidencelevel_mapping")
    def test_unapproved_therapeutic(self, mock_get_evidencelevel_mapping, graphkb_conn) -> None:
        mock_get_evidencelevel_mapping.return_value = {"other": "test"}

        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "therapeutic"
        statement["evidenceLevel"] = [{"@rid": "other", "displayName": "level"}]  # type: ignore

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "therapeutic"

    @patch("ipr.ipr.get_evidencelevel_mapping")
    def test_approved_therapeutic(self, mock_get_evidencelevel_mapping, graphkb_conn) -> None:
        mock_get_evidencelevel_mapping.return_value = {APPROVED_EVIDENCE_RIDS[0]: "test"}

        statement = base_graphkb_statement()
        statement["relevance"]["@rid"] = "therapeutic"
        statement["evidenceLevel"] = [{"@rid": APPROVED_EVIDENCE_RIDS[0], "displayName": "level"}]  # type: ignore

        result = convert_statements_to_alterations(
            graphkb_conn, [statement], "disease", {"variant_rid"}
        )
        assert len(result) == 1
        row = result[0]
        assert row["category"] == "therapeutic"


class TestKbmatchFilters:
    def test_germline_kb_matches(self):
        assert len(germline_kb_matches(GERMLINE_KB_MATCHES, GERMLINE_VARIANTS)) == len(
            GERMLINE_KB_MATCHES
        ), "Germline variant not matched to germline KB statement."
        assert not germline_kb_matches(
            GERMLINE_KB_MATCHES, SOMATIC_VARIANTS
        ), "Somatic variant matched to KB germline statement."
        assert len(germline_kb_matches(SOMATIC_KB_MATCHES, SOMATIC_VARIANTS)) == len(
            SOMATIC_KB_MATCHES
        ), "Somatic variant not matched to somatic KB statement."
        assert not germline_kb_matches(
            SOMATIC_KB_MATCHES, GERMLINE_VARIANTS
        ), "Germline variant matched to KB somatic statement."
