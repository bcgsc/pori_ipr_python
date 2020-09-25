import os
import pytest
from typing import Dict
from unittest.mock import MagicMock, patch

from ipr.connection import IprConnection
from ipr.inputs import read_tabbed_file
from ipr.main import create_report

from .constants import EXCLUDE_INTEGRATION_TESTS


def get_test_file(name: str) -> str:
    return os.path.join(os.path.dirname(__file__), 'test_data', name)


@pytest.fixture(scope='module')
def report_upload_content() -> Dict:
    mock = MagicMock()
    with patch.object(IprConnection, 'upload_report', new=mock):
        create_report(
            patient_id='PATIENT001',
            project='TEST',
            expression_variant_rows=read_tabbed_file(get_test_file('expression.tab')),
            small_mutation_rows=read_tabbed_file(get_test_file('small_mutations.short.tab')),
            copy_variant_rows=read_tabbed_file(get_test_file('copy_variants.tab')),
            structural_variant_rows=read_tabbed_file(get_test_file('fusions.tab')),
            username=os.environ['IPR_USER'],
            password=os.environ['IPR_PASS'],
            log_level='info',
            ipr_url='http://fake.url.ca',
            kb_disease_match='colorectal cancer',
            optional_content={
                'blargh': 'some fake content',
                'comparators': [
                    {'analysisRole': 'expression (disease)'},
                    {'analysisRole': 'expression (primary site)'},
                    {'analysisRole': 'expression (biopsy site)'},
                ],
            },
        )

    assert mock.called

    report_content = mock.call_args[0][0]
    return report_content


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
class TestCreateReport:
    def test_main_sections_present(self, report_upload_content: Dict) -> None:
        sections = set(report_upload_content.keys())

        for section in [
            'structuralVariants',
            'expressionVariants',
            'copyVariants',
            'smallMutations',
            'kbMatches',
            'genes',
        ]:
            assert section in sections

    def test_kept_low_quality_fusion(self, report_upload_content: Dict) -> None:
        fusions = [(sv['gene1'], sv['gene2']) for sv in report_upload_content['structuralVariants']]
        assert ('SARM1', 'SUZ12') in fusions

    def test_pass_through_content_added(self, report_upload_content: Dict) -> None:
        # check the passthorough content was added
        assert 'blargh' in report_upload_content

    def test_found_fusion_partner_gene(self, report_upload_content: Dict) -> None:
        genes = report_upload_content['genes']
        assert any([g.get('knownFusionPartner', False) for g in genes])

    def test_found_oncogene(self, report_upload_content: Dict) -> None:
        genes = report_upload_content['genes']
        assert any([g.get('oncogene', False) for g in genes])

    def test_found_tumour_supressor(self, report_upload_content: Dict) -> None:
        genes = report_upload_content['genes']
        assert any([g.get('tumourSuppressor', False) for g in genes])

    def test_found_cancer_related_gene(self, report_upload_content: Dict) -> None:
        genes = report_upload_content['genes']
        assert any([g.get('cancerRelated', False) for g in genes])
