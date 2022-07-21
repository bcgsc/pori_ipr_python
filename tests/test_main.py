import json
import os
import pandas as pd
import pytest
import sys
from typing import Dict
from unittest.mock import MagicMock, patch

from ipr.connection import IprConnection
from ipr.main import command_interface

from .constants import EXCLUDE_INTEGRATION_TESTS


def get_test_file(name: str) -> str:
    return os.path.join(os.path.dirname(__file__), 'test_data', name)


@pytest.fixture(scope='module')
def report_upload_content(tmp_path_factory) -> Dict:
    mock = MagicMock()
    json_file = tmp_path_factory.mktemp('inputs') / 'content.json'
    json_file.write_text(
        json.dumps(
            {
                'blargh': 'some fake content',
                'comparators': [
                    {'analysisRole': 'expression (disease)', 'name': '1'},
                    {'analysisRole': 'expression (primary site)', 'name': '2'},
                    {'analysisRole': 'expression (biopsy site)', 'name': '3'},
                    {'analysisRole': 'expression (internal pancancer cohort)', 'name': '4'},
                ],
                'patientId': 'PATIENT001',
                'project': 'TEST',
                'expressionVariants': pd.read_csv(
                    get_test_file('expression.tab'), sep='\t'
                ).to_dict('records'),
                'smallMutations': pd.read_csv(
                    get_test_file('small_mutations.short.tab'), sep='\t'
                ).to_dict('records'),
                'copyVariants': pd.read_csv(get_test_file('copy_variants.tab'), sep='\t').to_dict(
                    'records'
                ),
                'structuralVariants': pd.read_csv(get_test_file('fusions.tab'), sep='\t').to_dict(
                    'records'
                ),
                'kbDiseaseMatch': 'colorectal cancer',
            }
        )
    )

    with patch.object(
        sys,
        'argv',
        [
            'ipr',
            '--username',
            os.environ.get('IPR_USER', os.environ['USER']),
            '--password',
            os.environ['IPR_PASS'],
            '--ipr_url',
            'http://fake.url.ca',
            '--content',
            str(json_file),
            '--therapeutics',
        ],
    ):
        with patch.object(IprConnection, 'upload_report', new=mock):
            command_interface()

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
