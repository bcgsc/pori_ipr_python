import os
import pandas as pd
import pytest
from typing import Dict
from unittest.mock import MagicMock, patch

from ipr.connection import IprConnection
from ipr.main import create_report

from .constants import EXCLUDE_INTEGRATION_TESTS


def get_test_file(name: str) -> str:
    return os.path.join(os.path.dirname(__file__), 'test_data', name)


@pytest.fixture(scope='module')
def probe_upload_content() -> Dict:
    mock = MagicMock()
    with patch.object(IprConnection, 'upload_report', new=mock):
        create_report(
            content={
                'patientId': 'PATIENT001',
                'project': 'TEST',
                'smallMutations': pd.read_csv(
                    get_test_file('small_mutations_probe.tab'),
                    sep='\t',
                    dtype={'chromosome': 'string'},
                ).to_dict('records'),
                'structuralVariants': pd.read_csv(get_test_file('fusions.tab'), sep='\t').to_dict(
                    'records'
                ),
                'blargh': 'some fake content',
                'kbDiseaseMatch': 'colorectal cancer',
            },
            username=os.environ['IPR_USER'],
            password=os.environ['IPR_PASS'],
            log_level='info',
            ipr_url='http://fake.url.ca',
        )

    assert mock.called

    report_content = mock.call_args[0][0]
    return report_content


@pytest.mark.skipif(EXCLUDE_INTEGRATION_TESTS, reason="excluding long running integration tests")
class TestCreateReport:
    def test_found_probe_small_mutations(self, probe_upload_content: Dict) -> None:
        assert probe_upload_content['smallMutations']

    def test_found_probe_small_mutations_match(self, probe_upload_content: Dict) -> None:
        # verify each probe had a KB match
        for sm_probe in probe_upload_content['smallMutations']:
            match_list = [
                kb_match
                for kb_match in probe_upload_content['kbMatches']
                if kb_match['variant'] == sm_probe["key"]
            ]
            assert (
                match_list
            ), f"probe match failure: {sm_probe['gene']} {sm_probe['proteinChange']} key: {sm_probe['proteinChange']}"
