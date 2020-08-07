import os
from typing import Dict
from unittest.mock import MagicMock, patch

import pytest

from ipr.connection import IprConnection
from ipr.main import create_report
from ipr.inputs import read_tabbed_file

from .constants import EXCLUDE_INTEGRATION_TESTS


def get_test_file(name: str) -> str:
    return os.path.join(os.path.dirname(__file__), 'test_data', name)


@pytest.fixture(scope='module')
def probe_upload_content() -> Dict:
    mock = MagicMock()
    with patch.object(IprConnection, 'upload_report', new=mock):
        create_report(
            patient_id='PATIENT001',
            project='TEST',
            small_mutation_rows=read_tabbed_file(get_test_file('small_mutations_probe.tab')),
            structural_variant_rows=read_tabbed_file(get_test_file('fusions.tab')),
            username=os.environ['IPR_USER'],
            password=os.environ['IPR_PASS'],
            log_level='info',
            ipr_url='http://fake.url.ca',
            kb_disease_match='colorectal cancer',
            optional_content={'blargh': 'some fake content'},
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
