import os
from unittest.mock import patch, MagicMock
from argparse import Namespace
import tempfile

import pytest

from genomic_report.ipr import IprConnection
from genomic_report.main import main


def get_test_file(name):
    return os.path.join(os.path.dirname(__file__), 'test_data', name)


@pytest.fixture(scope='module')
def report_upload_content():
    mock = MagicMock()
    with tempfile.TemporaryDirectory() as tmp_dir:
        # can't use builtin tmpdir b/c not session scoped
        with patch.object(IprConnection, 'upload_report', new=mock):
            args = Namespace(
                expression_variants=get_test_file('expression.tab'),
                small_mutations=get_test_file('small_mutations.tab'),
                copy_variants=get_test_file('copy_variants.tab'),
                structural_variants=get_test_file('fusions.tab'),
                username=os.environ['USERNAME'],
                password=os.environ['PASSWORD'],
                output_json=os.path.join(tmp_dir, 'output.json'),
                log_level='info',
                ipr_url='http://fake.url.ca',
            )

            main(args, {'blargh': 'some fake content'})

    assert mock.called

    report_content = mock.call_args[0][0]
    return report_content


def test_main_sections_present(report_upload_content):
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


def test_pass_through_content_added(report_upload_content):
    # check the passthorough content was added
    assert 'blargh' in report_upload_content


def test_found_fusion_partner_gene(report_upload_content):
    genes = report_upload_content['genes']
    assert any([g.get('knownFusionPartner', False) for g in genes])


def test_found_oncogene(report_upload_content):
    genes = report_upload_content['genes']
    assert any([g.get('oncogene', False) for g in genes])


def test_found_tumour_supressor(report_upload_content):
    genes = report_upload_content['genes']
    assert any([g.get('tumourSuppressor', False) for g in genes])


def test_found_cancer_related_gene(report_upload_content):
    genes = report_upload_content['genes']
    assert any([g.get('cancerRelated', False) for g in genes])
