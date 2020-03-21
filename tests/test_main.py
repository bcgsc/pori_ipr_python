
import os
from unittest.mock import patch, MagicMock
from argparse import Namespace


from genomic_report.ipr import IprConnection
from genomic_report.main import main


def get_test_file(name):
    return os.path.join(
        os.path.dirname(__file__),
        'test_data',
        name
    )


def test_report_upload(tmpdir):
    mock = MagicMock()

    with patch.object(IprConnection, 'upload_report', new=mock):

        args = Namespace(
            expression_variants=get_test_file('expression.tab'),
            small_mutations=get_test_file('small_mutations.tab'),
            copy_variants=get_test_file('copy_variants.tab'),
            structural_variants=get_test_file('fusions.tab'),
            username=os.environ['USERNAME'],
            password=os.environ['PASSWORD'],
            output_json=os.path.join(tmpdir, 'output.json'),
            log_level='info',
            ipr_url='http://fake.url.ca'
        )

        main(args)

    assert mock.called

    report_content = mock.call_args[0][0]
    sections = set(report_content.keys())

    for section in ['structuralVariants', 'expressionVariants', 'copyVariants', 'smallMutations', 'kbMatches', 'genes']:
        assert section in sections
