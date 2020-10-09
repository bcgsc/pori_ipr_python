import os
import pytest
from unittest import mock

from ipr.connection import IprConnection

IMAGE_DIR = os.path.join(os.path.dirname(__file__), '../docs/images')


class TestPostImages:
    def test_no_images_ok(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'successful'}], raise_for_status=lambda: None
            )
            return m

        with mock.patch('ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            result = conn.post_images('report_id', files={}, data={})
            assert result is None

    def test_images_load_ok(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'successful'}], raise_for_status=lambda: None
            )
            return m

        with mock.patch('ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            result = conn.post_images(
                'report_id',
                files={
                    'expression.correlation': os.path.join(IMAGE_DIR, 'expression_correlation.png'),
                    'mixcr.circos_trb_vj_gene_usage': os.path.join(
                        IMAGE_DIR, 'mixcr.circos_trb_vj_gene_usage.png'
                    ),
                },
                data={},
            )
            assert result is None

    def test_images_with_data_load_ok(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'successful'}], raise_for_status=lambda: None
            )
            return m

        with mock.patch('ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            result = conn.post_images(
                'report_id',
                files={
                    'expression.correlation': os.path.join(IMAGE_DIR, 'expression_correlation.png'),
                    'mixcr.circos_trb_vj_gene_usage': os.path.join(
                        IMAGE_DIR, 'mixcr.circos_trb_vj_gene_usage.png'
                    ),
                },
                data={'expression.correlation.title': 'this is a title'},
            )
            assert result is None

    def test_bad_file(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'successful'}], raise_for_status=lambda: None
            )
            return m

        with mock.patch('ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            with pytest.raises(FileNotFoundError):
                conn.post_images(
                    'report_id',
                    files={'expression.correlation': 'thing/that/does/not/exist.png'},
                )

    def test_failed_image_load(self):
        def request(*args, **kwargs):
            m = mock.MagicMock(
                json=lambda: [{'upload': 'anything else', 'key': 'thing'}],
                raise_for_status=lambda: None,
            )
            return m

        with mock.patch('ipr.connection.requests.request', request):
            conn = IprConnection('user', 'pass')
            with pytest.raises(ValueError):
                conn.post_images(
                    'report_id',
                    {
                        'expression.correlation': os.path.join(
                            IMAGE_DIR, 'expression_correlation.png'
                        )
                    },
                )
