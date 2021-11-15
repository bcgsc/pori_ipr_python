import pytest

from ipr.util import create_variant_name, trim_empty_values


@pytest.mark.parametrize(
    'input,output_keys',
    [[{'key': 0}, ['key']], [{'key': None}, []], [{'key': ''}, []], [{'gene1': None}, ['gene1']]],
)
def test_trim_empty_values(input, output_keys):
    modified_object = trim_empty_values(input)
    assert sorted(modified_object.keys()) == sorted(output_keys)


@pytest.mark.parametrize(
    'variant,result',
    [
        [
            {'variantType': 'exp', 'gene': 'GENE', 'expressionState': 'increased expression'},
            'GENE (increased expression)',
        ],
        [
            {'variantType': 'cnv', 'gene': 'GENE', 'cnvState': 'amplification'},
            'GENE (amplification)',
        ],
        [
            {'variantType': 'other', 'variant': 'anything'},
            'anything',
        ],
    ],
)
def test_create_variant_name(variant, result):
    name = create_variant_name(variant)
    assert name == result
