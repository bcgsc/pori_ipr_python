import pytest

from ipr.util import create_variant_name_tuple, trim_empty_values


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
            'increased expression',
        ],
        [{'variantType': 'cnv', 'gene': 'GENE', 'cnvState': 'amplification'}, 'amplification'],
        [{'variantType': 'other', 'gene2': 'GENE', 'variant': 'GENE:anything'}, 'anything'],
    ],
)
def test_create_variant_name_tuple(variant, result):
    gene, name = create_variant_name_tuple(variant)
    assert name == result
    assert gene == 'GENE'
