import os

from genomic_report.inputs import load_small_mutations, load_copy_variants


DATA_DIR = os.path.join(os.path.dirname(__file__), 'test_data')


def test_load_small_mutations():
    records = load_small_mutations(os.path.join(DATA_DIR, 'small_mutations.tab'))
    assert records
    assert len(records) == 2614


def test_load_copy_variants():
    records = load_copy_variants(os.path.join(DATA_DIR, 'copy_variants.tab'))
    assert records
    assert len(records) == 4599
