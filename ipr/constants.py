VARIANT_CLASSES = {'Variant', 'CategoryVariant', 'PositionalVariant', 'CatalogueVariant'}


APPROVED_EVIDENCE_LEVELS = {
    # sourceIds of levels by source name
    'oncokb': ['1', 'r1'],
    'profyle': ['t1'],
    'cancer genome interpreter': [
        'cpic guidelines',
        'european leukemianet guidelines',
        'fda guidelines',
        'nccn guidelines',
        'nccn/cap guidelines',
    ],
}

DEFAULT_URL = 'https://iprstaging-api.bcgsc.ca/api'
