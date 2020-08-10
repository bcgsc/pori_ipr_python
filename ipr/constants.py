from graphkb.util import IterableNamespace


BASE_THERAPEUTIC_TERMS = ['therapeutic efficacy', 'eligibility']
BASE_DIAGNOSTIC_TERM = 'diagnostic indicator'
BASE_PROGNOSTIC_TERM = 'prognostic indicator'
BASE_BIOLOGICAL_TERMS = ['functional effect', 'tumourigenesis', 'predisposing']

VARIANT_CLASSES = {'Variant', 'CategoryVariant', 'PositionalVariant', 'CatalogueVariant'}


REPORT_KB_SECTIONS = IterableNamespace(
    therapeutic='therapeutic',
    prognostic='prognostic',
    biological='biological',
    unknown='unknown',
    novel='novel',
    diagnostic='diagnostic',
)

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
