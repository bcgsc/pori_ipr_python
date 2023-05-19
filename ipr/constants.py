DEFAULT_URL = 'https://iprstaging-api.bcgsc.ca/api'
GERMLINE_BASE_TERMS = ('pharmacogenomic', 'cancer predisposition')  # based on graphkb.constants
VARIANT_CLASSES = {'Variant', 'CategoryVariant', 'PositionalVariant', 'CatalogueVariant'}

# all possible values for review status are: ['pending', 'not required', 'passed', 'failed', 'initial']
FAILED_REVIEW_STATUS = 'failed'

TMB_HIGH = 10.0  # genomic mutations per mb - https://www.bcgsc.ca/jira/browse/GERO-296
TMB_HIGH_CATEGORY = 'high mutation burden'
