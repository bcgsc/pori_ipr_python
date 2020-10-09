"""
upload variant and report information to IPR
"""
import pandas
from typing import Dict, List

from graphkb import GraphKBConnection

from .types import KbMatch, IprVariant
from .util import (
    get_terms_set,
    get_preferred_drug_representation,
    find_variant,
    create_variant_name_tuple,
)


def create_therapeutic_options(
    graphkb_conn: GraphKBConnection, kb_matches: List[KbMatch], variants: List[IprVariant]
) -> List[Dict]:
    """
    Generate therapeutic options summary from the list of kb-matches
    """
    options = []
    resistance_markers = get_terms_set(graphkb_conn, ['no sensitivity'])

    for match in kb_matches:
        row_type = 'therapeutic'
        if match['category'] != 'therapeutic' or match['relevance'] == 'eligibility':
            continue
        if match['kbRelevanceId'] in resistance_markers:
            row_type = 'chemoresistance'
        variant = find_variant(variants, match['variantType'], match['variant'])
        drug = get_preferred_drug_representation(graphkb_conn, match['kbContextId'])

        gene, variant_string = create_variant_name_tuple(variant)

        options.append(
            {
                'gene': gene,
                'type': row_type,
                'therapy': drug['displayName'],
                'therapyGraphkbId': drug['@rid'],
                'context': match['relevance'],
                'contextGraphkbId': match['kbRelevanceId'],
                'variantGraphkbId': match['kbVariantId'],
                'variant': variant_string,
                'evidenceLevel': match['evidenceLevel'],
                'notes': match['kbStatementId'],
            }
        )
    if not options:
        return options
    options_df = pandas.DataFrame.from_records(options)

    def delimited_list(inputs, delimiter=' / ') -> str:
        return delimiter.join(sorted(list({i for i in inputs if i})))

    options_df = options_df.groupby(['gene', 'type', 'therapy', 'variant']).agg(
        {
            'evidenceLevel': delimited_list,
            'context': delimited_list,
            'notes': lambda x: delimited_list(x, ' '),
        }
    )
    options_df = options_df.reset_index()
    options = options_df.to_dict('records')
    therapeutic_rank = 0
    chemoresistance_rank = 0
    for option in options:
        if option['type'] == 'therapeutic':
            option['rank'] = therapeutic_rank
            therapeutic_rank += 1
        else:
            option['rank'] = chemoresistance_rank
            chemoresistance_rank += 1
    return options
