import argparse
import datetime
import json
import logging
import os
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
from graphkb import GraphKBConnection
from typing import Dict, List, Sequence

from .annotate import (
    annotate_copy_variants,
    annotate_expression_variants,
    annotate_msi,
    annotate_positional_variants,
    get_gene_information,
)
from .connection import IprConnection
from .constants import DEFAULT_URL
from .inputs import (
    check_comparators,
    check_variant_links,
    preprocess_copy_variants,
    preprocess_expression_variants,
    preprocess_small_mutations,
    preprocess_structural_variants,
    validate_report_content,
)
from .ipr import (
    create_key_alterations,
    filter_structural_variants,
    germline_kb_matches,
    select_expression_plots,
)
from .summary import summarize
from .therapeutic_options import create_therapeutic_options
from .types import IprVariant, KbMatch
from .util import LOG_LEVELS, logger, trim_empty_values

CACHE_GENE_MINIMUM = 5000


def file_path(path: str) -> str:
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f'{repr(path)} is not a valid filename. does not exist')
    return path


def timestamp() -> str:
    return datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')


def command_interface() -> None:
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
    req = parser.add_argument_group('required arguments')
    (req if not os.environ.get('USER') else parser).add_argument(
        '--username',
        required=not os.environ.get('USER'),
        default=os.environ.get('USER'),
        help='username to use connecting to graphkb/ipr',
    )
    req.add_argument('--password', required=True, help='password to use connecting to graphkb/ipr')
    req.add_argument(
        '-c', '--content', required=True, type=file_path, help="Report Content as JSON"
    )
    parser.add_argument('--ipr_url', default=DEFAULT_URL)
    parser.add_argument('--graphkb_url', default=None)
    parser.add_argument('--log_level', default='info', choices=LOG_LEVELS.keys())
    parser.add_argument(
        '--therapeutics', default=False, help='Generate therapeutic options', action='store_true'
    )
    parser.add_argument(
        '--skip_comments',
        default=False,
        action='store_true',
        help='Turn off generating the analyst comments section of the report',
    )
    parser.add_argument(
        '-o', '--output_json_path', help='path to a JSON to output the report upload body'
    )
    parser.add_argument(
        '-w',
        '--always_write_output_json',
        action="store_true",
        help='Write to output_json_path on successful IPR uploads instead of just when the upload fails',
    )

    args = parser.parse_args()

    with open(args.content, 'r') as fh:
        content = json.load(fh)

    create_report(
        username=args.username,
        password=args.password,
        content=content,
        ipr_url=args.ipr_url,
        graphkb_url=args.graphkb_url,
        log_level=args.log_level,
        output_json_path=args.output_json_path,
        always_write_output_json=args.always_write_output_json,
        generate_therapeutics=args.therapeutics,
        generate_comments=not args.skip_comments,
    )


def clean_unsupported_content(upload_content: Dict) -> Dict:
    """
    Remove unsupported content. This content is either added to facilitate creation
    or to support upcoming and soon to be supported content that we would like
    to implement but is not yet supported by the upload
    """
    drop_columns = ['variant', 'variantType', 'histogramImage']
    for variant_section in [
        'expressionVariants',
        'smallMutations',
        'copyVariants',
        'structuralVariants',
    ]:
        for variant in upload_content[variant_section]:
            for col in drop_columns:
                if col in variant:
                    del variant[col]

    for row in upload_content['kbMatches']:
        del row['kbContextId']
        del row['kbRelevanceId']
    return upload_content


def create_report(
    username: str,
    password: str,
    content: Dict,
    ipr_url: str = DEFAULT_URL,
    log_level: str = 'info',
    output_json_path: str = '',
    always_write_output_json: bool = False,
    ipr_upload: bool = True,
    interactive: bool = False,
    graphkb_url: str = '',
    generate_therapeutics: bool = False,
    generate_comments: bool = True,
    match_germline: bool = False,
    custom_kb_match_filter=None,
) -> Dict:
    """Run the matching and create the report JSON for upload to IPR.

    Args:
        username: the username for connecting to GraphKB and IPR
        password: the password for connecting to GraphKB and IPR
        ipr_url: base URL to use in connecting to IPR
        log_level: the logging level
        content: report content
        output_json_path: path to a JSON file to output the report upload body.
        always_write_output_json: with successful IPR upload
        ipr_upload: upload report to ipr
        interactive: progressbars for interactive users
        cache_gene_minimum: minimum number of genes required for gene name caching optimization
        generate_therapeutics: create therapeutic options for upload with the report
        generate_comments: create the analyst comments section for upload with the report
        match_germline: match only germline statements to germline events and non-germline statements to non-germline events.
        custom_kb_match_filter: function(List[kbMatch]) -> List[kbMatch]

    Returns:
        ipr_conn.upload_report return dictionary
    """
    # set the default logging configuration
    logging.basicConfig(
        level=LOG_LEVELS[log_level],
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        datefmt='%m-%d-%y %H:%M:%S',
    )
    # validate the JSON content follows the specification
    validate_report_content(content)
    kb_disease_match = content['kbDiseaseMatch']

    # validate the input variants
    small_mutations = preprocess_small_mutations(content.get('smallMutations', []))
    structural_variants = preprocess_structural_variants(content.get('structuralVariants', []))
    copy_variants = preprocess_copy_variants(content.get('copyVariants', []))
    expression_variants = preprocess_expression_variants(content.get('expressionVariants', []))
    if expression_variants:
        check_comparators(content, expression_variants)

    genes_with_variants = check_variant_links(
        small_mutations, expression_variants, copy_variants, structural_variants
    )

    # Setup connections
    ipr_conn = IprConnection(username, password, ipr_url)
    if graphkb_url:
        logger.info(f'connecting to graphkb: {graphkb_url}')
        graphkb_conn = GraphKBConnection(graphkb_url)
    else:
        graphkb_conn = GraphKBConnection()
    graphkb_conn.login(username, password)

    gkb_matches: List[KbMatch] = []

    # Signature category variants
    if 'tmburMutationBurden' in content.keys():
        logger.warning(
            'GERO-296 - not yet implemented - high tumour mutation burden category matching.'
        )

    msi = content.get('msi', [])
    msi_matches = []
    msi_variant: IprVariant = {}
    if msi:
        # only one msi variant per library
        if isinstance(msi, list):
            msi_cat = msi[0].get('kbCategory')
        elif isinstance(msi, str):
            msi_cat = msi
        else:
            msi_cat = msi.get('kbCategory')
            msi_variant = msi.copy()
        logger.info(f'Matching GKB msi {msi_cat}')
        msi_matches = annotate_msi(graphkb_conn, msi_cat, kb_disease_match)
        if msi_matches:
            msi_variant['kbCategory'] = msi_cat  # type: ignore
            msi_variant['variant'] = msi_cat
            msi_variant['key'] = msi_cat
            msi_variant['variantType'] = 'msi'
            logger.info(f"GERO-295 '{msi_cat}' matches {len(msi_matches)} msi statements.")
            gkb_matches.extend(msi_matches)
            logger.debug(f"\tgkb_matches: {len(gkb_matches)}")

    logger.info(f'annotating {len(small_mutations)} small mutations')
    gkb_matches.extend(
        annotate_positional_variants(
            graphkb_conn, small_mutations, kb_disease_match, show_progress=interactive
        )
    )
    logger.debug(f"\tgkb_matches: {len(gkb_matches)}")

    logger.info(f'annotating {len(structural_variants)} structural variants')
    gkb_matches.extend(
        annotate_positional_variants(
            graphkb_conn, structural_variants, kb_disease_match, show_progress=interactive
        )
    )
    logger.debug(f"\tgkb_matches: {len(gkb_matches)}")

    logger.info(f'annotating {len(copy_variants)} copy variants')
    gkb_matches.extend(
        annotate_copy_variants(
            graphkb_conn, copy_variants, kb_disease_match, show_progress=interactive
        )
    )
    logger.debug(f"\tgkb_matches: {len(gkb_matches)}")

    logger.info(f'annotating {len(expression_variants)} expression variants')
    gkb_matches.extend(
        annotate_expression_variants(
            graphkb_conn, expression_variants, kb_disease_match, show_progress=interactive
        )
    )
    logger.debug(f"\tgkb_matches: {len(gkb_matches)}")

    all_variants: Sequence[IprVariant]
    all_variants = expression_variants + copy_variants + structural_variants + small_mutations  # type: ignore
    if msi_matches:
        all_variants.append(msi_variant)  # type: ignore

    if match_germline:  # verify germline kb statements matched germline observed variants
        gkb_matches = germline_kb_matches(gkb_matches, all_variants)
        if gkb_matches:
            logger.info(f"Removing {len(gkb_matches)} germline events without medical matches.")

    if custom_kb_match_filter:
        logger.info(f'custom_kb_match_filter on {len(gkb_matches)} variants')
        gkb_matches = custom_kb_match_filter(gkb_matches)
        logger.info(f'\t custom_kb_match_filter left {len(gkb_matches)} variants')

    key_alterations, variant_counts = create_key_alterations(gkb_matches, all_variants)

    logger.info('fetching gene annotations')
    gene_information = get_gene_information(graphkb_conn, sorted(genes_with_variants))

    if generate_therapeutics:
        logger.info('generating therapeutic options')
        targets = create_therapeutic_options(graphkb_conn, gkb_matches, all_variants)
    else:
        targets = []

    logger.info('generating analyst comments')
    if generate_comments:
        comments = {
            'comments': summarize(
                graphkb_conn, gkb_matches, disease_name=kb_disease_match, variants=all_variants
            )
        }
    else:
        comments = {'comments': ''}

    # thread safe deep-copy the original content
    output = json.loads(json.dumps(content))
    output.update(
        {
            'kbMatches': [trim_empty_values(a) for a in gkb_matches],
            'copyVariants': [
                trim_empty_values(c) for c in copy_variants if c['gene'] in genes_with_variants
            ],
            'smallMutations': [trim_empty_values(s) for s in small_mutations],
            'expressionVariants': [
                trim_empty_values(e)
                for e in expression_variants
                if e['gene'] in genes_with_variants
            ],
            'kbDiseaseMatch': kb_disease_match,
            'kbUrl': graphkb_conn.url,
            'kbVersion': timestamp(),
            'structuralVariants': [
                trim_empty_values(s)
                for s in filter_structural_variants(
                    structural_variants, gkb_matches, gene_information
                )
            ],
            'genes': gene_information,
            'genomicAlterationsIdentified': key_alterations,
            'variantCounts': variant_counts,
            'analystComments': comments,
            'therapeuticTarget': targets,
        }
    )
    output.setdefault('images', []).extend(select_expression_plots(gkb_matches, all_variants))

    output = clean_unsupported_content(output)
    ipr_result = None
    upload_error = None

    if ipr_upload:
        try:
            logger.info(f'Uploading to IPR {ipr_conn.url}')
            ipr_result = ipr_conn.upload_report(output)
            logger.info(ipr_result)
            output.update(ipr_result)
        except Exception as err:
            upload_error = err
            logger.error(f"ipr_conn.upload_report failed: {err}", exc_info=True)
    if output_json_path:
        if always_write_output_json or not ipr_result:
            logger.info(f'Writing IPR upload json to: {output_json_path}')
            with open(output_json_path, 'w') as fh:
                fh.write(json.dumps(output))
    logger.info(f'made {graphkb_conn.request_count} requests to graphkb')
    logger.info(f'average load {int(graphkb_conn.load or 0)} req/s')
    if upload_error:
        raise upload_error
    return output
