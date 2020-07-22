import argparse
import datetime
import json
import logging
import os
from typing import Dict, List, Optional

from argparse_env import Action, ArgumentParser

from graphkb import GraphKBConnection
from graphkb.match import cache_missing_features

from . import ipr
from .annotate import annotate_category_variants, annotate_positional_variants, get_gene_information
from .inputs import (
    check_variant_links,
    load_copy_variants,
    load_expression_variants,
    load_small_mutations,
    load_structural_variants,
)
from .types import KbMatch
from .util import LOG_LEVELS, logger, trim_empty_values
from .summary import summarize

CACHE_GENE_MINIMUM = 5000


def file_path(path: str) -> str:
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f'{repr(path)} is not a valid filename. does not exist')
    return path


def timestamp() -> str:
    return datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')


def command_interface() -> None:
    parser = ArgumentParser()
    parser.add_argument(
        '--username',
        env=True,
        action=Action,
        required=True,
        help='username to use connecting to graphkb/ipr',
    )
    parser.add_argument(
        '--password',
        env=True,
        action=Action,
        required=True,
        sensitive=True,
        help='password to use connecting to graphkb/ipr',
    )
    parser.add_argument('-c', '--copy_variants', required=False, type=file_path)
    parser.add_argument('-m', '--small_mutations', required=False, type=file_path)
    parser.add_argument('-s', '--structural_variants', required=False, type=file_path)
    parser.add_argument('-e', '--expression_variants', required=False, type=file_path)
    parser.add_argument(
        '-d',
        '--kb_disease_match',
        required=True,
        help='Disease name to be used in matching to GraphKB',
    )
    parser.add_argument('--ipr_url', default=ipr.DEFAULT_URL)
    parser.add_argument('--graphkb_url', default=None)
    parser.add_argument('--log_level', default='info', choices=LOG_LEVELS.keys())
    parser.add_argument('--patient_id', required=True, help='The patient ID for this report')
    parser.add_argument('--project', default='TEST', help='The project to upload this report to')
    parser.add_argument(
        '-o', '--output_json_path', help='path to a JSON to output the report upload body',
    )
    parser.add_argument(
        '-w',
        '--always_write_output_json',
        action="store_true",
        help='Write to output_json_path on successful IPR uploads',
    )

    args = parser.parse_args()

    create_report(
        username=args.username,
        password=args.password,
        patient_id=args.patient_id,
        project=args.project,
        kb_disease_match=args.kb_disease_match,
        ipr_url=args.ipr_url,
        graphkb_url=args.graphkb_url,
        log_level=args.log_level,
        expression_variants_file=args.expression_variants,
        structural_variants_file=args.structural_variants,
        copy_variants_file=args.copy_variants,
        small_mutations_file=args.small_mutations,
        output_json_path=args.output_json_path,
        always_write_output_json=args.always_write_output_json,
    )


def clean_unsupported_content(upload_content: Dict) -> Dict:
    """
    Remove unsupported content. This content is either added to facilitate creation
    or to support upcoming and soon to be supported content that we would like
    to implement but is not yet supported by the upload
    """
    drop_columns = [
        'variant',
        'variantType',
        'histogramImage',
        'hgvs_protein',
        'hgvs_cds',
        'hgvs_genomic',
    ]
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
    return upload_content


def create_report(
    username: str,
    password: str,
    patient_id: str,
    kb_disease_match: str,
    project: str = 'TEST',
    ipr_url: str = ipr.DEFAULT_URL,
    log_level: str = 'info',
    expression_variants_file: str = None,
    structural_variants_file: str = None,
    copy_variants_file: str = None,
    small_mutations_file: str = None,
    optional_content: Optional[Dict] = None,
    output_json_path: str = None,
    always_write_output_json: bool = False,
    ipr_upload: bool = True,
    interactive: bool = False,
    cache_gene_minimum: int = CACHE_GENE_MINIMUM,
    graphkb_url: str = '',
) -> Optional[Dict]:
    """
    Run the matching and create the report JSON for upload to IPR

    Args:
        username: the username for connecting to GraphKB and IPR
        password: the password for connecting to GraphKB and IPR
        kb_disease_match: disease name to be used in matching to GraphKB
        ipr_url: base URL to use in connecting to IPR
        log_level: the logging level
        expression_variants_file: path to the expression variants input file
        structural_variants_file: path to the structural variants input file
        copy_variants_file: path to the copy number variants input file
        small_mutations_file: path to the small mutations input file
        optional_content: pass-through content to include in the JSON upload
        output_json_path: path to a JSON file to output the report upload body.
        always_write_output_json: with successful IPR upload
        ipr_upload: upload report to ipr
        interactive: progressbars for interactive users
        cache_gene_minimum: minimum number of genes required for gene name caching optimization
    Returns:
        ipr_conn.upload_report return dictionary
    """
    # set the default logging configuration
    logging.basicConfig(
        level=LOG_LEVELS[log_level],
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        datefmt='%m-%d-%y %H:%M:%S',
    )
    ipr_conn = ipr.IprConnection(username, password, ipr_url)
    if graphkb_url:
        print('connecting to graphkb', graphkb_url)
        graphkb_conn = GraphKBConnection(graphkb_url)
    else:
        graphkb_conn = GraphKBConnection()
    graphkb_conn.login(username, password)

    if copy_variants_file:
        logger.info(f'loading copy variants from: {copy_variants_file}')
        copy_variants = load_copy_variants(copy_variants_file)
        logger.info(f'loaded {len(copy_variants)}')
    else:
        logger.info("no copy variants given")
        copy_variants = []

    if small_mutations_file:
        logger.info(f'loading small mutations from: {small_mutations_file}')
        small_mutations = load_small_mutations(small_mutations_file)
        logger.info(f'loaded {len(small_mutations)} small mutations from: {small_mutations_file}')
    else:
        logger.info("no small mutations given")
        small_mutations = []

    if expression_variants_file:
        logger.info(f'loading expression variants from: {expression_variants_file}')
        expression_variants = load_expression_variants(expression_variants_file)
        logger.info(
            f'loaded {len(expression_variants)} expression variants from: {expression_variants_file}'
        )
    else:
        logger.info("no expression given")
        expression_variants = []

    if structural_variants_file:
        f'loading structural variants from: {structural_variants_file}'
        structural_variants = load_structural_variants(structural_variants_file)
        logger.info(
            f'loaded {len(structural_variants)} structural variants from: {structural_variants_file}'
        )
    else:
        logger.info("no structural variants given")
        structural_variants = []

    genes_with_variants = check_variant_links(
        small_mutations, expression_variants, copy_variants, structural_variants
    )

    # cache of the graphkb gene names speeds up calculation for large
    # numbers of genes, but has significant overhead and slows down
    # calculations on small numbers of genes.
    if len(genes_with_variants) > cache_gene_minimum:
        logger.info('caching genes to improve matching speed')
        cache_missing_features(graphkb_conn)

    # filter excess variants not required for extra gene information
    logger.info(f'annotating small mutations from: {small_mutations_file}')
    alterations: List[KbMatch] = annotate_positional_variants(
        graphkb_conn, small_mutations, kb_disease_match, show_progress=interactive
    )

    logger.info(f'annotating structural variants from: {structural_variants_file}')
    alterations.extend(
        annotate_positional_variants(
            graphkb_conn, structural_variants, kb_disease_match, show_progress=interactive
        )
    )

    logger.info(f'annotating copy variants from {copy_variants_file}')
    alterations.extend(
        annotate_category_variants(
            graphkb_conn, copy_variants, kb_disease_match, show_progress=interactive
        )
    )

    logger.info(f'annotating expression variants from: {expression_variants_file}')
    alterations.extend(
        annotate_category_variants(
            graphkb_conn,
            expression_variants,
            kb_disease_match,
            copy_variant=False,
            show_progress=interactive,
        )
    )

    logger.info('fetching gene annotations')
    gene_information = get_gene_information(graphkb_conn, genes_with_variants)

    output = optional_content or dict()

    key_alterations, variant_counts = ipr.create_key_alterations(
        alterations, expression_variants + copy_variants + structural_variants + small_mutations
    )

    output.update(
        {
            'patientId': patient_id,
            'project': project,
            'kbMatches': [trim_empty_values(a) for a in alterations],
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
                for s in ipr.filter_structural_variants(
                    structural_variants, alterations, gene_information
                )
            ],
            'genes': gene_information,
            'genomicAlterationsIdentified': key_alterations,
            'variantCounts': variant_counts,
        }
    )
    output.setdefault('images', []).extend(
        ipr.select_expression_plots(
            alterations, expression_variants + copy_variants + structural_variants + small_mutations
        )
    )
    for section in output:
        section_content_type = 'rows' if not isinstance(output[section], str) else 'characters'
        logger.info(f'section {section} has {len(output[section])} {section_content_type}')

    ipr_result = None
    report_id = None
    comments = summarize(
        graphkb_conn,
        alterations,
        disease_name=kb_disease_match,
        variants=expression_variants + copy_variants + structural_variants + small_mutations,
    )
    output = clean_unsupported_content(output)
    if ipr_upload:
        try:
            logger.info(f'Uploading to IPR {ipr_conn.url}')
            ipr_result = ipr_conn.upload_report(output)
            report_id = ipr_result['ident']
            logger.info(ipr_result)
            output.update(ipr_result)
        except Exception as err:
            logger.error(f"ipr_conn.upload_report failed: {err}", exc_info=True)
    if output_json_path:
        if always_write_output_json or not ipr_result:
            logger.info(f'Writing IPR upload json to: {output_json_path}')
            with open(output_json_path, 'w') as fh:
                fh.write(json.dumps(output))
    if report_id:
        try:
            ipr_conn.set_analyst_comments(report_id, {'comments': comments})
            logger.info(f'report {report_id} was annotated with generated comments')
        except Exception as err:
            logger.error(f"ipr_conn.set_analyst_comments failed: {err}", exc_info=True)
    logger.info(f'made {graphkb_conn.request_count} requests to graphkb')
    logger.info(f'average load {int(graphkb_conn.load or 0)} req/s')
    return output
