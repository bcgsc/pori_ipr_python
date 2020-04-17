import argparse
import logging
import os
import datetime
from typing import Dict, Optional

from argparse_env import ArgumentParser, Action
from graphkb import GraphKBConnection

from .inputs import (
    load_copy_variants,
    load_small_mutations,
    load_expression_variants,
    load_structural_variants,
    check_variant_links,
)
from .annotate import annotate_category_variants, annotate_positional_variants, get_gene_information
from .util import logger, LOG_LEVELS, trim_empty_values
from . import ipr


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
    parser.add_argument('--log_level', default='info', choices=LOG_LEVELS.keys())
    parser.add_argument('--patient_id', required=True, help='The patient ID for this report')
    parser.add_argument('--project', default='TEST', help='The project to upload this report to')
    parser.add_argument(
        '-w',
        '--write_to_json',
        help='path to a JSON to output the report upload body to on failure to upload',
    )

    args = parser.parse_args()

    create_report(
        username=args.username,
        password=args.password,
        patient_id=args.patient_id,
        project=args.project,
        kb_disease_match=args.kb_disease_match,
        ipr_url=args.ipr_url,
        log_level=args.log_level,
        expression_variants_file=args.expression_variants,
        structural_variants_file=args.structural_variants,
        copy_variants_file=args.copy_variants,
        small_mutations_file=args.small_mutations,
        write_to_json=args.write_to_json,
    )


def clean_unsupported_content(upload_content: Dict) -> Dict:
    for variant_section in [
        'expressionVariants',
        'smallMutations',
        'copyVariants',
        'structuralVariants',
    ]:
        for variant in upload_content[variant_section]:
            if 'variant' in variant:
                del variant['variant']
            if 'variantType' in variant:
                del variant['variantType']
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
    write_to_json: str = None,
    optional_content: Optional[Dict] = None,
) -> None:
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
        write_to_json: path to a JSON file to output the report upload body to if given on failure to upload
        optional_content: Pass-through content to include in the JSON upload
    """
    # set the default logging configuration
    logging.basicConfig(
        level=LOG_LEVELS[log_level],
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        datefmt='%m-%d-%y %H:%M:%S',
    )
    ipr_conn = ipr.IprConnection(username, password, ipr_url)
    graphkb_conn = GraphKBConnection()
    graphkb_conn.login(username, password)

    copy_variants = load_copy_variants(copy_variants_file) if copy_variants_file else []
    logger.info(f'loaded {len(copy_variants)} copy variants from: {copy_variants_file}')

    small_mutations = load_small_mutations(small_mutations_file) if small_mutations_file else []
    logger.info(f'loaded {len(small_mutations)} small mutations from: {small_mutations_file}')

    expression_variants = (
        load_expression_variants(expression_variants_file) if expression_variants_file else []
    )
    logger.info(
        f'loaded {len(expression_variants)} expression variants from: {expression_variants_file}'
    )

    structural_variants = (
        load_structural_variants(structural_variants_file) if structural_variants_file else []
    )
    logger.info(
        f'loaded {len(structural_variants)} structural variants from: {structural_variants_file}'
    )

    genes_with_variants = check_variant_links(
        small_mutations, expression_variants, copy_variants, structural_variants
    )

    # filter excess variants not required for extra gene information
    logger.info(f'annotating small mutations from: {small_mutations_file}')
    alterations = annotate_positional_variants(graphkb_conn, small_mutations, kb_disease_match)

    logger.info(f'annotating structural variants from: {structural_variants_file}')
    alterations.extend(
        annotate_positional_variants(graphkb_conn, structural_variants, kb_disease_match)
    )

    logger.info(f'annotating copy variants from {copy_variants_file}')
    alterations.extend(annotate_category_variants(graphkb_conn, copy_variants, kb_disease_match))

    logger.info(f'annotating expression variants from: {expression_variants_file}')
    alterations.extend(
        annotate_category_variants(graphkb_conn, expression_variants, kb_disease_match, False)
    )
    logger.info('fetching gene annotations')
    gene_information = get_gene_information(graphkb_conn, genes_with_variants)

    output = optional_content or dict()

    key_alterations, variant_counts = ipr.create_key_alterations(
        alterations, expression_variants, copy_variants, structural_variants, small_mutations
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
            'structuralVariants': [trim_empty_values(s) for s in structural_variants],
            'genes': gene_information,
            'genomicAlterationsIdentified': key_alterations,
            'variantCounts': variant_counts,
        }
    )
    for section in output:
        section_content_type = 'rows' if not isinstance(output[section], str) else 'characters'
        logger.info(f'section {section} has {len(output[section])} {section_content_type}')

    logger.info(f'made {graphkb_conn.request_count} requests to graphkb')

    output = clean_unsupported_content(output)

    try:
        result = ipr_conn.upload_report(output)
        logger.info(result)
    except Exception as err:
        if write_to_json:
            logging.info(f'writing report upload content to file: {write_to_json}')
            with open(write_to_json, 'w') as fh:
                import json

                fh.write(json.dumps(output))
        raise err
