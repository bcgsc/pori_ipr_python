import argparse
import json
import logging
import os

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
from .util import logger, LOG_LEVELS


def file_path(path):
    if not os.path.exists(path):
        raise argparse.ArgumentTypeError(f'{repr(path)} is not a valid filename. does not exist')
    return path


def command_interface():
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
    parser.add_argument('--log_level', default='info', choices=LOG_LEVELS.keys())

    # TODO: upload JSON to IPR instead of writing output
    parser.add_argument(
        '-o',
        '--output_json',
        default='ipr_input.json',
        help='file path to write the output json content to',
    )

    args = parser.parse_args()

    main(args)


def main(args):
    # set the default logging configuration
    logging.basicConfig(
        level=LOG_LEVELS[args.log_level],
        format='%(asctime)s %(name)s %(levelname)s %(message)s',
        datefmt='%m-%d-%y %H:%M:%S',
    )

    graphkb_conn = GraphKBConnection()
    graphkb_conn.login(args.username, args.password)
    disease_name = 'colorectal cancer'

    copy_variants = load_copy_variants(args.copy_variants) if args.copy_variants else []
    logger.info(f'loaded {len(copy_variants)} copy variants')

    small_mutations = load_small_mutations(args.small_mutations) if args.small_mutations else []
    logger.info(f'loaded {len(small_mutations)} small mutations')

    expression_variants = (
        load_expression_variants(args.expression_variants) if args.expression_variants else []
    )
    logger.info(f'loaded {len(expression_variants)} expression variants')

    structural_variants = (
        load_structural_variants(args.structural_variants) if args.structural_variants else []
    )
    logger.info(f'loaded {len(structural_variants)} structural variants')

    genes_with_variants = check_variant_links(
        small_mutations, expression_variants, copy_variants, structural_variants
    )

    # filter excess variants not required for extra gene information
    alterations = []  # annotate_positional_variants(graphkb_conn, small_mutations, disease_name)
    logger.info('annotating structural variants')
    alterations.extend(
        annotate_positional_variants(graphkb_conn, structural_variants, disease_name)
    )

    logger.info('annotating copy variants')
    alterations.extend(annotate_category_variants(graphkb_conn, copy_variants, disease_name))

    logger.info('annotating expression variants')
    alterations.extend(
        annotate_category_variants(graphkb_conn, expression_variants, disease_name, False)
    )
    logger.info('fetching gene annotations')
    gene_information = get_gene_information(graphkb_conn, genes_with_variants)
    # TODO: Append gene level information to each variant type (until IPR does this itself?)

    logger.info(f'writing: {args.output_json}')
    with open(args.output_json, 'w') as fh:
        output = {
            'alterations': alterations,
            'cnv': [c for c in copy_variants if c['gene'] in genes_with_variants],
            'smallMutations': small_mutations,
            'outliers': [e for e in expression_variants if e['gene'] in genes_with_variants],
            'sv': structural_variants,
            'genes': gene_information,
        }
        for section in output:
            logger.info(f'section {section} has {len(output[section])} rows')
        fh.write(json.dumps(output, indent='  ', sort_keys=True,))
    logger.info(f'made {graphkb_conn.request_count} requests to graphkb')
    # TODO: upload to IPR
