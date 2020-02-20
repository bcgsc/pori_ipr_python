import os
import json
import argparse

from argparse_env import ArgumentParser, Action
from graphkb import GraphKBConnection

from .inputs import load_copy_variants, load_small_mutations
from .annotate import annotate_variants


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
    conn = GraphKBConnection()
    conn.login(args.username, args.password)

    copy_variant_genes = set()
    genes_with_variants = set()  # filter excess copy variants

    copy_variants = load_copy_variants(args.copy_variants) if args.copy_variants else []
    print(f'loaded {len(copy_variants)} copy variants')

    small_mutations = load_small_mutations(args.small_mutations) if args.small_mutations else []
    print(f'loaded {len(small_mutations)} small mutations')

    # filter excess variants not required for extra gene information

    for variant in copy_variants:
        if variant['variant']:
            genes_with_variants.add(variant['gene'])
        copy_variant_genes.add(variant['gene'])

    for variant in small_mutations:
        gene = variant['gene']
        if gene not in copy_variant_genes:
            raise ValueError(
                f'gene ({gene}) has a small mutation but is missing copy number information'
            )
        genes_with_variants.add(gene)

    alterations = annotate_variants(conn, small_mutations=small_mutations)
    alterations.extend(annotate_variants(conn, copy_variants=copy_variants))

    # TODO: Append gene level information to each variant type (until IPR does this itself)

    print(f'writing: {args.output_json}')
    with open(args.output_json, 'w') as fh:
        fh.write(
            json.dumps(
                {
                    'alterations': alterations,
                    'cnv': [c for c in copy_variants if c['gene'] in genes_with_variants],
                    'smallMutations': small_mutations,
                },
                indent='  ',
                sort_keys=True,
            )
        )
    print(f'made {conn.request_count} requests to graphkb')
    # TODO: upload to IPR
