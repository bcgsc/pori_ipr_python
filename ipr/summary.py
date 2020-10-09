from typing import Sequence, Dict, List, Set, Tuple
import base64
import json
from urllib.parse import urlencode


from graphkb.types import Statement, Record
from graphkb import GraphKBConnection
from graphkb.vocab import get_term_tree
from graphkb.constants import RELEVANCE_BASE_TERMS
from graphkb.statement import categorize_relevance
from graphkb.util import convert_to_rid_list

from .types import KbMatch, IprVariant
from .util import (
    convert_to_rid_set,
    get_preferred_drug_representation,
    get_preferred_gene_name,
    generate_ontology_preference_key,
)


OTHER_DISEASES = 'other disease types'
ENTREZ_GENE_URL = 'https://www.ncbi.nlm.nih.gov/gene'
# TODO: https://www.bcgsc.ca/jira/browse/DEVSU-1181
GRAPHKB_GUI = 'https://graphkb.bcgsc.ca'


def filter_by_record_class(
    record_list: List[Dict], *record_classes, exclude: bool = False
) -> List[Dict]:
    """
    Given a list of records, return the subset matching a class or list of classes
    """

    def check(name: str) -> bool:
        if exclude:
            return name not in record_classes
        else:
            return name in record_classes

    return [rec for rec in record_list if check(rec['@class'])]


def natural_join(word_list: List[str]) -> str:
    if len(word_list) > 1:
        return ', '.join(word_list[:-1]) + ', and ' + word_list[-1]
    return ''.join(word_list)


def natural_join_records(records: List[Dict], covert_to_word=lambda x: x['displayName']) -> str:
    word_list = sorted(list({covert_to_word(rec) for rec in records}))
    return natural_join(word_list)


def create_graphkb_link(
    record_ids: List[str],
    record_class: str = 'Statement',
) -> str:
    """
    Create a link for a set of statements to the GraphKB client
    """
    record_ids = sorted(list(set(record_ids)))
    if len(record_ids) == 1:
        return f'{GRAPHKB_GUI}/view/{record_class}/{record_ids[0].replace("#", "")}'
    complex_param = base64.b64encode(json.dumps({'target': record_ids}).encode("utf-8"))
    search_params = {'complex': complex_param, '@class': record_class}
    return f'{GRAPHKB_GUI}/data/table?{urlencode(search_params)}'


def substitute_sentence_template(
    template: str,
    conditions: List[Record],
    subjects: List[Record],
    relevance: Record,
    evidence: List[Record],
    statement_rids: List[str] = [],
    disease_matches: Set[str] = set(),
) -> str:
    """
    Create the filled-in sentence template for a given template and list of substitutions
    which may be the result of the aggregation of 1 or more statements
    """
    disease_conditions = filter_by_record_class(conditions, 'Disease')
    variant_conditions = filter_by_record_class(
        conditions, 'CategoryVariant', 'CatalogueVariant', 'PositionalVariant'
    )
    other_conditions = filter_by_record_class(
        conditions,
        'CategoryVariant',
        'CatalogueVariant',
        'PositionalVariant',
        'Disease',
        exclude=True,
    )
    result = template.replace(r'{relevance}', relevance['displayName'])

    def merge_diseases(diseases: List[Dict]) -> str:
        if len(convert_to_rid_set(diseases) - disease_matches) >= 2 and all(
            [d['@class'] == 'Disease' for d in diseases]
        ):
            words = sorted(
                list(set([s['displayName'] for s in diseases if s['@rid'] in disease_matches]))
            )
            words.append(OTHER_DISEASES)
            return natural_join(words)
        else:
            return natural_join_records(diseases)

    if r'{subject}' in template:
        # remove subject from the conditions replacements
        subjects_ids = convert_to_rid_set(subjects)
        disease_conditions = [d for d in disease_conditions if d['@rid'] not in subjects_ids]
        variant_conditions = [d for d in variant_conditions if d['@rid'] not in subjects_ids]
        other_conditions = [d for d in other_conditions if d['@rid'] not in subjects_ids]

        result = result.replace(r'{subject}', merge_diseases(subjects))

    if r'{conditions:disease}' in template:
        result = result.replace(r'{conditions:disease}', merge_diseases(disease_conditions))
    else:
        other_conditions.extend(disease_conditions)

    if r'{conditions:variant}' in template:
        result = result.replace(r'{conditions:variant}', natural_join_records(variant_conditions))
    else:
        other_conditions.extend(variant_conditions)

    result = result.replace(r'{conditions}', natural_join_records(other_conditions))

    link_url = create_graphkb_link(statement_rids) if statement_rids else ''

    if r'{evidence}' in template:

        evidence_str = ', '.join(sorted(list({e['displayName'] for e in evidence})))
        if link_url:
            evidence_str = f'<a href="{link_url}" target="_blank" rel="noopener">{evidence_str}</a>'
        result = result.replace(r'{evidence}', evidence_str)

    return result


def aggregate_statements(
    graphkb_conn: GraphKBConnection,
    template: str,
    statements: List[Statement],
    disease_matches: Set[str],
) -> Dict[str, str]:
    """
    Group Statements that only differ in disease conditions and evidence
    """
    hash_other: Dict[Tuple, List[Statement]] = {}

    def generate_key(statement: Statement) -> Tuple:
        result = [
            cond['displayName']
            for cond in filter_by_record_class(statement['conditions'], 'Disease', exclude=True)
            if cond['@rid'] != statement['subject']['@rid']
        ]
        if statement.get('subject', {}).get('@class', 'Disease') != 'Disease':
            subject = statement['subject']
            if subject['@class'] == 'Therapy':
                alt = get_preferred_drug_representation(graphkb_conn, subject['@rid'])
                statement['subject'] = alt
            result.append(statement['subject']['displayName'])
        result.append(statement['relevance']['displayName'])
        result.append(statement['displayNameTemplate'])
        return tuple(sorted(set(result)))

    for statement in statements:
        key = generate_key(statement)
        hash_other.setdefault(key, []).append(statement)

    result = {}
    for key, group in hash_other.items():
        conditions = []
        subjects = []
        evidence = []
        relevance = group[0]['relevance']
        template = group[0]['displayNameTemplate']
        for statement in group:
            conditions.extend(statement['conditions'])
            evidence.extend(statement['evidence'])
            subjects.append(statement['subject'])

        sentence = substitute_sentence_template(
            template,
            conditions,
            subjects,
            relevance,
            evidence,
            statement_rids=convert_to_rid_list(group),
            disease_matches=disease_matches,
        )

        for statement in group:
            result[statement['@rid']] = sentence
    return result


def display_variants(gene_name: str, variants: List[IprVariant]):
    def display_variant(variant: IprVariant):
        if 'gene' in variant and 'proteinChange' in variant:
            return f'{variant["gene"]}:{variant["proteinChange"]}'
        if 'gene1' in variant and 'gene2' in variant:
            return f'({variant["gene1"]},{variant["gene2"]}):fusion(e.{variant.get("exon1", "?")},e.{variant.get("exon2", "?")})'
        if 'kbCategory' in variant and variant['kbCategory']:
            return f'{variant["kbCategory"]} of {variant["gene"]}'

        raise ValueError(f'Unable to form display_variant of {variant["variant"]}: {variant}')

    result = sorted(list({v for v in [display_variant(e) for e in variants] if gene_name in v}))
    variants_text = natural_join(result)
    if len(result) > 1:
        return (
            f'Multiple variants of the gene {gene_name} were observed in this case: {variants_text}'
        )
    elif result:
        return f'{variants_text[0].upper()}{variants_text[1:]} was observed in this case.'
    return ''


def create_section_html(
    graphkb_conn: GraphKBConnection,
    gene_name: str,
    sentences_by_statement_id: Dict[str, str],
    statements: Dict[str, Statement],
    exp_variants: List[IprVariant],
) -> str:
    """
    Generate HTML for a gene section of the comments
    """
    output = [f'<h2>{gene_name}</h2>']

    sentence_categories: Dict[str, str] = {}

    for statement_id, sentence in sentences_by_statement_id.items():
        relevance = statements[statement_id]['relevance']['@rid']
        category = categorize_relevance(
            graphkb_conn, relevance, RELEVANCE_BASE_TERMS + [('resistance', ['no sensitivity'])]
        )
        sentence_categories[sentence] = category

    # get the entrez gene description
    genes = sorted(
        graphkb_conn.query(
            {
                'target': 'Feature',
                'filters': {
                    'AND': [
                        {'source': {'target': 'Source', 'filters': {'name': 'entrez gene'}}},
                        {'name': gene_name},
                        {'biotype': 'gene'},
                    ]
                },
            },
        ),
        key=generate_ontology_preference_key,
    )

    variants_text = display_variants(gene_name, exp_variants)
    if not variants_text:
        # exclude sections where they are not linked to an experimental variant. this can occur when there are co-occurent statements collected
        return ''
    if genes and genes[0].get('description', ''):
        description = '. '.join(genes[0]['description'].split('. ')[:2])
        sourceId = genes[0]['sourceId']

        output.append(
            f'''
<blockquote class="entrez_description" cite="{ENTREZ_GENE_URL}/{sourceId}">
    {description}.
</blockquote>
<p>
    {variants_text}
</p>
'''
        )

    sentences_used: Set[str] = set()

    for section in [
        {s for (s, v) in sentence_categories.items() if v == 'diagnostic'},
        {s for (s, v) in sentence_categories.items() if v == 'biological'},
        {s for (s, v) in sentence_categories.items() if v in ['therapeutic', 'prognostic']},
        {
            s
            for (s, v) in sentence_categories.items()
            if v not in ['diagnostic', 'biological', 'therapeutic', 'prognostic', 'resistance']
        },
        {s for (s, v) in sentence_categories.items() if v == 'resistance'},
    ]:

        content = '. '.join(sorted(list(section - sentences_used)))
        sentences_used.update(section)
        output.append(f'<p>{content}</p>')
    return '\n'.join(output)


def section_statements_by_genes(
    graphkb_conn: GraphKBConnection, statements: Sequence[Statement]
) -> Dict[str, Set[str]]:
    """
    Determine the statements associated with each gene name
    """
    genes: Dict[str, Set[str]] = {}

    for statement in statements:
        for condition in statement['conditions']:
            if condition.get('biotype', '') == 'gene':
                gene = get_preferred_gene_name(graphkb_conn, condition['@rid'])
                genes.setdefault(gene, set()).add(statement['@rid'])
            else:
                if condition.get('reference1', ''):
                    gene = get_preferred_gene_name(graphkb_conn, condition['reference1'])
                    genes.setdefault(gene, set()).add(statement['@rid'])

                if condition.get('reference2', ''):
                    gene = get_preferred_gene_name(graphkb_conn, condition['reference2'])
                    genes.setdefault(gene, set()).add(statement['@rid'])
    return genes


def summarize(
    graphkb_conn: GraphKBConnection,
    matches: Sequence[KbMatch],
    disease_name: str,
    variants: List[IprVariant],
) -> str:
    """
    Given a list of GraphKB matches generate a text summary to add to the report
    """
    templates: Dict[str, List[Statement]] = {}
    statements: Dict[str, Statement] = {}
    variants_by_keys = {v['key']: v for v in variants}
    variant_keys_by_statement_ids: Dict[str, Set[str]] = {}

    for match in matches:
        rid = match['kbStatementId']
        exp_variant = match['variant']
        variant_keys_by_statement_ids.setdefault(rid, set()).add(exp_variant)

    exp_variants_by_statements: Dict[str, List[IprVariant]] = {}
    for rid, keys in variant_keys_by_statement_ids.items():
        exp_variants_by_statements[rid] = [variants_by_keys[key] for key in keys]

    disease_matches = convert_to_rid_set(
        get_term_tree(graphkb_conn, disease_name, ontology_class='Disease')
    )

    # get details for statements
    for match in matches:
        rid = match['kbStatementId'].replace('#', '')
        result = graphkb_conn.request(f'/statements/{rid}?neighbors=1')['result']

        templates.setdefault(result['displayNameTemplate'], []).append(result)
        statements[result['@rid']] = result

    # aggregate similar sentences
    sentences = {}
    for template, group in templates.items():
        sentences.update(aggregate_statements(graphkb_conn, template, group, disease_matches))

    # section statements by genes
    statements_by_genes = section_statements_by_genes(graphkb_conn, list(statements.values()))

    output: List[str] = [
        '<h3>The comments below were automatically generated from matches to GraphKB and have not been manually reviewed</h3>'
    ]

    for section, statement_rids in sorted(
        statements_by_genes.items(), key=lambda x: len(x[1]), reverse=True
    ):
        exp_variants = {}
        for variant_list in [exp_variants_by_statements[r] for r in statement_rids]:
            for variant in variant_list:
                exp_variants[variant['key']] = variant

        output.append(
            create_section_html(
                graphkb_conn,
                section,
                {r: sentences[r] for r in statement_rids},
                {r: statements[r] for r in statement_rids},
                list(exp_variants.values()),
            )
        )

    return '\n'.join(output)
