from typing import Optional, Union

try:
    from typing import TypedDict  # type: ignore
except ImportError:
    from typing_extensions import TypedDict


class KbMatch(TypedDict):
    variant: str
    variantType: str
    approvedTherapy: bool
    category: str
    context: str
    kbContextId: str
    disease: str
    evidenceLevel: str
    kbStatementId: str
    kbVariant: str
    kbVariantId: str
    matchedCancer: bool
    reference: str
    relevance: str
    kbRelevanceId: str
    externalSource: str
    externalStatementId: str
    reviewStatus: str


class IprGene(TypedDict):
    name: str
    cancerRelated: Optional[bool]
    knownFusionPartner: Optional[bool]
    knownSmallMutation: Optional[bool]
    tumourSuppressor: Optional[bool]
    oncogene: Optional[bool]
    therapeuticAssociated: Optional[bool]


class IprGeneVariant(TypedDict):
    gene: str
    key: str
    variantType: str
    variant: str


class IprCopyVariant(IprGeneVariant):
    # variantType == 'cnv'
    kbCategory: str
    cnvState: str


class IprExprVariant(IprGeneVariant):
    # variantType == 'exp'
    kbCategory: str
    expressionState: str
    histogramImage: Optional[str]


class IprSmallMutationVariant(IprGeneVariant):
    germline: Optional[bool]
    hgvsCds: Optional[str]
    hgvsGenomic: Optional[str]
    hgvsProtein: Optional[str]
    startPosition: Optional[int]
    endPosition: Optional[int]
    normalAltCount: Optional[int]
    normalDepth: Optional[int]
    normalRefCount: Optional[int]
    rnaAltCount: Optional[int]
    rnaDepth: Optional[int]
    rnaRefCount: Optional[int]
    tumourAltCount: Optional[int]
    tumourDepth: Optional[int]
    tumourRefCount: Optional[int]


class IprFusionVariant(TypedDict):
    key: str
    variantType: str
    variant: str
    gene1: str
    gene2: str
    exon1: int
    exon2: int
    svg: Optional[str]  # path to svg image of fusion


class ImageDefinition(TypedDict):
    key: str
    path: str


IprStructuralVariant = Union[IprSmallMutationVariant, IprFusionVariant]
IprVariant = Union[IprCopyVariant, IprExprVariant, IprStructuralVariant]
