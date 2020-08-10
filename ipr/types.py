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


class IprGeneVariant(TypedDict):
    gene: str
    key: str
    variantType: str
    variant: str


class IprGene(TypedDict):
    name: str
    cancerRelated: Optional[bool]
    knownFusionPartner: Optional[bool]
    knownSmallMutation: Optional[bool]
    tumourSuppressor: Optional[bool]
    oncogene: Optional[bool]
    therapeuticAssociated: Optional[bool]


class IprStructuralVariant(TypedDict):
    key: str
    variantType: str
    variant: str
    gene1: str
    gene2: str
    exon1: int
    exon2: int


class ImageDefinition(TypedDict):
    key: str
    path: str


IprVariant = Union[IprGeneVariant, IprStructuralVariant]
