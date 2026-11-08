"""Every external data source Koel uses, with its licence and access terms.

Rule: only sources marked ``redistributable=True`` may be bundled into outputs
that Koel publishes. Anything else stays local and is used only for analysis
under its own terms.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Source:
    key: str
    group: str
    url: str
    filename: str
    licence: str
    redistributable: bool
    approx_size: str
    notes: str = ""
    tags: tuple[str, ...] = field(default_factory=tuple)


CIVIC_BASE = "https://civicdb.org/downloads/nightly"

SOURCES: list[Source] = [
    Source(
        key="civic-evidence",
        group="civic",
        url=f"{CIVIC_BASE}/nightly-ClinicalEvidenceSummaries.tsv",
        filename="civic/ClinicalEvidenceSummaries.tsv",
        licence="CC0 1.0",
        redistributable=True,
        approx_size="~20 MB",
        notes="Expert-curated evidence items, including oncogenicity evidence.",
    ),
    Source(
        key="civic-variants",
        group="civic",
        url=f"{CIVIC_BASE}/nightly-VariantSummaries.tsv",
        filename="civic/VariantSummaries.tsv",
        licence="CC0 1.0",
        redistributable=True,
        approx_size="~5 MB",
        notes="Variant coordinates and protein changes.",
    ),
    Source(
        key="civic-molecular-profiles",
        group="civic",
        url=f"{CIVIC_BASE}/nightly-MolecularProfileSummaries.tsv",
        filename="civic/MolecularProfileSummaries.tsv",
        licence="CC0 1.0",
        redistributable=True,
        approx_size="~5 MB",
        notes="Links evidence items to variants.",
    ),
    Source(
        key="civic-assertions",
        group="civic",
        url=f"{CIVIC_BASE}/nightly-AssertionSummaries.tsv",
        filename="civic/AssertionSummaries.tsv",
        licence="CC0 1.0",
        redistributable=True,
        approx_size="<1 MB",
        notes="AMP/ASCO/CAP and ClinGen/CGC/VICC-style assertions.",
    ),
    Source(
        key="alphamissense-aa",
        group="alphamissense",
        url="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz",
        filename="alphamissense/AlphaMissense_aa_substitutions.tsv.gz",
        licence="CC BY-NC-SA 4.0",
        redistributable=False,
        approx_size="~1.2 GB",
        notes=(
            "Scores for every possible single amino-acid substitution in the "
            "canonical human proteome. Not validated or approved for clinical use. "
            "The downloaded file header (checked 8 Nov 2026) states CC BY-NC-SA 4.0: "
            "non-commercial use only, and anything derived from the scores must be "
            "shared under the same licence. Koel outputs that include these scores "
            "must carry that licence."
        ),
    ),
]

SOURCES.append(
    Source(
        key="uniprot-mapping",
        group="uniprot",
        url="https://rest.uniprot.org/uniprotkb/search",
        filename="uniprot/gene_to_uniprot.tsv",
        licence="CC BY 4.0",
        redistributable=True,
        approx_size="<1 MB",
        notes="Gene symbol to reviewed human UniProt accession, queried per gene.",
    )
)

SOURCES.append(
    Source(
        key="gnomad-v4-api",
        group="gnomad",
        url="https://gnomad.broadinstitute.org/api",
        filename="gnomad/<GENE>.json",
        licence="gnomAD: released without access restrictions (verify current terms)",
        redistributable=True,
        approx_size="~1-20 MB per gene",
        notes="Population allele counts per gene via the public GraphQL API; paced to respect rate limits.",
    )
)

SOURCES.append(
    Source(
        key="esm1b-llr",
        group="esm1b",
        url="https://huggingface.co/spaces/ntranoslab/esm_variants/resolve/main/ALL_hum_isoforms_ESM1b_LLR.zip",
        filename="esm1b/ALL_hum_isoforms_ESM1b_LLR.zip",
        licence="CC BY-NC 4.0",
        redistributable=False,
        approx_size="~1.3 GB",
        notes=(
            "Precomputed ESM1b log-likelihood ratios for all human protein isoforms "
            "(Brandes et al., Nature Genetics 2023). Non-commercial; licence from the "
            "Hugging Face space card, checked 8 Nov 2026. Lower LLR = more damaging."
        ),
    )
)

SOURCES.append(
    Source(
        key="intogen-drivers",
        group="intogen",
        url="https://www.intogen.org/download?file=IntOGen-Drivers-20240920.zip",
        filename="intogen/IntOGen-Drivers-20240920.zip",
        licence="CC0 1.0",
        redistributable=True,
        approx_size="~1 MB",
        notes="Compendium of cancer driver genes with mode of action (Act = oncogene, LoF = tumour suppressor).",
    )
)

# Sources still to be chosen, with the licence question that must be settled first.
PENDING = {
    "structures": "AlphaFold DB (CC BY 4.0) for residue coordinates and pLDDT.",
}


def by_group(group: str) -> list[Source]:
    return [s for s in SOURCES if s.group == group or s.key == group]
