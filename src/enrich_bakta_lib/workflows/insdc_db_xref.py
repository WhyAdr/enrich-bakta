"""Database names in the INSDC cross-reference vocabulary snapshot."""

VOCABULARY_ID = "insdc-db-xref-2026-11-25"
VOCABULARY_DATE = "2026-11-25"
FEATURE_TABLE_VERSION = "11.4"
FEATURE_TABLE_URL = "https://www.insdc.org/submitting-standards/feature-table/"
VOCABULARY_URL = (
    "https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/"
)
NCBI_DB_XREF_URL = "https://www.ncbi.nlm.nih.gov/genbank/collab/db_xref/"

# The INSDC table's qualifier examples give the namespace when display names differ.
# The NCBI table and validator registry corroborate InterPro, TAIR, and Phytozome.
# The registry also confirms BOLD despite the mixed-case example in both tables.
DATABASE_PREFIXES = frozenset(
    """AFTOL
APHIDBASE
ASAP
ATCC
ATCC(dna)
ATCC(in host)
AceView/WormGenes
AntWeb
ApiDB
ApiDB_CryptoDB
ApiDB_PlasmoDB
ApiDB_ToxoDB
Araport
Axeldb
BDGP_EST
BDGP_INS
BEEBASE
BEETLEBASE
BEI
BGD
BOLD
CABRI
CCAP
CDD
CGD
ECOCYC
ENSEMBL
EPD
ERIC
ESTLIB
EcoGene
EnsemblGenomes
FANTOM_DB
FBOL
FLYBASE
Fungorum
GABI
GDB
GI
GO
GOA
GRIN
GeneDB
GeneID
Greengenes
H-InvDB
HGNC
HMP
HOMD
HPM
HSSP
IKMC
IMGT/GENE-DB
IMGT/HLA
IMGT/LIGM
IRD
ISFinder
ISHAM-ITS
InterPro
IntrepidBio
JCM
JGIDB
LocusID
MGI
MIM
MaizeGDB
MarpolBase
MedGen
MycoBank
NBRC
NMPDR
NRESTdb
NextDB
OrthoMCL
Osa1
PBmice
PDB
PFAM
PGN
PIR
PSEUDO
Pathema
Phytozome
PomBase
PseudoCap
RAP-DB
RATMAP
RBGE_garden
RBGE_herbarium
RFAM
RGD
RNAcentral
RZPD
RiceGenes
SEED
SGD
SGN
SK-FST
SRPDB
SoyBase
SubtiList
TAIR
TIGRFAM
TubercuList
UNILIB
UNITE
UniProtKB/Swiss-Prot
UniProtKB/TrEMBL
UniSTS
VBASE2
VGNC
VectorBase
ViPR
WorfDB
WormBase
Xenbase
ZFIN
dbEST
dbProbe
dbSNP
dbSTS
dictyBase
miRBase
niaEST
taxon""".splitlines()
)
