import argparse
import pandas as pd

MIN_COSMIC_COUNT = 3

GENE = 'Hugo_Symbol'
PROTEIN_CHANGE = 'Protein_Change'
VARIANT_CLASSIFICATION = 'Variant_Classification'
CHROMOSOME = 'Chromosome'
START_POSITION = 'Start_position'
REF_ALLELE = 'Reference_Allele'
ALT_ALLELE = 'Tumor_Seq_Allele2'
REF_COUNT = 'REF_COUNT'
ALT_COUNT = 'ALT_COUNT'

COSMIC_GENE = 'Gene name'
COSMIC_AA = 'Mutation AA'

EXAC_CHR = 'CHROM'
EXAC_POS = 'POS'
EXAC_REF = 'REF'
EXAC_ALT = 'ALT'
EXAC_AC = 'AC'

_maf_column_map = {
    GENE: 'gene',
    PROTEIN_CHANGE: 'protein_change',
    VARIANT_CLASSIFICATION: 'variant_classification',
    CHROMOSOME: 'chromosome',
    START_POSITION: 'start_position',
    REF_ALLELE: 'ref_allele',
    ALT_ALLELE: 'alt_allele'
}

_cosmic_column_map = {
    COSMIC_GENE: 'gene',
    COSMIC_AA: 'protein_change',
    'count': 'count'
}

_exac_column_map = {
    EXAC_CHR: 'chromosome',
    EXAC_POS: 'start_position',
    EXAC_REF: 'ref_allele',
    EXAC_ALT: 'exac_alt_allele',
    EXAC_AC: 'exac_ac'
}


def check_column_names(df, map):
    for column_name in map.keys():
        assert column_name in df.columns, \
            'Expected column %s not found among %s' % (column_name, df.columns)


def read(handle, **kwargs):
    return pd.read_csv(handle, sep='\t', comment='#', dtype='object', **kwargs)


def standard_read(handle, column_map, **kwargs):
    check_column_names(read(handle, nrows=3), column_map)
    return read(handle, usecols=column_map.keys(), **kwargs).rename(columns=column_map)


parser = argparse.ArgumentParser()
parser.add_argument('--maf', help='MAF to annotate and filter', required=True)
parser.add_argument('--cosmic', help='COSMIC mutation export', required=True)
parser.add_argument('--cosmic_version', help='COSMIC version', required=True)
parser.add_argument('--exac', help='formatted exac', required=True)
parser.add_argument('--filter_syn', help='Removes syn variants. True/False', default=False)
parser.add_argument('--filter_read_depth', help='Filters based on specified read depth. Int.', default=0)
args = parser.parse_args()

inputs = {
    'maf_handle': args.maf,
    'cosmic_handle': args.cosmic,
    'cosmic_version': args.cosmic_version,
    'exac': args.exac,
    'filter_syn': args.filter_syn,
    'min_depth': args.filter_read_depth
}

df = read(inputs['maf_handle'], low_memory=False).rename(columns=_maf_column_map)
df_cosmic = standard_read(inputs['cosmic_handle'], _cosmic_column_map, low_memory=False)
df_exac = standard_read(inputs['exac'], _exac_column_map, low_memory=False)

if inputs['filter_syn']:
    variants = ['Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site']
    df = df[df['variant_classification'].isin(variants)]

if int(inputs['min_depth']) != 0:
    df['read_depth'] = df[ALT_COUNT].astype(int) + df[REF_COUNT].astype(int)
    df = df[df['read_depth'].astype(int) >= inputs['min_depth']]

df.index = range(0, len(df))

df = df.merge(df_cosmic,
              on=[_cosmic_column_map[COSMIC_GENE], _cosmic_column_map[COSMIC_AA]], how='left')

df = df.merge(df_exac,
              on=[_exac_column_map[EXAC_CHR], _exac_column_map[EXAC_POS], _exac_column_map[EXAC_REF]], how='left')

df.loc[:, _exac_column_map[EXAC_ALT]] = df.loc[:, _exac_column_map[EXAC_ALT]].fillna('')
df.loc[:, 'count'] = df.loc[:, 'count'].fillna(0)


for idx in df.index:
    boolean = df.loc[idx, _maf_column_map[ALT_ALLELE]] in df.loc[idx, _exac_column_map[EXAC_ALT]]
    if not boolean:
        df.loc[idx, _exac_column_map[EXAC_AC]] = 0

cosmic_count_col = 'COSMIC_v' + inputs['cosmic_version'] + '_counts'
df = df.rename(columns={'count': cosmic_count_col})

idx_cosmic = (df[cosmic_count_col].astype(int) >= int(MIN_COSMIC_COUNT))
idx_exac = (df[_exac_column_map[EXAC_AC]].astype(str) == str(0))
idx_filter = idx_exac | idx_cosmic

df_pass_filter = df.loc[idx_filter, :]
df_fail_filter = df.loc[~idx_filter, :]

df_pass_filter.to_csv('df_filter_pass.txt', sep='\t', index=False)
df_fail_filter.to_csv('df_filter_fail.txt', sep='\t', index=False)