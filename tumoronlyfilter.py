import argparse
import pandas as pd

GENE = 'Hugo_Symbol'
PROTEIN = 'Protein_Change'
CHROMOSOME = 'Chromosome'
START_POSITION = 'Start_position'
REF_ALLELE = 'Reference.Allele'
ALT_ALLELE = 'Alternate.Allele'
REF_COUNT = 'REF_COUNT'
ALT_COUNT = 'ALT_COUNT'

EXAC_CHR = 'CHROM'
EXAC_POS = 'POS'
EXAC_REF = 'REF'
EXAC_ALT = 'ALT'
EXAC_AF = 'AF'
EXAC_AC = 'AC'
EXAC_AC_AFR = 'AC_AFR'
EXAC_AC_AMR = 'AC_AMR'
EXAC_AC_EAS = 'AC_EAS'
EXAC_AC_FIN = 'AC_FIN'
EXAC_AC_NFE = 'AC_NFE'
EXAC_AC_OTH = 'AC_OTH'
EXAC_AC_SAS = 'AC_SAS'
EXAC_AN = 'AN'
EXAC_AN_AFR = 'AN_AFR'
EXAC_AN_AMR = 'AN_AMR'
EXAC_AN_EAS = 'AN_EAS'
EXAC_AN_FIN = 'AN_FIN'
EXAC_AN_NFE = 'AN_NFE'
EXAC_AN_OTH = 'AN_OTH'
EXAC_AN_SAS = 'AN_SAS'

MAPPED_GENE = 'gene'
MAPPED_CHR = 'chromosome'
MAPPED_REF = 'ref_allele'
MAPPED_ALT = 'alt_allele'
MAPPED_POS = 'start_position'
MAPPED_AA = 'protein_change'

maf_column_map = {
    GENE: MAPPED_GENE,
    CHROMOSOME: MAPPED_CHR,
    PROTEIN: MAPPED_AA,
    START_POSITION: MAPPED_POS,
    REF_ALLELE: MAPPED_REF,
    ALT_ALLELE: MAPPED_ALT
}

output_column_map = {v: k for k, v in maf_column_map.items()}

exac_column_map = {
    EXAC_CHR: MAPPED_CHR,
    EXAC_POS: MAPPED_POS,
    EXAC_REF: MAPPED_REF,
    EXAC_ALT: MAPPED_ALT,
    EXAC_AF: 'exac_af',
    EXAC_AC: 'exac_ac',
    EXAC_AC_AFR: 'exac_ac_afr',
    EXAC_AC_AMR: 'exac_ac_amr',
    EXAC_AC_EAS: 'exac_ac_eas',
    EXAC_AC_FIN: 'exac_ac_fin',
    EXAC_AC_NFE: 'exac_ac_nfe',
    EXAC_AC_OTH: 'exac_ac_oth',
    EXAC_AC_SAS: 'exac_ac_sas',
    EXAC_AN: 'exac_an',
    EXAC_AN_AFR: 'exac_an_afr',
    EXAC_AN_AMR: 'exac_an_amr',
    EXAC_AN_EAS: 'exac_an_eas',
    EXAC_AN_FIN: 'exac_an_fin',
    EXAC_AN_NFE: 'exac_an_nfe',
    EXAC_AN_OTH: 'exac_an_oth',
    EXAC_AN_SAS: 'exac_an_sas',

}

whitelist_column_map = {0: MAPPED_CHR, 1: MAPPED_POS, 2: 'End_position', 3:'alteration'}

population_keys = [EXAC_AC_AFR, EXAC_AC_AMR, EXAC_AC_EAS, EXAC_AC_FIN, EXAC_AC_NFE, EXAC_AC_OTH, EXAC_AC_SAS]
populations = [exac_column_map[x] for x in population_keys]

def check_column_names(df, map):
    for column_name in map.keys():
        assert column_name in df.columns, \
            'Expected column %s not found among %s' % (column_name, df.columns)


def read(handle, **kwargs):
    return pd.read_csv(handle, sep='\t', comment='#', dtype='object', **kwargs)


def standard_read(handle, column_map, **kwargs):
    check_column_names(read(handle, nrows=3), column_map)
    return read(handle, **kwargs).rename(columns=column_map)


parser = argparse.ArgumentParser()
parser.add_argument('--maf', help='MAF to annotate and filter', required=True)
parser.add_argument('--exac', help='formatted exac', required=True)
parser.add_argument('--whitelist', help='whitelist for somatic sites', default=False)
parser.add_argument('--min_exac_ac',
                    help='Minimum allele count across any population to be considered germline', default=10)
parser.add_argument('--filter_syn', help='Removes syn variants. True/False', action='store_true')
parser.add_argument('--filter_read_depth', help='Filters based on specified read depth. Int.', default=0)
args = parser.parse_args()

inputs = {
    'maf_handle': args.maf,
    'exac_handle': args.exac,
    'whitelist_handle': args.whitelist,
    'min_exac_ac': args.min_exac_ac,
    'filter_syn': args.filter_syn,
    'min_depth': args.filter_read_depth
}

df = standard_read(inputs['maf_handle'], maf_column_map, low_memory=False)
exac = standard_read(inputs['exac_handle'], exac_column_map, low_memory=False)

idx_original = df.index

if inputs['filter_syn']:
    variants = ['Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site']
    df = df[df['variant_classification'].isin(variants)]

if int(inputs['min_depth']) != 0:
    df['read_depth'] = df[ALT_COUNT].astype(int) + df[REF_COUNT].astype(int)
    df = df[df['read_depth'].astype(int) >= inputs['min_depth']]

merge_cols = [MAPPED_CHR, MAPPED_POS, MAPPED_REF, MAPPED_ALT]
df = df.merge(exac, on=merge_cols, how='left')

idx_exac = df[(df.loc[:, populations].fillna(0).astype(int) > inputs['min_exac_ac']).sum(axis=1) == int(0)].index

if inputs['whitelist_handle']:
    whitelist = read(inputs['whitelist_handle'], header=-1).rename(columns=whitelist_column_map)
    df[whitelist_column_map[3]] = df[MAPPED_GENE].astype(str) + ':' + df[MAPPED_AA].str.split('p.', expand=True).loc[:,1].astype(str)
    df[whitelist_column_map[3]] = df[whitelist_column_map[3]].fillna('')
    idx_whitelist = df[df[whitelist_column_map[3]].isin(whitelist[whitelist_column_map[3]])].index

    idx_passed = idx_exac.union(idx_whitelist)
    df = df.drop([whitelist_column_map[3]], axis=1)
else:
    idx_passed = idx_exac

idx_failed = idx_original.difference(idx_passed)


df = df.rename(columns=output_column_map)
df = df.fillna('')

df_passed = df.loc[idx_passed, :]
df_failed = df.loc[idx_failed, :]

df_passed.to_csv('tumoronlyfilter_passed.txt', sep='\t', index=False)
df_failed.to_csv('tumoronlyfilter_failed.txt', sep='\t', index=False)
