import argparse
import pandas as pd
import numpy as np

GENE = 'Hugo_Symbol'
PROTEIN = 'Protein_Change'
CHROMOSOME = 'Chromosome'
ALT = 'Alteration'
START_POSITION = 'Start_position'
END_POSITION = 'End_position'
REF_ALLELE = 'Reference_Allele'
ALT_ALLELE = 'Tumor_Seq_Allele2'
REF_COUNT = 't_ref_count'
ALT_COUNT = 't_alt_count'
VAR_CLASS = 'Variant_Classification'

maf_handle = 'maf_handle'
exac_handle = 'exac_handle'
whitelist_handle = 'whitelist_handle'
filter_syn = 'filter_syn'
min_exac_ac = 'min_exac_ac'
min_depth = 'min_depth'

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
MAPPED_VAR_CLASS = 'variant_classification'
MAPPED_REF_COUNT = 'ref_count'
MAPPED_ALT_COUNT = 'alt_count'

EXAC_COMMON = 'exac_common'
WL = 'whitelist'
DEPTH = 'read_depth'
LOW_DEPTH = 'low_read_depth'
CODING = 'coding'
COMMON = 'common_variant'

maf_column_map = {
    GENE: MAPPED_GENE,
    CHROMOSOME: MAPPED_CHR,
    PROTEIN: MAPPED_AA,
    START_POSITION: MAPPED_POS,
    REF_ALLELE: MAPPED_REF,
    ALT_ALLELE: MAPPED_ALT,
    VAR_CLASS: MAPPED_VAR_CLASS,
    REF_COUNT: MAPPED_REF_COUNT,
    ALT_COUNT: MAPPED_ALT_COUNT
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

whitelist_column_map = {0: MAPPED_CHR, 1: MAPPED_POS, 2: END_POSITION, 3:ALT}

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
    return read(handle, encoding='latin-1', **kwargs).rename(columns=column_map)

def apply_str(x):
    try:
        return x.astype(int).astype(str)
    except ValueError:
        return x.astype(str)

def annotate_read_depth(series_alt_count, series_ref_count):
    return series_alt_count.astype(int).add(series_ref_count.astype(int))

def get_idx_low_depth(series_depth, min_depth):
    return series_depth[series_depth.astype(int).le(int(min_depth))].index

def get_idx_coding_classifications(series_classification):
    coding_classifications = [
        'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site',
        'Frame_Shift_Ins', 'Frame_Shift_Del', 'In_Frame_Ins', 'In_Frame_Del']
    return series_classification[series_classification.isin(coding_classifications)].index

def main(inputs):
    df = standard_read(inputs[maf_handle], maf_column_map, low_memory=False)
    exac = standard_read(inputs[exac_handle], exac_column_map, low_memory=False)
    df.drop(df.columns[df.columns.str.lower().str.contains('exac')], axis=1, inplace=True)

    df.loc[:, LOW_DEPTH] = np.nan
    df.loc[:, CODING] = np.nan
    df.loc[:, WL] = np.nan
    df.loc[:, EXAC_COMMON] = 0.0
    idx_original = df.index

    df[DEPTH] = annotate_read_depth(df[MAPPED_ALT_COUNT], df[MAPPED_REF_COUNT])
    df.loc[:, LOW_DEPTH] = 0.0
    idx_read_depth = get_idx_low_depth(df[DEPTH], inputs[min_depth])

    df.loc[:, CODING] = 0.0
    idx_coding = get_idx_coding_classifications(df[MAPPED_VAR_CLASS])
    idx_noncoding = idx_original.difference(idx_coding)

    if inputs[whitelist_handle]:
        df.loc[:, WL] = 0.0

        whitelist = read(inputs_dict[whitelist_handle], header=-1).rename(columns=whitelist_column_map)
        df[whitelist_column_map[3]] = df[MAPPED_GENE].astype(str) + ':' + \
                                      df[MAPPED_AA].str.split('p.', expand=True).loc[:, 1].astype(str)
        df[whitelist_column_map[3]] = df[whitelist_column_map[3]].fillna('')
        idx_whitelist = df[df[whitelist_column_map[3]].isin(whitelist[whitelist_column_map[3]])].index
    else:
        idx_whitelist = pd.DataFrame().index

    merge_cols = [MAPPED_CHR, MAPPED_POS, MAPPED_REF, MAPPED_ALT]
    df = df.merge(exac, on=merge_cols, how='left')

    df.loc[:, populations] = df.loc[:, populations].fillna(0.0)
    idx_common_exac = df[(df.loc[:, populations].astype(float) > 10.0).sum(axis=1) != 0].index

    df.loc[idx_read_depth, LOW_DEPTH] = 1.0
    df.loc[idx_coding, CODING] = 1.0
    df.loc[idx_whitelist, WL] = 1.0
    df.loc[idx_common_exac, EXAC_COMMON] = 1.0

    idx_reject = idx_read_depth.union(idx_noncoding).union(idx_common_exac).difference(idx_whitelist)
    idx_common = idx_common_exac.difference(idx_whitelist)
    idx_pass = idx_original.difference(idx_reject)

    df[COMMON] = 0
    df.loc[idx_common, COMMON] = 1

    df.drop(whitelist_column_map[3], axis=1, inplace=True)

    df_pass = df.loc[idx_pass, :]
    df_reject = df.loc[idx_reject, :]

    df_pass.to_csv('pass.txt', sep='\t', index=False)
    df_reject.to_csv('reject.txt', sep='\t', index=False)
    df.to_csv('testoutput.txt', sep='\t')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--maf', help='MAF to annotate and filter', required=True)
    parser.add_argument('--exac', help='formatted exac', required=True)
    parser.add_argument('--whitelist', help='whitelist for somatic sites', default=False)
    parser.add_argument('--min_exac_ac', help='Minimum allele count across any population to filter', default=10)
    parser.add_argument('--filter_syn', help='Removes syn variants. True/False', action='store_true', default=False)
    parser.add_argument('--filter_read_depth', help='Filters based on specified read depth. Int.', default=0)
    args = parser.parse_args()

    inputs_dict = {
        maf_handle: args.maf,
        exac_handle: args.exac,
        whitelist_handle: args.whitelist,
        filter_syn: args.filter_syn,
        min_exac_ac: args.min_exac_ac,
        min_depth: args.filter_read_depth
    }

    main(inputs_dict)
