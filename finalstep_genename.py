import pandas as pd

# a lot of similar code to heatmap_Pvalue
def filter_snps(df, p_value_cutoff, or_cutoff):
    return df[(df['P'] <= p_value_cutoff) & (df['OR'] >= or_cutoff)]['SNP']

def prepare_snp_list(file_atypical, file_intermediate, file_typical, p_value_cutoff, or_cutoff):
    df_atypical = pd.read_csv(file_atypical, sep=r'\s+')
    df_intermediate = pd.read_csv(file_intermediate, sep=r'\s+')
    df_typical = pd.read_csv(file_typical, sep=r'\s+')
    
    snps_atypical = filter_snps(df_atypical, p_value_cutoff, or_cutoff)
    snps_intermediate = filter_snps(df_intermediate, p_value_cutoff, or_cutoff)
    snps_typical = filter_snps(df_typical, p_value_cutoff, or_cutoff)
    
    union_snps = pd.Series(list(set(snps_atypical) | set(snps_intermediate) | set(snps_typical)))
    
    return union_snps

def map_p_values(union_snps, file):
    df = pd.read_csv(file, sep=r'\s+')
    
    df_filtered = df[df['SNP'].isin(union_snps)][['SNP', 'P']]
    
    missing_snps = union_snps[~union_snps.isin(df_filtered['SNP'])]
    df_missing = pd.DataFrame({'SNP': missing_snps, 'P': 1})
    
    df_combined = pd.concat([df_filtered, df_missing])
    
    return df_combined

def combine_p_values(df_atypical, df_intermediate, df_typical):
    combined_df = pd.concat([df_atypical, df_intermediate, df_typical]).drop_duplicates(subset=['SNP'])
    return combined_df

def map_snp_to_gene(heatmap_data_file, gene_mapping_file):
    heatmap_df = pd.read_csv(heatmap_data_file)
    
    # read gene mapping data
    gene_mapping_df = pd.read_csv(gene_mapping_file, sep=r'\s+', header=None, names=['gene_name', 'SNP'])
    
    # merge heatmap data with gene mapping data
    merged_df = pd.merge(heatmap_df, gene_mapping_df, on='SNP', how='left')
    
    # group by gene name
    grouped_df = merged_df.groupby('gene_name').sum().reset_index()
    
    return grouped_df

if __name__ == "__main__":
    file_atypical = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_atypical_control_filtered_high_fisher.assoc.fisher'
    file_intermediate = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_intermediate_control_filtered_high_fisher.assoc.fisher'
    file_typical = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_typical_control_filtered_high_fisher.assoc.fisher'
    gene_mapping_file = 'ampadwgs_coding_high_rs.SETID'
    
    p_value_cutoff = 0.05
    or_cutoff = 1.5
    
    union_snps = prepare_snp_list(file_atypical, file_intermediate, file_typical, p_value_cutoff, or_cutoff)
    
    print("Unique SNPs (rsIDs):")
    print(union_snps)
    
    # Map P values
    df_atypical = map_p_values(union_snps, file_atypical)
    df_intermediate = map_p_values(union_snps, file_intermediate)
    df_typical = map_p_values(union_snps, file_typical)
    
    combined_df = combine_p_values(df_atypical, df_intermediate, df_typical)
    
    combined_df.to_csv('combined_p_values.csv', index=False)
    
    # map SNPs to gene names and group by gene name
    grouped_gene_df = map_snp_to_gene('combined_p_values.csv', gene_mapping_file)
    
    # Print grouped data by gene name
    print("Grouped Data by Gene Name:")
    print(grouped_gene_df)
    
    # save the grouped data
    grouped_gene_df.to_csv('grouped_gene_data.csv', index=False)
