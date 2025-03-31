import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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

def map_frequency(union_snps, file, subtype):
    df = pd.read_csv(file, sep=r'\s+')
    
    # select rows from the input dataframe that is there in the union_snps 
    # keep only the SNPs and keep P values
    df_filtered = df[df['SNP'].isin(union_snps)][['SNP', 'F_A']]
    
    # if there are missing SNPs, make P value 1 because theres no association
    missing_snps = union_snps[~union_snps.isin(df_filtered['SNP'])]
    df_missing = pd.DataFrame({'SNP': missing_snps, 'F_A': 0})
    
    # combine and add subtype column
    df_combined = pd.concat([df_filtered, df_missing])
    df_combined['Subtype'] = subtype
    
    return df_combined




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
    df_atypical = map_frequency(union_snps, file_atypical, 'atypical')
    df_intermediate = map_frequency(union_snps, file_intermediate, 'intermediate')
    df_typical = map_frequency(union_snps, file_typical, 'typical')
    
    combined_df = pd.concat([df_atypical, df_intermediate, df_typical])
    
    
    gene_mapping_df = pd.read_csv(gene_mapping_file, sep=r'\s+', header=None, names=['gene_name', 'SNP'])
    
    # merge heatmap data with gene mapping data
    merged_df = pd.merge(combined_df, gene_mapping_df, on='SNP', how='left')
    
    # group by gene name
    grouped_df = merged_df.groupby(['Subtype','gene_name']).sum().reset_index()

 
    pivot_df = grouped_df.pivot(index='gene_name', columns='Subtype', values='F_A')
    plt.figure()
    sns.heatmap(pivot_df, annot=True, cmap="RdPu")
    plt.title("Allele Frequencies by Gene Name Across Alzheimer's Disease Subtypes", fontsize=16)
    plt.show()
