import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# here filter SNPs based on P value and OR cutoffs
def filter_snps(df, p_value_cutoff, or_cutoff):
    return df[(df['P'] <= p_value_cutoff) & (df['OR'] >= or_cutoff)]['SNP']

# read the files
def prepare_snp_list(file_atypical, file_intermediate, file_typical, p_value_cutoff, or_cutoff):
    df_atypical = pd.read_csv(file_atypical, sep=r'\s+')
    df_intermediate = pd.read_csv(file_intermediate, sep=r'\s+')
    df_typical = pd.read_csv(file_typical, sep=r'\s+')
    
    # use the cutoffs and get SNP IDs
    snps_atypical = filter_snps(df_atypical, p_value_cutoff, or_cutoff)
    snps_intermediate = filter_snps(df_intermediate, p_value_cutoff, or_cutoff)
    snps_typical = filter_snps(df_typical, p_value_cutoff, or_cutoff)
    
    # get the union of SNP IDs
    union_snps = pd.Series(list(set(snps_atypical) | set(snps_intermediate) | set(snps_typical)))
    
    return union_snps

# map the P values
def map_p_values(union_snps, file, subtype):
    df = pd.read_csv(file, sep=r'\s+')
    
    # select rows from the input dataframe that is there in the union_snps 
    # keep only the SNPs and keep P values
    df_filtered = df[df['SNP'].isin(union_snps)][['SNP', 'P']]
    
    # if there are missing SNPs, make P value 1 because theres no association
    missing_snps = union_snps[~union_snps.isin(df_filtered['SNP'])]
    df_missing = pd.DataFrame({'SNP': missing_snps, 'P': 1})
    
    # combine and add subtype column
    df_combined = pd.concat([df_filtered, df_missing])
    df_combined['Subtype'] = subtype
    
    return df_combined

# make P values into different format
def combine_p_values(df_atypical, df_intermediate, df_typical, format_option='rows'):
    if format_option == 'columns':
        # merge dataframes
        combined_df = pd.merge(df_atypical, df_intermediate, on='SNP', how='outer', suffixes=('_atypical', '_intermediate'))
        combined_df = pd.merge(combined_df, df_typical, on='SNP', how='outer')
        combined_df.rename(columns={'P': 'P_typical'}, inplace=True)
    else:
        # combime multiple dataframes into one dataframe
        combined_df = pd.concat([df_atypical, df_intermediate, df_typical])
    return combined_df

# heatmap work
def generate_heatmap(df, format_option='rows', cmap='RdPu'):
    if format_option == 'columns':
        heatmap_data = df.set_index('SNP')
        sns.heatmap(heatmap_data, annot=True)
        plt.title("P Values of SNPs Across Alzheimer's Disease Subtypes", fontsize = 16)

        plt.show()
    else:
        pivot_df = df.pivot(index='SNP', columns='Subtype', values='P')
        sns.heatmap(pivot_df, annot=True, cmap="RdPu")
        plt.title("P Values of SNPs Across Alzheimer's Disease Subtypes", fontsize = 16)
        plt.show()

if __name__ == "__main__":
    file_atypical = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_atypical_control_filtered_high_fisher.assoc.fisher'
    file_intermediate = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_intermediate_control_filtered_high_fisher.assoc.fisher'
    file_typical = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_typical_control_filtered_high_fisher.assoc.fisher'
    
    # cutoffs
    p_value_cutoff = 0.05
    or_cutoff = 1.5
    
    # prep SNP list
    union_snps = prepare_snp_list(file_atypical, file_intermediate, file_typical, p_value_cutoff, or_cutoff)
   
    # mapping the P values
    df_atypical = map_p_values(union_snps, file_atypical, 'atypical')
    df_intermediate = map_p_values(union_snps, file_intermediate, 'intermediate')
    df_typical = map_p_values(union_snps, file_typical, 'typical')
    
    # combining P values
    format_option = 'rows'  # Change to 'columns' if preferred
    combined_df = combine_p_values(df_atypical, df_intermediate, df_typical, format_option=format_option)
    
    # save the heatmap data
    combined_df.to_csv('combined_p_values.csv', index=False)
    
    # make the heatmap
    generate_heatmap(combined_df, format_option=format_option, cmap='RdPu')

