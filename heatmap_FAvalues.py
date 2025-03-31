import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Function to filter SNPs based on P value and OR cutoffs
def filter_snps(df, p_value_cutoff, or_cutoff):
    return df[(df['P'] <= p_value_cutoff) & (df['OR'] >= or_cutoff)]['SNP']

# Step 1: Prepare the SNP list
def prepare_snp_list(file_atypical, file_intermediate, file_typical, p_value_cutoff, or_cutoff):
    # Read files
    df_atypical = pd.read_csv(file_atypical, delim_whitespace=True)
    df_intermediate = pd.read_csv(file_intermediate, delim_whitespace=True)
    df_typical = pd.read_csv(file_typical, delim_whitespace=True)
    
    # Apply cutoffs and get SNP IDs
    snps_atypical = filter_snps(df_atypical, p_value_cutoff, or_cutoff)
    snps_intermediate = filter_snps(df_intermediate, p_value_cutoff, or_cutoff)
    snps_typical = filter_snps(df_typical, p_value_cutoff, or_cutoff)
    
    # Take union of SNP IDs
    union_snps = pd.Series(list(set(snps_atypical) | set(snps_intermediate) | set(snps_typical)))
    
    return union_snps

# Step 2: Map the F_A values
def map_fa_values(union_snps, file, subtype):
    # Read file
    df = pd.read_csv(file, delim_whitespace=True)
    
    # Select union SNPs and keep F_A values
    df_filtered = df[df['SNP'].isin(union_snps)][['SNP', 'F_A']]
    
    # Add missing SNPs with F_A value 0
    missing_snps = union_snps[~union_snps.isin(df_filtered['SNP'])]
    df_missing = pd.DataFrame({'SNP': missing_snps, 'F_A': 0})
    
    # Combine and add subtype column
    df_combined = pd.concat([df_filtered, df_missing])
    df_combined['Subtype'] = subtype
    
    return df_combined

# Step 3: Combine F_A values into desired format
def combine_fa_values(df_atypical, df_intermediate, df_typical, format_option='rows'):
    if format_option == 'columns':
        # Merge dataframes on SNP
        combined_df = pd.merge(df_atypical, df_intermediate, on='SNP', how='outer', suffixes=('_atypical', '_intermediate'))
        combined_df = pd.merge(combined_df, df_typical, on='SNP', how='outer')
        combined_df.rename(columns={'F_A': 'F_typical'}, inplace=True)
    else:
        # Concatenate dataframes
        combined_df = pd.concat([df_atypical, df_intermediate, df_typical])
    
    return combined_df

# Step 4: Generate the heatmap
def generate_heatmap(df, format_option='rows'):
    if format_option == 'columns':
        heatmap_data = df.set_index('SNP')
        sns.heatmap(heatmap_data, annot=True, cmap="RdPu")
        plt.title("F_A Values of SNPs Across Alzheimer's Disease Subtypes", fontsize = 16)
        plt.show()
    else:
        pivot_df = df.pivot(index='SNP', columns='Subtype', values='F_A')
        sns.heatmap(pivot_df, annot=True, cmap="RdPu")
        plt.title("F_A Values of SNPs Across Alzheimer's Disease Subtypes", fontsize = 16)
        plt.show()

# Main execution
if __name__ == "__main__":
    # File paths
    file_atypical = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_atypical_control_filtered_high_fisher.assoc.fisher'
    file_intermediate = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_intermediate_control_filtered_high_fisher.assoc.fisher'
    file_typical = 'AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_typical_control_filtered_high_fisher.assoc.fisher'
    
    # Cutoffs
    p_value_cutoff = 0.05
    or_cutoff = 1.5
    
    # Prepare SNP list
    union_snps = prepare_snp_list(file_atypical, file_intermediate, file_typical, p_value_cutoff, or_cutoff)
    
    # Map F_A values
    df_atypical = map_fa_values(union_snps, file_atypical, 'atypical')
    df_intermediate = map_fa_values(union_snps, file_intermediate, 'intermediate')
    df_typical = map_fa_values(union_snps, file_typical, 'typical')
    
    # Combine F_A values (option: 'columns' or 'rows')
    format_option = 'rows'  # Change to 'columns' if preferred
    combined_df = combine_fa_values(df_atypical, df_intermediate, df_typical, format_option=format_option)
    
    # Save the combined data
    combined_df.to_csv('combined_fa_values.csv', index=False)
    
    # Generate the heatmap
    generate_heatmap(combined_df, format_option=format_option)
