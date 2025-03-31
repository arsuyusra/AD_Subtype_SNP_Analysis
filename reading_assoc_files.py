id_list=[] #list to keep selected rsIDs

with open("AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_atypical_control_filtered_high_fisher.assoc.fisher", "r") as myfile:
    header = myfile.readline() #skip the first line
    print(header)
    for line in myfile: #read each line
        data=line.split()
        chr=data[0]
        rsID=data[1]
        if data[4]!="NA": #if value doesn't exist, assign 0
          ad_freq=float(data[4]) #frequency in AD samples
        else:
          ad_freq=0
        if data[5]!="NA":
          ctrl_freq=float(data[5]) #frequency in controls
        else:
          ctrl_freq=0
        pvalue=float(data[7])
        if pvalue<0.05: #check if difference between AD and control samples is statistically significant 
            id_list.append(rsID)

#write the list into a file
with open('atypical_control_rsID.txt', 'w') as f:
    for id in id_list:
        f.write(f"{id}\n")

id_list1=[] #list to keep selected rsIDs

with open("AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_typical_control_filtered_high_fisher.assoc.fisher", "r") as myfile:
    header = myfile.readline() #skip the first line
    print(header)
    for line in myfile: #read each line
        data=line.split()
        chr=data[0]
        rsID=data[1]
        if data[4]!="NA": #if value doesn't exist, assign 0
          ad_freq=float(data[4]) #frequency in AD samples
        else:
          ad_freq=0
        if data[5]!="NA":
          ctrl_freq=float(data[5]) #frequency in controls
        else:
          ctrl_freq=0
        pvalue=float(data[7])
        if pvalue<0.05: #check if difference between AD and control samples is statistically significant 
            id_list1.append(rsID)

#write the list into a file
with open('typical_control_rsID.txt', 'w') as f:
    for id in id_list1:
        f.write(f"{id}\n")

id_list2=[] #list to keep selected rsIDs

with open("AMPADWGS_WGS_all_rsIDs_dbSNP153_020520_intermediate_control_filtered_high_fisher.assoc.fisher", "r") as myfile:
    header = myfile.readline() #skip the first line
    print(header)
    for line in myfile: #read each line
        data=line.split()
        chr=data[0]
        rsID=data[1]
        if data[4]!="NA": #if value doesn't exist, assign 0
          ad_freq=float(data[4]) #frequency in AD samples
        else:
          ad_freq=0
        if data[5]!="NA":
          ctrl_freq=float(data[5]) #frequency in controls
        else:
          ctrl_freq=0
        pvalue=float(data[7])
        if pvalue<0.05: #check if difference between AD and control samples is statistically significant 
            id_list2.append(rsID)

#write the list into a file
with open('intermediate_control_rsID.txt', 'w') as f: #w means write
    for id in id_list2:
        f.write(f"{id}\n")

text_file = open("typical_control_rsID.txt", "r") #r means read
ListToSort = text_file.read().splitlines()
text_file.close()

text_file1 = open("atypical_control_rsID.txt", "r")
ListToSort1 = text_file1.read().splitlines()
text_file1.close()

text_file2 = open("intermediate_control_rsID.txt", "r")
ListToSort2 = text_file2.read().splitlines()
text_file2.close()

print(ListToSort)
print(ListToSort1)
print(ListToSort2)

intersection_set = set(ListToSort) & set(ListToSort1) & set(ListToSort2) #sets cant have any duplicates
#wanted unique identifiers from the list 

print(intersection_set)

with open('intersection_set.txt', 'w') as file:
    for item in intersection_set:
        file.write(f"{item}\n")


atypical_specific = set(ListToSort1) - set(ListToSort) - set(ListToSort2)
intermediate_specific = set(ListToSort2) - set(ListToSort) - set(ListToSort1)
typical_specific = set(ListToSort) - set(ListToSort1) - set(ListToSort2)

# Write specific rsIDs to files
with open('atypical_specific_rsID.txt', 'w') as file:
    for item in atypical_specific:
        file.write(f"{item}\n")

with open('intermediate_specific_rsID.txt', 'w') as file:
    for item in intermediate_specific:
        file.write(f"{item}\n")

with open('typical_specific_rsID.txt', 'w') as file:
    for item in typical_specific:
        file.write(f"{item}\n")

import pandas as pd
import os

output_dir = 'Documents'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    

# Lists of rsIDs
atypical_specific = ['rs375386696', 'rs3131044', 'rs2229519', 'rs61737148', 'rs7447815', 'rs41286592', 'rs17602729', 'rs73345458', 'rs10529064', 'rs2758628', 'rs2857667', 'rs8192646', 'rs1799853', 'rs1513629', 'rs1799990', 'rs11681642', 'rs35778930', 'rs3892097', 'rs201093021', 'rs12409540']
intermediate_specific = ['rs762557740', 'rs13894', 'rs1849674', 'rs200978296', 'rs41268500', 'rs35068483', 'rs73299503', 'rs1638213', 'rs74830030', 'rs200265097', 'rs4639077', 'rs801175', 'rs8050680', 'rs77865995', 'rs73537037', 'rs3215610', 'rs6550714', 'rs17492695', 'rs3117232', 'rs112215666']
typical_specific = ['rs12692392', 'rs62538278', 'rs885985', 'rs61743546', 'rs1475768', 'rs9886752', 'rs142805492', 'rs181983734']

max_length = max(len(atypical_specific), len(intermediate_specific), len(typical_specific))

# even out the length
atypical_specific += [''] * (max_length - len(atypical_specific))
intermediate_specific += [''] * (max_length - len(intermediate_specific))
typical_specific += [''] * (max_length - len(typical_specific))

# create the dataframe
data = {'Atypical-specific rsIDs': atypical_specific,
    'Intermediate-specific rsIDs': intermediate_specific,
    'Typical-specific rsIDs': typical_specific
}

df = pd.DataFrame(data)

# this is gonna save it to a CSV file
output_file = os.path.join(output_dir, 'specific_rsIDs_table.csv')
df.to_csv(output_file, index=False)

# this is the DataFrame
print(df)
