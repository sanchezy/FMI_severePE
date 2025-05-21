#!/usr/bin/env python

#This script aims to:
#Concatenate in a single file all the outputs from scrublet. The file will be loaded into Seurat as metadata. Add a column with the name of the sample it comes from.
# Run as: python scrublet-oneoutput-file-5thRun.py
import pandas as pd

# Give a list of samples ID (these are the directory names from cellranger). Note the order must be the same as the order of aggregate from cellranger.
# In this list, eliminate the first element
samplesID = [
"F1758LJ-CAM-sn", 
"F1676VQ-CAM-sn",
"F1678CM-CAM-sn",
"F1682RH-CAM-sn", 
"F1686GS-CAM-sn",
"F1668RK-CAM-sn"]

# This are local directories
input_dir = '../FMI-pilot/scrublet-using-CellBender-input-5thRun/'
output_dir = '../FMI-pilot/scrublet-using-CellBender-input-5thRun/'

# Initialise a variable that will add a string to the barcodes.
i = 1 

#open the first one
df1 = pd.read_csv(input_dir + 'F1758LJ-CAM-sn_scrublet_output_table.csv')
df1['sample'] = "F1758LJ-CAM-sn"

# Iterate over the number of libraries and concatenate the results
for ID in samplesID:
    
    df2 = pd.read_csv(input_dir + ID + '_scrublet_output_table.csv')
    df2['sample'] = ID
    # Get rid off the last string. It is always 1.
    df2.barcodes = df2['barcodes'].str[:-1] 
    # Add a string in the barcodes, this will allow to avoid repetions when concatenating textfiles from different libraries
    df2.barcodes = df2.barcodes + str(i)
    final = pd.concat([df1,df2],ignore_index=True)
    df1 = final
    i = i+1 
    print('Number of cells:')
    print(final.shape[0])
    print('Done for sample ' + ID)

print('Number of cells:')
print(final.shape[0])

# Save all the file into a final table
final.to_csv(output_dir + 'scrublet_output_0.25-20230914.csv', index=False)


# Print number of doublets detected
n_doublets = final.groupby('predicted_doublet').count()
print(n_doublets)
print('Done')