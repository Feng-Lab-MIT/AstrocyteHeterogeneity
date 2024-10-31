#From ChatGPT 4.0 on 10/31/24

#!/bin/bash

# Requires xlsxwriter library in Python
# Install by running: pip install XlsxWriter

# Name of the output Excel file
output_file="combined_output.xlsx"

# Create a temporary Python script to merge the text files into an Excel workbook
echo 'import sys
import glob
import pandas as pd
from xlsxwriter import Workbook
import re
import numpy as np

# Create a workbook and set nan_inf_to_errors=True to handle NaN/Inf values
workbook = Workbook(sys.argv[1], {"nan_inf_to_errors": True})
for txt_file in glob.glob("*.txt"):
    # Extract the tab name using the pattern "RM one-way ANOVA of <name>.txt"
    match = re.search(r"RM one-way ANOVA of (.+?)\.txt", txt_file)
    if match:
        tab_name = match.group(1)[:31]  # Limit tab names to 31 characters
    else:
        tab_name = txt_file.rsplit(".", 1)[0][:31]  # Use filename if pattern not found

    # Read the text file into a DataFrame and replace NaN/Inf with empty strings
    df = pd.read_csv(txt_file, delimiter="\t")  # Assuming tab-delimited text files
    df = df.replace([np.inf, -np.inf, np.nan], "")  # Replace NaN/Inf with empty strings

    # Write the DataFrame to the worksheet
    worksheet = workbook.add_worksheet(tab_name)
    for i, col_name in enumerate(df.columns):
        worksheet.write(0, i, col_name)
        for j, value in enumerate(df[col_name]):
            worksheet.write(j + 1, i, value)
workbook.close()' > merge_txt_to_excel.py

# Run the Python script with the desired output file
python3 merge_txt_to_excel.py "$output_file"

# Clean up the temporary Python script
rm merge_txt_to_excel.py

