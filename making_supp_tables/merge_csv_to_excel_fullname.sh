#Written by ChatGPT 4.0 on 10/31/24

#!/bin/bash

# Requires xlsxwriter library in Python
# Install by running: pip install XlsxWriter

# Name of the output Excel file
output_file="combined_output.xlsx"

# Create a temporary Python script to merge the CSVs into an Excel workbook
echo 'import sys
import glob
import pandas as pd
from xlsxwriter import Workbook

# Create a workbook and add a worksheet for each CSV
workbook = Workbook(sys.argv[1])
for csv_file in glob.glob("*.csv"):
    # Use the full filename (without .csv extension) as the tab name
    tab_name = csv_file.rsplit(".", 1)[0]
    df = pd.read_csv(csv_file)
    worksheet = workbook.add_worksheet(tab_name[:31])  # Limit tab names to 31 chars
    for i, col_name in enumerate(df.columns):
        worksheet.write(0, i, col_name)
        for j, value in enumerate(df[col_name]):
            worksheet.write(j + 1, i, value)
workbook.close()' > merge_csv_to_excel.py

# Run the Python script with the desired output file
python3 merge_csv_to_excel.py "$output_file"

# Clean up the temporary Python script
rm merge_csv_to_excel.py
