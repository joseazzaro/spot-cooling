# -*- coding: utf-8 -*-
import pdfplumber
import json

pdf_files = [
    "ASHRAE-D-HO-2667.pdf",
    "ASHRAE-D-HO-2668.pdf"
]

for pdf_file in pdf_files:
    print(f"\n{'='*80}")
    print(f"Extracting from: {pdf_file}")
    print(f"{'='*80}\n")
    
    try:
        with pdfplumber.open(pdf_file) as pdf:
            print(f"Total pages: {len(pdf.pages)}\n")
            
            # Extract text from all pages
            for i, page in enumerate(pdf.pages):
                print(f"\n--- PAGE {i+1} ---")
                text = page.extract_text()
                print(text)
                
                # Try to extract tables
                tables = page.extract_tables()
                if tables:
                    print(f"\nTables found on page {i+1}:")
                    for j, table in enumerate(tables):
                        print(f"\nTable {j+1}:")
                        for row in table:
                            print(row)
                
    except Exception as e:
        print(f"Error reading {pdf_file}: {e}")
