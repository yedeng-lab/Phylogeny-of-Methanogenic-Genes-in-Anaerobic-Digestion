import csv
import os
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description="Split a GSC abundance file by gene names.")
parser.add_argument('--input_file', type=str, help='Path to the input aundance file.', required=True)
parser.add_argument('--output_dir', type=str, help='Directory to save the output abundance files.', required=True)

args = parser.parse_args()

gene_data = defaultdict(list)
with open(args.input_file, 'r', newline='', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    header = next(reader)
    for row in reader:
        gene_name = row[0].split('_')[0]
        gene_data[gene_name].append(row)

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

for gene_name, rows in gene_data.items():
    output_file = os.path.join(args.output_dir, f'{gene_name}.txt')
    with open(output_file, 'w', newline='') as txtfile:
        writer = csv.writer(txtfile, delimiter='\t')
        writer.writerow(header)
        writer.writerows(rows)

print(f"Finished，generated abundance tables for {len(gene_data)} functional genes in {args.output_dir}")

