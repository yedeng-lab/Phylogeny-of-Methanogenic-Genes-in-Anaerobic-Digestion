import argparse
from Bio import SeqIO

def split_fasta_by_gene(input_fasta, output_dir):
        gene_dict = {}
            
        for record in SeqIO.parse(input_fasta, "fasta"):
            parts = record.id.split('_')
            if len(parts) > 1:
                gene_label = parts[0]
                if gene_label not in gene_dict:
                    gene_dict[gene_label] = []
                gene_dict[gene_label].append(record)
        
        for gene_label, records in gene_dict.items():
            output_path = f"{output_dir}/{gene_label}.fasta"
            #SeqIO.write(records, output_path, "fasta")
            with open(output_path, 'w') as output_file:
                for record in records:
                    output_file.write(f">{record.id}\n")
                    output_file.write(f"{str(record.seq)}\n")
            print(f"Written {len(records)} sequences to {output_path}")

if __name__ == "__main__":
        parser = argparse.ArgumentParser(description="Split FASTA file by gene labels.")
        parser.add_argument("input_fasta", help="Input FASTA file path.")
        parser.add_argument("output_dir", help="Output directory for split FASTA files.")
                        
        args = parser.parse_args()
        split_fasta_by_gene(args.input_fasta, args.output_dir)

