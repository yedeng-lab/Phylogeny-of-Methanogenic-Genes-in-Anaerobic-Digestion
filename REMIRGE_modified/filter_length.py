import os
import numpy as np
from Bio import SeqIO
import argparse
from scipy import stats

def write_single_line_fasta(sequences, output_file):
    with open(output_file, 'w') as out_f:
        for record in sequences:
            out_f.write(f">{record.id}\n{record.seq}\n")

def process_gene_sequences(input_dir, output_dir, n):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    summary_stats = []

    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta"):
            gene_name = filename[:-6]  # Remove '.fasta' from the filename to get the gene name
            input_path = os.path.join(input_dir, filename)
            sequences = list(SeqIO.parse(input_path, "fasta"))
            lengths = [len(seq) for seq in sequences]

            # Calculate statistics
            mean_length = np.mean(lengths)
            max_length = np.max(lengths)
            min_length = np.min(lengths)
            median_length = np.median(lengths)
            upper_quartile = np.percentile(lengths, 75)
            lower_quartile = np.percentile(lengths, 25)
            above_mean_count = sum(length > mean_length for length in lengths)
            above_median_count = sum(length > median_length for length in lengths)
            above_upper_quartile_count = sum(length > upper_quartile for length in lengths)
            above_lower_quartile_count = sum(length > lower_quartile for length in lengths)
            between_lower_upper_count = sum(length > lower_quartile and length < upper_quartile for length in lengths)

            # Filter sequences above mean length
            if max_length / min_length > n:
                filtered_sequences = [seq for seq in sequences if len(seq) >= mean_length]
                filtered_count = len(filtered_sequences)
            else:
                filtered_sequences = [seq for seq in sequences ]
                filtered_count = len(filtered_sequences)

            # Output filtered sequences to new file
            output_path = os.path.join(output_dir, f"{gene_name}.fasta")
            write_single_line_fasta(filtered_sequences, output_path)

            # Append statistics to summary
            summary_stats.append({
                "Gene": gene_name,
                "Total_Sequences": len(sequences),
                "Filtered_Sequences": filtered_count,
                "Proportion_Filtered": filtered_count / len(sequences) if len(sequences) > 0 else 0,
                "Mean_Length": mean_length,
                "Filtered_Mean_Length": np.mean([len(seq) for seq in filtered_sequences]) if filtered_sequences else 0,
                "Max_Length": max_length,
                "Min_Length": min_length,
                "Median_Length": median_length,
                "Upper_Quartile": upper_quartile,
                "Lower_Quartile": lower_quartile,
                "Above_Mean_Count": above_mean_count,
                "Above_Median_Count": above_median_count,
                "Above_Upper_Quartile_Count": above_upper_quartile_count,
                "Above_Lower_Quartile_Count": above_lower_quartile_count,
                "Between_Lower_and_Upper_Quartile_Count": between_lower_upper_count,
            })

    # Write summary statistics to a file
    summary_file = os.path.join(output_dir, "summary_statistics.txt")
    with open(summary_file, 'w') as summary_f:
        summary_f.write("Gene\tTotal_Sequences\tFiltered_Sequences\tProportion_Filtered\tMean_Length\tFiltered_Mean_Length\tMax_Length\tMin_Length\tMedian_Length\tUpper_Quartile\tLower_Quartile\tAbove_Mean_Count\tAbove_Median_Count\tAbove_Upper_Quartile_Count\tAbove_Lower_Quartile_Count\tBetween_Lower_and_Upper_Quartile_Count\n")
        for stats in summary_stats:
            summary_f.write(f"{stats['Gene']}\t{stats['Total_Sequences']}\t{stats['Filtered_Sequences']}\t{stats['Proportion_Filtered']:.2f}\t{stats['Mean_Length']:.2f}\t{stats['Filtered_Mean_Length']:.2f}\t{stats['Max_Length']}\t{stats['Min_Length']}\t{stats['Median_Length']}\t{stats['Upper_Quartile']}\t{stats['Lower_Quartile']}\t{stats['Above_Mean_Count']}\t{stats['Above_Median_Count']}\t{stats['Above_Upper_Quartile_Count']}\t{stats['Above_Lower_Quartile_Count']}\t{stats['Between_Lower_and_Upper_Quartile_Count']}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process gene sequences in FASTA files.')
    parser.add_argument('input_dir', type=str, help='Directory containing input FASTA files')
    parser.add_argument('output_dir', type=str, help='Directory to save output files and statistics')
    parser.add_argument('n', type=float, help='Sequence length variation threshold to filter nucleotide sequences i.e. max_length/min_length')

    args = parser.parse_args()
    process_gene_sequences(args.input_dir, args.output_dir, args.n)
