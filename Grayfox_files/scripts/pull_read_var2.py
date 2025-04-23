import argparse
import os
from statistics import variance
import pysam
from collections import defaultdict

def extract_unmapped_read_ids(bam_file, output_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(output_file, "w") as out:
        for read in bam:
            if read.is_unmapped:
                out.write(read.query_name + "\n")

def count_reads(bam_file):
    mapped = 0
    unmapped = 0
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue  # Skip non-primary alignments
            if read.is_unmapped:
                unmapped += 1
            else:
                mapped += 1
    return mapped, unmapped

def calculate_variance_across_assemblies(output_dir, summary_file, bam_files):
    samples = defaultdict(dict)
    mapped_read_counts = defaultdict(dict)
    unmapped_read_counts = defaultdict(dict)

    for filepath in bam_files:
        base = os.path.basename(filepath)
        parts = base.split(".")
        if len(parts) < 2:
            continue
        sample = parts[0]
        genome = parts[1]
        mapped_count, unmapped_count = count_reads(filepath)
        mapped_read_counts[sample][genome] = mapped_count
        unmapped_read_counts[sample][genome] = unmapped_count

    for filepath in os.listdir(output_dir):
        if not filepath.endswith("_unmapped.txt"):
            continue
        base = os.path.basename(filepath)
        parts = base.split(".")
        sample = parts[0]
        genome = parts[1]
        with open(os.path.join(output_dir, filepath)) as f:
            reads = set(line.strip() for line in f if line.strip())
        samples[sample][genome] = reads

    with open(summary_file, "w") as summary:
        summary.write("Sample\tMappedReads(gf:cfam3:cfam4:arctic4)\tUnmappedReads(gf:cfam3:cfam4:arctic4)\tVariance\tUnmappedIntersection\tTotalMapped\tJaccardIndices(gf:cfam3:cfam4:arctic4)\n")

        for sample, genome_reads in samples.items():
            read_counts = {asm: len(reads) for asm, reads in genome_reads.items()}
            unmapped_intersection = set.intersection(*genome_reads.values()) if len(genome_reads) == 4 else set()
            total_mapped = sum(mapped_read_counts[sample].values())

            summary.write(f"{sample}\t")
            summary.write(", ".join([f"{asm}:{mapped_read_counts[sample].get(asm, 0)}" for asm in ["gf", "cfam3", "cfam4", "arctic4"]]))
            summary.write("\t")
            summary.write(", ".join([f"{asm}:{unmapped_read_counts[sample].get(asm, 0)}" for asm in ["gf", "cfam3", "cfam4", "arctic4"]]))

            if len(read_counts) == 4:
                var = variance(read_counts.values())

                jaccard_vals = {}
                assemblies = list(genome_reads.keys())
                for i in range(len(assemblies)):
                    for j in range(i + 1, len(assemblies)):
                        a, b = assemblies[i], assemblies[j]
                        set_a, set_b = genome_reads[a], genome_reads[b]
                        intersection = len(set_a & set_b)
                        union = len(set_a | set_b)
                        jaccard = intersection / union if union > 0 else 0
                        jaccard_vals[f"{a}-{b}"] = jaccard

                summary.write(f"\t{var:.2f}\t{len(unmapped_intersection)}\t{total_mapped}\t")
                summary.write(", ".join([f"{k}:{v:.2f}" for k, v in jaccard_vals.items()]))
            else:
                summary.write("\tNA\tNA\tNA")

            summary.write("\n")
            print(f"Processed {sample}: variance and similarities written to summary.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract unmapped read IDs and calculate variance across assemblies.")
    parser.add_argument("-l", "--list_file", required=True, help="Path to a text file containing a list of BAM files")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output directory")
    parser.add_argument("--calc_variance", action="store_true", help="Calculate variance across genome assemblies")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.list_file, "r") as f:
        bam_files = [line.strip() for line in f if line.strip()]

    for bam_file in bam_files:
        output_file = os.path.join(args.output_dir, os.path.basename(bam_file).replace(".bam", "_unmapped.txt"))
        if not os.path.exists(output_file):
            try:
                extract_unmapped_read_ids(bam_file, output_file)
                print(f"Unmapped read IDs from {bam_file} have been saved to {output_file}")
            except Exception as e:
                print(f"Error processing {bam_file}: {e}")
        else:
            print(f"Output file {output_file} already exists. Skipping extraction.")

    if args.calc_variance:
        summary_file = os.path.join(args.output_dir, "assembly_comparison_summary.tsv")
        calculate_variance_across_assemblies(args.output_dir, summary_file, bam_files)
