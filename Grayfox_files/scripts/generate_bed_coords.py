import pysam
import argparse
import os
from collections import defaultdict

def load_unmapped_ids(file_path):
    with open(file_path, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def find_conditionally_unmapped_reads(gf_file, other_files):
    gf_unmapped = load_unmapped_ids(gf_file)
    result = {}
    for asm, filepath in other_files.items():
        other_unmapped = load_unmapped_ids(filepath)
        result[asm] = other_unmapped - gf_unmapped
    return result

def extract_mapped_coordinates(bam_file, read_ids, bed_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(bed_file, "w") as out:
        for read in bam:
            if read.query_name in read_ids and not read.is_unmapped:
                chrom = read.reference_name
                start = read.reference_start
                end = read.reference_end
                name = read.query_name
                out.write(f"{chrom}\t{start}\t{end}\t{name}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract coordinates of conditionally unmapped reads for multiple samples.")
    parser.add_argument("-l", "--list_file", required=True, help="Path to a text file containing a list of BAM files")
    parser.add_argument("-u", "--unmapped_dir", required=True, help="Directory containing unmapped read ID files")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory to save output BED files")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.list_file, 'r') as f:
        bam_files = [line.strip() for line in f if line.strip()]

    for bam_file in bam_files:
        bam_basename = os.path.basename(bam_file)
        sample = bam_basename.replace(".gf.sorted.bam", "")

        unmapped_files = {
            'gf': os.path.join(args.unmapped_dir, f"{sample}.gf.sorted_unmapped.txt"),
            'cfam3': os.path.join(args.unmapped_dir, f"{sample}.cfam3.sorted_unmapped.txt"),
            'cfam4': os.path.join(args.unmapped_dir, f"{sample}.cfam4.sorted_unmapped.txt"),
            'arctic4': os.path.join(args.unmapped_dir, f"{sample}.arctic4.sorted_unmapped.txt")
        }

        if not all(os.path.exists(f) for f in unmapped_files.values()):
            print(f"Skipping {sample}, not all unmapped files found.")
            continue

        other_assemblies = {k: v for k, v in unmapped_files.items() if k != 'gf'}
        selected_reads = find_conditionally_unmapped_reads(unmapped_files['gf'], other_assemblies)

        for asm, ids in selected_reads.items():
            bed_path = os.path.join(args.output_dir, f"{sample}_{asm}_vs_gf.bed")
            extract_mapped_coordinates(bam_file, ids, bed_path)
            print(f"Extracted coordinates for {sample} ({asm} not in gf) to: {bed_path}")

