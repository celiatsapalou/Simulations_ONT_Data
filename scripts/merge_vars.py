import argparse

parser = argparse.ArgumentParser(description='Merge nearby VCF variants.')
parser.add_argument('input_vcf', type=str, help='Input VCF file.')
parser.add_argument('output_vcf', type=str, help='Output VCF file.')

args = parser.parse_args()

input_vcf = args.input_vcf
output_vcf = args.output_vcf


def write_variant(fout, variant):
    if variant:
        if len(variant["alts"]) > 1:
            alt = ','.join(variant["alts"])
        else:
            alt = variant["alts"][0]
        fields = [variant["chrom"], variant["pos"], ".", variant["ref"], alt]
        merged_line = '\t'.join(map(str, fields)) + '\t' + '\t'.join(variant["rest"])
        fout.write(merged_line + '\n')

with open(input_vcf, 'r') as fin, open(output_vcf, 'w') as fout:
    current_variant = None
    for line in fin:
        # Copy header lines directly to the output
        if line.startswith("#"):
            fout.write(line)
            continue

        fields = line.strip().split("\t")
        chrom, pos, _, ref, alt, *rest = fields
        pos = int(pos)

        # If there's no current variant being processed
        if not current_variant:
            current_variant = {"chrom": chrom, "pos": pos, "ref": ref, "alts": [alt], "rest": rest}
        else:
            # Check if the current line's variant is within 10 base pairs of the current variant
            if chrom == current_variant["chrom"] and abs(pos - current_variant["pos"]) <= 20:
                # Add the alt of the current line to the current variant's alts
                current_variant["alts"].append(alt)
            else:
                # Write out the current_variant
                write_variant(fout, current_variant)

                # Start processing the next variant
                current_variant = {"chrom": chrom, "pos": pos, "ref": ref, "alts": [alt], "rest": rest}

    # Handle the last variant
    write_variant(fout, current_variant)

