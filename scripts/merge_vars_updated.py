import argparse

parser = argparse.ArgumentParser(description='Merge nearby VCF variants.')
parser.add_argument('input_vcf', type=str, help='Input VCF file.')
parser.add_argument('output_vcf', type=str, help='Output VCF file.')

args = parser.parse_args()

input_vcf = args.input_vcf
output_vcf = args.output_vcf


def merge_info_fields(info_list):
    # For simplicity, this function will just concatenate the INFO fields
    # You can add more advanced merging logic if necessary
    return ";".join(info_list)


def write_variant(fout, variant):
    if variant:
        # Ensure unique alt values and merge INFO fields
        unique_alts = list(set(variant["alts"]))
        merged_info = merge_info_fields([rest[6] for rest in variant["rest"]])

        if len(unique_alts) > 1:
            alt = ','.join(unique_alts)
        else:
            alt = unique_alts[0]

        rest_of_the_fields = variant["rest"][0][:6] + [merged_info] + variant["rest"][0][7:]
        fields = [variant["chrom"], variant["pos"], ".", variant["ref"], alt] + rest_of_the_fields

        fout.write('\t'.join(map(str, fields)) + '\n')


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
            current_variant = {"chrom": chrom, "pos": pos, "ref": ref, "alts": [alt], "rest": [rest]}
        else:
            # Check if the current line's variant is within 20 base pairs of the current variant
            if chrom == current_variant["chrom"] and abs(pos - current_variant["pos"]) <= 20:
                # Add the alt of the current line to the current variant's alts and the rest of the fields to the rest
                current_variant["alts"].append(alt)
                current_variant["rest"].append(rest)
            else:
                # Write out the current_variant
                write_variant(fout, current_variant)

                # Start processing the next variant
                current_variant = {"chrom": chrom, "pos": pos, "ref": ref, "alts": [alt], "rest": [rest]}

    # Handle the last variant
    write_variant(fout, current_variant)
