#!/usr/bin/env python3
"""Convert a single ExomeDepth per-sample CNV CSV to VCF.

Variants with Bayes Factor below --lowqc_bf (default 9.0) are tagged
FILTER=lowQC; all others receive FILTER=PASS.  The BF value is also
written as the VCF QUAL field so downstream tools can use it directly.
"""
import sys
import csv
import argparse


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('csv_file',   help='ExomeDepth CNV calls CSV')
    ap.add_argument('sample_id',  help='Sample identifier (used for output filename)')
    ap.add_argument(
        '--lowqc_bf', type=float, default=9.0, metavar='BF',
        help='Bayes Factor threshold: calls below this value receive FILTER=lowQC (default: %(default)s)',
    )
    return ap.parse_args()


def main():
    args = parse_args()

    out_file = f'{args.sample_id}.exomedepth.cnv.vcf'

    vcf_header = [
        '##fileformat=VCFv4.2',
        '##source=ExomeDepth',
        '##FILTER=<ID=PASS,Description="All filters passed">',
        (f'##FILTER=<ID=lowQC,Description="Low-confidence CNV call '
         f'(ExomeDepth Bayes Factor < {args.lowqc_bf})">'),
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of the variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=BF,Number=1,Type=Float,Description="ExomeDepth Bayes Factor supporting the CNV call">',
        '##INFO=<ID=NEXONS,Number=1,Type=Integer,Description="Number of target regions spanned by the CNV">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    ]

    with open(args.csv_file) as fh, open(out_file, 'w') as out:
        out.write('\n'.join(vcf_header) + '\n')
        for row in csv.DictReader(fh):
            chrom  = row['chromosome']
            start  = int(row['start'])
            end    = int(row['end'])
            ctype  = row['type'].lower()
            svtype = 'DEL' if 'del' in ctype else 'DUP'
            svlen  = (end - start) if svtype == 'DUP' else -(end - start)
            bf_raw = row['BF']
            nexons = row['nexons']
            info   = f'SVTYPE={svtype};END={end};SVLEN={svlen};BF={bf_raw};NEXONS={nexons}'

            try:
                bf_val = float(bf_raw)
                qual   = f'{bf_val:.2f}'
                filt   = 'lowQC' if bf_val < args.lowqc_bf else 'PASS'
            except (ValueError, TypeError):
                qual = '.'
                filt = 'PASS'

            out.write(f'{chrom}\t{start}\t.\tN\t<{svtype}>\t{qual}\t{filt}\t{info}\n')


if __name__ == '__main__':
    main()
