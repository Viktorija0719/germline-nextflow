#!/usr/bin/env python3
"""
Convert BED produced by `bedtools merge -c 8 -o collapse` back to a minimal VCF.

Input columns: CHROM, START (0-based), END, INFO-collapse (comma-joined from multiple variants)
Output: standard VCF with SVTYPE, END, SVLEN in INFO.
"""
import sys
import argparse
import re


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('-i', '--input',  required=True, help='Input BED file')
    ap.add_argument('-o', '--output', default='-',   help='Output VCF (default: stdout)')
    ap.add_argument('-v', '--vcf',    action='store_true', help='Write VCF format (always true)')
    return ap.parse_args()


def get_svtype(info):
    m = re.search(r'SVTYPE=([^;,\s]+)', info)
    return m.group(1) if m else 'BND'


def get_svlen(info, start, end):
    m = re.search(r'SVLEN=(-?\d+)', info)
    return m.group(1) if m else str(int(end) - int(start))


def main():
    args = parse_args()
    out  = open(args.output, 'w') if args.output != '-' else sys.stdout

    out.write('##fileformat=VCFv4.2\n')
    out.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    out.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">\n')
    out.write('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">\n')
    out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    with open(args.input) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            first_info = parts[3].split(',')[0] if len(parts) > 3 else ''
            svtype = get_svtype(first_info)
            svlen  = get_svlen(first_info, start, end)
            pos    = str(int(start) + 1)
            info   = f'SVTYPE={svtype};END={end};SVLEN={svlen}'
            out.write(f'{chrom}\t{pos}\t.\tN\t<{svtype}>\t.\tPASS\t{info}\n')

    if args.output != '-':
        out.close()


if __name__ == '__main__':
    main()
