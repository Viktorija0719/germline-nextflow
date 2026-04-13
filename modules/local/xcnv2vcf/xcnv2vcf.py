#!/usr/bin/env python3
"""Convert a XHMM DATA.xcnv file to one VCF per sample."""
import sys
import csv
from collections import defaultdict


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: xcnv2vcf.py <DATA.xcnv>")

    xcnv_file = sys.argv[1]
    samples = defaultdict(list)

    with open(xcnv_file) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            samples[row['SAMPLE']].append(row)

    vcf_header = [
        "##fileformat=VCFv4.1",
        "##source=XHMM",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of the variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=Q_SOME,Number=1,Type=Float,Description="XHMM Q_SOME quality score">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    ]

    for sample, calls in samples.items():
        calls.sort(key=lambda c: (c['CHR'], int(c['INTERVAL'].split(':')[1].split('-')[0])))
        with open('%s.xcnv.vcf' % sample, 'w') as out:
            out.write('\n'.join(vcf_header) + '\n')
            for c in calls:
                chrom = c['CHR']
                _, span = c['INTERVAL'].split(':', 1)
                start_s, end_s = span.split('-')
                start  = int(start_s)
                end    = int(end_s)
                svtype = 'DUP' if c['CNV'] == 'DUP' else 'DEL'
                svlen  = (end - start) if svtype == 'DUP' else -(end - start)
                info   = 'SVTYPE=%s;END=%d;SVLEN=%d;Q_SOME=%s' % (svtype, end, svlen, c['Q_SOME'])
                out.write('%s\t%d\t.\tN\t<%s>\t.\tPASS\t%s\n' % (chrom, start, svtype, info))


if __name__ == '__main__':
    main()
