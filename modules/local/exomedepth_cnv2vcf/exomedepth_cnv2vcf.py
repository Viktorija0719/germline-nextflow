#!/usr/bin/env python3
"""Convert a single ExomeDepth per-sample CNV CSV to VCF."""
import sys
import csv


def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: exomedepth_cnv2vcf.py <calls.csv> <sample_id>")

    csv_file  = sys.argv[1]
    sample_id = sys.argv[2]
    out_file  = '%s.exomedepth.cnv.vcf' % sample_id

    vcf_header = [
        "##fileformat=VCFv4.1",
        "##source=ExomeDepth",
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End coordinate of the variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
        '##INFO=<ID=BF,Number=1,Type=Float,Description="ExomeDepth Bayes Factor supporting the CNV call">',
        '##INFO=<ID=NEXONS,Number=1,Type=Integer,Description="Number of target regions spanned by the CNV">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO',
    ]

    with open(csv_file) as fh, open(out_file, 'w') as out:
        out.write('\n'.join(vcf_header) + '\n')
        for row in csv.DictReader(fh):
            chrom  = row['chromosome'].strip('"')
            start  = int(row['start'].strip('"'))
            end    = int(row['end'].strip('"'))
            ctype  = row['type'].strip('"').lower()
            svtype = 'DEL' if 'del' in ctype else 'DUP'
            svlen  = (end - start) if svtype == 'DUP' else -(end - start)
            bf     = row['BF'].strip('"')
            nexons = row['nexons'].strip('"')
            info   = 'SVTYPE=%s;END=%d;SVLEN=%d;BF=%s;NEXONS=%s' % (svtype, end, svlen, bf, nexons)
            out.write('%s\t%d\t.\tN\t<%s>\t.\tPASS\t%s\n' % (chrom, start, svtype, info))


if __name__ == '__main__':
    main()
