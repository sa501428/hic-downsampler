import sys

import strawC
import struct

def get_chromosomes(filepath: str):
    chrom_dot_sizes = strawC.getChromosomes(filepath)
    chromosomes = []
    for chromosome in chrom_dot_sizes:
        chrom = chromosome.name
        if chrom.lower() == 'all' or chromosome.index == 0:
            print("Skipping chromosome", chromosome.name, chromosome.index)
            continue
        chromosomes.append(chromosome)
    chromosomes.sort(key=lambda chrom: chrom.index)
    names = []
    for chrom in chromosomes:
        names.append(chrom.name)
    return names


def write_short_format_line(outfile, chr1, x1, chr2, y1, score):
    outfile.write("{0} {1} {2} {3} {4}\n".format(chr1, x1, chr2, y1, int(score)))


def write_short_binary_format_line(outfile, chr1, binX, chr2, binY, counts):
    for x in [int(chr1), int(binX), int(chr2), int(binY), float(counts)]:
        outfile.write(struct.pack('b', x))


def write_all_contacts_to_binary_format(outfile, chr1, chr2, result):
    for contact in result:
        write_short_binary_format_line(outfile, chr1, contact.binX, chr2, contact.binY, contact.counts)


def write_all_contacts_to_short_format(outfile, chr1, chr2, result):
    for contact in result:
        write_short_format_line(outfile, chr1, contact.binX, chr2, contact.binY, contact.counts)


def extract_all_raw_contacts(hicfile, resolution, outfile_name):
    chroms = get_chromosomes(hicfile)
    print(chroms)
    use_short_binary_output = outfile_name.endswith(".bn")

    outfile = open(outfile_name, 'w')
    for i in range(len(chroms)):
        chr1 = chroms[i]
        for j in range(i, len(chroms)):
            chr2 = chroms[j]
            result = strawC.strawC("observed", "NONE", hicfile, str(chr1), str(chr2), 'BP', resolution)
            if use_short_binary_output:
                write_all_contacts_to_binary_format(outfile, i + 1, j + 1, result)
            else:
                write_all_contacts_to_short_format(outfile, chr1, chr2, result)
    outfile.close()


# just_extract.py input.hic resolution outfile
if __name__ == "__main__":
    extract_all_raw_contacts(sys.argv[1], int(sys.argv[2]), sys.argv[3])
