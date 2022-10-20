import argparse
import sys
import os
import csv
import re


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50, help="Minimum gene length to consider")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes (shine box not included).")
    parser.add_argument('-p', dest='predicted_genes_file', type=str,
                        default=os.curdir + os.sep +"predict_genes.csv",
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=str,
                        default=os.curdir + os.sep + "genes.fna",
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file):
    """Extract the complete genome sequence as a single string
    """
    with open(fasta_file, "rt") as  monfich:
        ephem = ""
        for line in monfich:
            if line.startswith(">"):
                continue
            ephem += line
            ephem = ephem.strip()
        return ephem.upper()

def find_start(start_regex, sequence, start, stop):
    """Find the start codon
    """
    soluce = start_regex.search(sequence, start, stop)
    if  soluce != None:
        return soluce.start(0)
    return None


def find_stop(stop_regex, sequence, start):
    """Find the stop codon
    """
    soluce = stop_regex.finditer(sequence, start)
    for i in soluce:
        if (i.start(0) - start)%3 == 0:
            return i.start(0)
    return None

def has_shine_dalgarno(shine_regex, sequence, start, max_shine_dalgarno_distance):
    """Find a shine dalgarno motif before the start codon
    """
    fin = start - max_shine_dalgarno_distance
    debut = start - 6
    soluce = shine_regex.search(sequence)

    if soluce is not None:
        if debut >= soluce.start(0)  + 1 >= fin and \
            debut >= soluce.end(0) + 1 >= fin:
            return True

    return False

def predict_genes(sequence, start_regex, stop_regex, shine_regex,
                  min_gene_len, max_shine_dalgarno_distance, min_gap):
    """Predict most probable genes
    """
    pos = 0
    mem = []
    while len(sequence) - pos >= min_gap:
        pos = find_start(start_regex, sequence, pos, len(sequence))
        if pos is not None:
            stop = find_stop(stop_regex, sequence, pos)
            if stop is not None:
                taille = stop - pos
                if taille >= min_gene_len:
                    if has_shine_dalgarno(shine_regex, sequence, pos, max_shine_dalgarno_distance):
                        mem.append([pos + 1, stop + 3])
                        pos = pos + stop + 2 + min_gap
                    else:
                        pos += 1
                else:
                    pos += 1
            else:
                pos += 1
        else:
            pos += 1
    return mem


def write_genes_pos(predicted_genes_file, probable_genes):
    """Write list of gene positions
    """
    try:
        with open(predicted_genes_file, "wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_genes(fasta_file, sequence, probable_genes, sequence_rc, probable_genes_comp):
    """Write gene sequence in fasta format
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for index, gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    index + 1, os.linesep,
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
                index = index+1
            for j, gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            index + 1 + j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(kmer):
    """Get the reverse complement"""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in kmer[::-1]])


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """

    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')

    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    args = get_arguments()
    seq_a = read_fasta(args.genome_file)
    result_a = predict_genes(seq_a, start_regex, stop_regex, shine_regex,
                args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)
    seq_b = reverse_complement(seq_a)
    result_b = predict_genes(seq_b, start_regex, stop_regex, shine_regex,
                args.min_gene_len, args.max_shine_dalgarno_distance, args.min_gap)
    write_genes_pos("predict.csv", result_a)
    write_genes("predica.csv", seq_a, result_a, seq_b, result_b)

if __name__ == '__main__':
    main()