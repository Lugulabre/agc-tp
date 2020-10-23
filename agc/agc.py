#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Maël Pretet"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Pretet"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Pretet"
__email__ = "pretetmael@gmail.com"
__status__ = "Developpement"


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


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    if amplicon_file.endswith("gz"):
        with gzip.open(amplicon_file, "rb") as filin:
            seq = b""
            for line in filin:
                if line.startswith(b">"):
                    if len(seq) >= minseqlen:
                        yield seq.decode('ascii')
                    seq = b""
                else:
                    seq += line.strip()

            yield seq.decode('ascii')
    else:
        with open(amplicon_file, "r") as filin:
            seq = ""
            for line in filin:
                if line.startswith(">"):
                    if len(seq) >= minseqlen:
                        yield seq
                    seq = ""
                else:
                    seq += line.strip()

            yield seq


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    dict_amplicon = {}
    for amplicon in read_fasta(amplicon_file, minseqlen):
        if not amplicon in dict_amplicon.keys():
            dict_amplicon[amplicon] = 1
        else:
            dict_amplicon[amplicon] += 1

    for key, value in sorted(dict_amplicon.items(), key=lambda x: x[1], reverse = True):
        if value >= mincount:
            yield [key, value]


def get_chunks(sequence, chunk_size):
    list_sub_seq = []
    for i in range(0, len(sequence), chunk_size):
        if i+chunk_size < len(sequence):
            list_sub_seq.append(sequence[i:i+chunk_size])

    if len(list_sub_seq) >= 4:
        return list_sub_seq
    else:
        raise ValueError


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def cut_kmer(sequence, kmer_size):
    for i in range(0, len(sequence)-kmer_size+1):
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    for kmer in cut_kmer(sequence, kmer_size):
        if not kmer in kmer_dict.keys():
            kmer_dict[kmer] = [id_seq]
        else:
            kmer_dict[kmer].append(id_seq)

    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    return [i[0] for i in Counter([ids 
        for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict 
        for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list):
    count_base = 0
    for i in range(0, len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            count_base += 1
    return count_base / len(alignment_list[0]) * 100


def detect_chimera(perc_identity_matrix):
    sum_stdev = 0
    bool1_seq = False
    bool2_seq = False
    for i in range(0, len(perc_identity_matrix)):
        sum_stdev += statistics.stdev(perc_identity_matrix[i])
        if perc_identity_matrix[i][0] > perc_identity_matrix[i][1]:
            bool1_seq = True
        if perc_identity_matrix[i][1] > perc_identity_matrix[i][0]:
            bool2_seq = True

    if sum_stdev / len(perc_identity_matrix) > 5 and bool1_seq and bool2_seq:
        return True
    
    return False


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    list_non_chimere = []
    kmer_dict = {}
    id_ref = 0

    for i, amplicon in enumerate(dereplication_fulllength(amplicon_file, minseqlen, mincount)):
        is_chim = True
        chunks = get_chunks(amplicon[0], chunk_size)
        mates = [search_mates(kmer_dict, amplicon, kmer_size) for chunk in chunks]
        common_id = common(mates[0], mates[1])
        for same in range(2, len(mates)):
            common_id = common(common_id, mates[same])

        if len(common_id) > 1:
            perc_identity_matrix = [[] for g in range(len(chunks))]
            for m in common_id[0:2]:
                db_seq = get_chunks(list_non_chimere[m], chunk_size)
                for j, chunk in enumerate(chunk_list):
                    #print(chunk, db_seq[j])
                    perc_identity_matrix[j].append(get_identity(
                        nw.global_align(chunk, db_seq[j], 
                            gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '../agc')) + "/MATCH")))
            is_chim = detect_chimera(perc_identity_matrix)
        else:
            is_chim = False

        if not is_chim:
            kmer_dict = get_unique_kmer(kmer_dict, amplicon[0], id_ref, kmer_size)
            list_non_chimere.append(amplicon[0])
            id_ref += 1
            yield(amplicon)


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    amplicon_list = chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    #amplicon_list = amplicon_file
    greedy_cluster = []
    for i, amplicon in enumerate(amplicon_list):
        print("Tour n°",i, "\r", end = "")
        if not greedy_cluster:
            greedy_cluster.append((amplicon[0], amplicon[1]))
        else:
            flag_dico = True
            for tupple in greedy_cluster:
                align_tmp = nw.global_align(amplicon[0], tupple[0], 
                            gap_open=-1, gap_extend=-1,
                            matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),
                            '../agc')) + "/MATCH")
                if get_identity(align_tmp) > 97:
                    #tupple[1] += amplicon[1]
                    flag_dico = False

            if flag_dico:
                greedy_cluster.append((amplicon[0], amplicon[1]))

    return greedy_cluster


def fill(text, width=80):
    '''Diviser le texte avec un retour à la ligne pour respecter le format fasta
    '''
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def write_OTU(OTU_list, output_file):
    '''Ecriture des OTU dans un fichier de sortie
    '''
    with open(output_file, "w") as filout:
        for i in range(0,len(OTU_list)):
            filout.write(">OTU_{} occurrence:{}\n{}\n".format(i+1, OTU_list[i][1], fill(OTU_list[i][0])))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    OTU_cluster = abundance_greedy_clustering(args.amplicon_file, args.minseqlen,
                  args.mincount, args.chunk_size, args.kmer_size)
    write_OTU(OTU_cluster, args.output_file)

if __name__ == '__main__':
    main()