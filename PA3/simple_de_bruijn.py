from os.path import join
import sys
import time
from collections import defaultdict, Counter
import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
# from BIOINFO_M260B.helpers import read_reads
from helpers import read_reads

"""
def read_assembly_reads(read_fn):
    reads = read_reads(read_fn)
    # POTENTIAL
    output_reads = [_[0] for _ in reads]
    # Only taking one end of the read works okay, but
    # this is an obvious area for improvement.
    return output_reads

def simple_de_bruijn(sequence_reads, k):

    # Creates A simple DeBruijn Graph with nodes that correspond to k-mers of size k.
    # :param sequence_reads: A list of reads from the genome
    # :param k: The length of the k-mers that are used as nodes of the DeBruijn graph
    # :return: A DeBruijn graph where the keys are k-mers and the values are the set
    #             of k-mers that

    de_bruijn_counter = defaultdict(Counter)
    # You may also want to check the in-degree and out-degree of each node
    # to help you find the beginnning and end of the sequence.
    for read in sequence_reads:
        # Cut the read into k-mers
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])

    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    # THE COVERAGE is 30x, so maybe at least 20 times not 2 times
    de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 1}
                       for key in de_bruijn_counter}

    # This line removes the empty nodes from the DeBruijn graph
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}
    # print(de_bruijn_graph.items())
    return de_bruijn_graph

def de_bruijn_reassemble(de_bruijn_graph):

    # Traverses the DeBruijn Graph created by simple_de_bruijn and
    # returns contigs that come from it.
    # :param de_bruijn_graph: A De Bruijn Graph
    # :return: a list of the


    assembled_strings = []
    while True:
        n_values = sum([len(de_bruijn_graph[k]) for k in de_bruijn_graph])
        print(n_values)
        if n_values == 0:
            break
        # outgoing_degs = list(de_bruijn_graph.keys())
        # incoming_degs = [_[1] for _ in de_bruijn_graph.items()]
        good_starts = [k for k in de_bruijn_graph if de_bruijn_graph[k]]
        # good_starts = [k for k in de_bruijn_graph if de_bruijn_graph[k] and outgoing_degs.count(k)>incoming_degs.count(k)]
        # You may want to find a better start
        # position by looking at in and out-degrees
        # but this will work okay.
        # LOOK FOR nodes with more outgoing edges than incoming edges
        current_point = good_starts[0]
        assembled_string = current_point
        while True:
            try:
                next_values = de_bruijn_graph[current_point]
                # pop() might be too random. Choose the most possible one.
                next_edge = next_values.pop()
                assembled_string += next_edge[-1]
                de_bruijn_graph[current_point] = next_values
                current_point = next_edge
            except KeyError:
                assembled_strings.append(assembled_string)
                break
    return assembled_strings

if __name__ == "__main__":
    # chr_name = 'practice_A_2'
    chr_name = 'hw3all_A_3'
    input_folder = './{}/{}'.format("data",chr_name)
    reads_fn = join(input_folder, 'reads_{}_chr_1.txt'.format(chr_name))
    # reads_fn = join(input_folder, 'reads_{}.txt'.format(chr_name))
    reads = read_assembly_reads(reads_fn)
    db_graph = simple_de_bruijn(reads, 25)
    for k in list(db_graph.keys()):
        print (k, db_graph[k])

    output = de_bruijn_reassemble(db_graph)
    output_fn_end = 'assembled_{}.txt'.format(chr_name)
    output_fn = join(input_folder, output_fn_end)
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + chr_name + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(output))
"""

def read_assembly_reads(read_fn):
    reads = read_reads(read_fn)
    output_reads = [_[0] for _ in reads] 
    output_reads += [_[1] for _ in reads]
    # Only taking one end of the read works okay, but
    # this is an obvious area for improvement.
    return output_reads


def simple_de_bruijn(sequence_reads, k):
    """
    Creates A simple DeBruijn Graph with nodes that correspond to k-mers of size k.
    :param sequence_reads: A list of reads from the genome
    :param k: The length of the k-mers that are used as nodes of the DeBruijn graph
    :return: A DeBruijn graph where the keys are k-mers and the values are the set
                of k-mers that
    """
    de_bruijn_counter = defaultdict(Counter)
    # You may also want to check the in-degree and out-degree of each node
    # to help you find the beginnning and end of the sequence.
    for read in sequence_reads:
        # Cut the read into k-mers
        kmers = [read[i: i + k] for i in range(len(read) - k)]
        for i in range(len(kmers) - 1):
            pvs_kmer = kmers[i]
            next_kmer = kmers[i + 1]
            de_bruijn_counter[pvs_kmer].update([next_kmer])

    # This line removes the nodes from the DeBruijn Graph that we have not seen enough.
    de_bruijn_graph = {key: {val for val in de_bruijn_counter[key] if de_bruijn_counter[key][val] > 2}
                       for key in de_bruijn_counter}

    # This line removes the empty nodes from the DeBruijn graph
    de_bruijn_graph = {key: de_bruijn_graph[key] for key in de_bruijn_graph if de_bruijn_graph[key]}
    return de_bruijn_graph


def de_bruijn_reassemble(de_bruijn_graph):
    """
    Traverses the DeBruijn Graph created by simple_de_bruijn and
    returns contigs that come from it.
    :param de_bruijn_graph: A De Bruijn Graph
    :return: a list of the
    """
    assembled_strings = []
    while True:
        n_values = sum([len(de_bruijn_graph[k]) for k in de_bruijn_graph])
        print(n_values)
        if n_values == 0:
            break
        outgoing_degs = list(de_bruijn_graph.keys())
        incoming_degs = [_[1] for _ in de_bruijn_graph.items()]
        good_start = ''
        for k in de_bruijn_graph:
            if de_bruijn_graph[k] and (outgoing_degs.count(k) > incoming_degs.count(k)) :
                good_start = k
                break
        # You may want to find a better start
        # position by looking at in and out-degrees,
        # but this will work okay.
#        current_point = good_starts[0]
        current_point = good_start
        assembled_string = current_point
        while True:
            try:
                next_values = de_bruijn_graph[current_point]
#                print(next_values)
#                next_edge = next_values.pop()
                next_edge = list(next_values)[0] 
                for k in next_values:
                    if k[:len(k)-1] == current_point[1:]:
#                        print(current_point, k)
                        next_edge = k
                next_values.remove(next_edge)
                assembled_string += next_edge[-1]
                de_bruijn_graph[current_point] = next_values
                current_point = next_edge
#            except KeyError:
            except (IndexError, KeyError) as e:
                assembled_strings.append(assembled_string)
                break
    return assembled_strings

if __name__ == "__main__":
    chr_name = 'hw3all_A_3_chr_1'
#    chr_name = 'practice_A_2_chr_1'
    input_folder = '../{}/{}'.format("data", chr_name)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(chr_name))
    reads = read_assembly_reads(reads_fn)
    db_graph = simple_de_bruijn(reads, 30)
    for k in list(db_graph.keys())[:40]:
        print (k, db_graph[k])

    output = de_bruijn_reassemble(db_graph)
    output_fn_end = 'assembled_{}.txt'.format(chr_name)
    output_fn = join(input_folder, output_fn_end)
    with open(output_fn, 'w') as output_file:
        output_file.write('>' + chr_name + '\n')
        output_file.write('>ASSEMBLY\n')
        output_file.write('\n'.join(output))