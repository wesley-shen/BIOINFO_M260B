import sys
import os
sys.path.insert(0, os.path.abspath(".."))
sys.path.insert(0, os.path.abspath("../.."))
import numpy as np
from os.path import join
import time
from BIOINFO_M260B.helpers import read_reads, read_reference, pretty_print_aligned_reads_with_ref
from collections import defaultdict

def trivial_algorithm(paired_end_reads, ref):
    """
    This is a functional aligner, but it's a huge simplification that
    generate a LOT of potential bugs.  It's also very slow.

    Read the spec carefully; consider how the paired-end reads are
    generated, and ideally, write your own algorithm
    instead of trying to tweak this one (which isn't very good).

    :param paired_end_reads: Paired-end reads generated from read_reads
    :param ref: A reference genome generated from read_reference
    :return: 2 lists:
                1) a list of alignment locations for each read (all_alignment_locations).
                    The list gives the starting position of the minimum-mismatch alignment of both reads.
                2) a list of the paired-end reads set so that both reads are in their optimal orientation
                   with respect to the reference genome.
    """
    all_read_alignment_locations = []
    output_read_pairs = []
    count = 0
    start = time.clock()

    # build up the hash table for the reference genome
    oligomer_len = 17
    ref_dict = defaultdict(set)
    for i in range(len(ref)-oligomer_len):
        ref_dict[ref[i:i+oligomer_len]].add(i)

    for read_pair in paired_end_reads:
        count += 1
        read_alignment_locations = []
        output_read_pair = []
        if count % 10 == 0:
            time_passed = (time.clock()-start)/60
            print ('{} reads aligned'.format(count), 'in {:.3} minutes'.format(time_passed))
            remaining_time = time_passed/count*(len(paired_end_reads)-count)
            print ('Approximately {:.3} minutes remaining'.format(remaining_time))
        for read in read_pair:
            # break the read into three segments
            read_segs = [read[0:oligomer_len], read[oligomer_len:oligomer_len*2], read[oligomer_len*2:len(read)]]
            min_mismatches = oligomer_len
            min_mismatch_location = -1
            for i in range(len(read_segs)):
                targeted_locations = ref_dict[read_segs[i]]
                min_mismatches = oligomer_len*2
                min_mismatches_location = -1
                if targeted_locations != []:                    
                    if (i == 0):
                        for j in targeted_locations:
                            mismatches = [1 if read[oligomer_len+k] != ref[j+oligomer_len+k] for k in range(0, len(read)-oligomer_len)]
		    elif (i == 1):
		        for j in targeted_location:
			    mismatch_1 = [1 if read[k] != ref[j-oligomer+k] for k in range(0, oligomer_len)]
			    mismatch_2 = [1 if read[oligomer_len*2-1+k] != ref[j+oligomer_len+k] for k in range(0, oligomer_len)]
			    mismatches = mismatch_1+mismatch
		    else:
		        for j in targeted_location:
			    mismatches = [1 if read[k] != ref[j-oligomer*2+1] for k in range(0, len(read)-oligomer_len)]
		else:
                    continue
            # for i in range(len(ref) - len(read)):
            #     mismatches = [1 if read[j] != ref[i + j] else 0 for j in range(len(read))]
            #     n_mismatches = sum(mismatches)
            #     # The above line should be familiar to Python users, but bears  some explanation for
            #     # people who are getting started with it. The "mismatches = ..." line
            #     # is called a "list comprehension. Basically, this is a short way of writing the loop:
            #     #
            #     # n_mismatches = 0
            #     # for j in range(len(read)):
            #     # if read[j] != ref[i+j]:
            #     #         n_mismatches += 1
            #     #
            #     # The first line creates a list which has a 1 for every mismatch and a 0 for every match.
            #     # The second line sums the list created by the first line, which counts the number of mismatches.
            #     if n_mismatches < min_mismatches:
            #         min_mismatches = n_mismatches
            #         min_mismatch_location = i

            # reversed_read = read[::-1]
            # for i in range(len(ref) - 50):
            #     mismatches = [1 if reversed_read[j] != ref[i + j] else 0 for j in range(len(read))]
            #     n_mismatches = sum(mismatches)
            #     if n_mismatches < min_mismatches:
            #         min_mismatches = n_mismatches
            #         min_mismatch_location = i
            #         read = reversed_read
            # if (mismatches < 3):
            #     read_alignment_locations.append(min_mismatch_location)
            #     output_read_pair.append(read)
            # # Note that there are some huge potential problems here.

        all_read_alignment_locations.append(read_alignment_locations)
        output_read_pairs.append(output_read_pair)
    return all_read_alignment_locations, output_read_pairs


if __name__ == "__main__":
    data_folder = 'practice_W_1'
    input_folder = join('../data/', data_folder)
    f_base = '{}_chr_1'.format(data_folder)
    reads_fn = join(input_folder, 'reads_{}.txt'.format(f_base))
    start = time.clock()
    input_reads = read_reads(reads_fn)
    # This will take a while; you can use an array slice for example:
    #
    #   input_reads = reads[:300]
    #
    # to generate some data quickly.

    reference_fn = join(input_folder, 'ref_{}.txt'.format(f_base))
    reference = read_reference(reference_fn)
    alignments, reads = trivial_algorithm(input_reads, reference)
    print (alignments)
    print (reads)
    output_str = pretty_print_aligned_reads_with_ref(reads, alignments, reference)
    output_fn = join(input_folder, 'aligned___{}.txt'.format(f_base))
    with(open(output_fn, 'w')) as output_file:
        output_file.write(output_str)
