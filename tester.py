# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 23:28:12 2017

@author: Wesley
"""

from complex_pileup import *

if __name__ == "__main__":
    donor = "GGTTGACTA"
    ref = "TGTTACGG"
    edit_distance_matrix(ref, donor)
    identify_changes(ref, donor, offset=0)