# -*- coding: utf-8 -*-

"""
@package    pyBlast
@brief      Test functions to be used with python package pytest
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# IMPORTS #########################################################################################

# Standard library packages import
import sys, os, string, tempfile
from random import randint as ri
from random import uniform as rf
from random import choice as rc

# Third party packages import
import pytest

# Import the current working dir in the python path to allow local package imports
sys.path.append(os.getcwd())

# local package imports
from BlastHit import BlastHit
from Blastn import Blastn

# HELPER FUNCTIONS AND CLASSES ####################################################################

"""Generate a random string"""
def rs (length):
    return ''.join(rc(string.ascii_lowercase + string.digits) for _ in range(length))

"""Generate a random DNA string"""
def rDNA (length):
    return ''.join(rc(["A","T","C","G"]) for _ in range(length))

"""Helper class to generates fasta file on the fly"""
class rand_subject_query_files (object):
    DNA_COMPLEMENT = {"A":"T","T":"A","C":"G","G":"C"}

    def __init__(self, len_subject, len_query, random_query=False):
        # Generate a random subject fasta file
        self.s_seq = rDNA(len_subject)
        self.s_path = tempfile.mkstemp()[1]
        with open (self.s_path, "w") as fp:
            fp.write (">ref\n{}\n".format(self.s_seq))

        # Generate either a random query sequence in forward or reverse orientation from the subject sequence
        if random_query:
            self.q_seq = rDNA(len_query)
        elif ri(0,1):
            self.s_orient = "+"
            self.s_start = ri(0, len_subject-len_query)
            self.s_end = self.s_start + len_query
            self.q_seq = self.s_seq [self.s_start:self.s_end]
        else:
            self.s_orient = "-"
            self.s_end = ri(0, len_subject-len_query)
            self.s_start = self.s_end + len_query
            self.q_seq = "".join ([self.DNA_COMPLEMENT[base] for base in self.s_seq[self.s_start-1:self.s_end-1:-1]])

        self.q_path = tempfile.mkstemp()[1]
        with open (self.q_path, "w") as fp:
            fp.write (">query\n{}\n".format(self.q_seq))

    def __enter__(self):
        print ("Create reference and query fasta files")
        return self

    def __exit__(self, type, value, traceback):
        print ("Destroy fasta files")
        os.remove(self.s_path)
        os.remove(self.q_path)

# TESTS BLAST HIT #################################################################################

# Define parameters for the test_BlastHit function with a pytest decorator
@pytest.mark.parametrize("identity, length, mis, gap, q_start, q_end, s_start, s_end, evalue, bscore", [
    (rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,10), rf(0,100)),
    pytest.mark.xfail((-1, ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100))),
    pytest.mark.xfail((rf(0,100), -1, ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100))),
    pytest.mark.xfail((rf(0,100), ri(1,100), -1, ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100))),
    pytest.mark.xfail((rf(0,100), ri(1,100), ri(0,100), -1, ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100))),
    pytest.mark.xfail((rf(0,100), ri(1,100), ri(0,100), ri(0,100), -1, ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100))),
    pytest.mark.xfail((rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), -1, ri(0,100), ri(0,100), rf(0,100), rf(0,100))),
    pytest.mark.xfail((rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), -1, ri(0,100), rf(0,100), rf(0,100))),
    pytest.mark.xfail((rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), -1, rf(0,100), rf(0,100))),
    pytest.mark.xfail((rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), -1, rf(0,100))),
    pytest.mark.xfail((rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), -1))
    ])

# Test BlastHit success and with failure with various parameters
def test_BlastHit(identity, length, mis, gap, q_start, q_end, s_start, s_end, evalue, bscore):
    BlastHit(identity=identity, length=length, mis=mis, gap=gap, q_start=q_start, q_end=q_end,
        s_start=s_start, s_end=s_end, evalue=evalue, bscore=bscore)

# TESTS BLASTN ####################################################################################

@pytest.fixture (params=['blastn', 'blastn-short', 'dc-megablast', 'megablast', 'rmblastn'])
def task (request):
    return request.param

@pytest.fixture (params=[False, pytest.mark.xfail(True)], ids=["Queries from Subject", "Random queries"])
def random_query (request):
    return request.param

"""Test Blastn class with simulated datasets when query are generated from the subject sequence"""
def test_Blastn(task, random_query):
    # Loop to try different random combinations
    for _ in range (5):
        with rand_subject_query_files(500, 50, random_query=random_query) as r:

            # Instantiate the database and perform a blastn
            with Blastn(r.s_path) as blastn:
                hit_list = blastn(query_path=r.q_path, task=task, best_query_hit=True, evalue=1)

                # Test values in hit_list
                assert len(hit_list) >= 1
                assert hit_list[0].s_orient == r.s_orient
                assert hit_list[0].s_start == r.s_start
                assert hit_list[0].s_end == r.s_end
                assert hit_list[0].q_seq == r.q_seq
