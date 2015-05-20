# -*- coding: utf-8 -*-

"""
@package    pyBlast
@brief      Test function to be used with python package pytest
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library packages import
import sys, os, string, tempfile
from random import randint as ri
from random import uniform as rf
from random import choice as rc

# Third party packages import
import pytest

# Import the current working dir in the python path to allow local package imports
sys.path.append(os.getcwd())

# HELPER FUNCTIONS AND CLASSES ####################################################################

def rs (length):
    """Generate a random string"""
    return ''.join(rc(string.ascii_lowercase + string.digits) for _ in range(length))

def rDNA (length):
    """Generate a random DNA string"""
    return ''.join(rc(["A","T","C","G"]) for _ in range(length))


# TESTS BLAST HIT #################################################################################

# Define parameters for the test_BlastHit function with a pytest decorator
@pytest.mark.parametrize("q_id, s_id, identity, length, mis, gap, q_start, q_end, s_start, s_end, evalue, bscore, q_seq", [
    (rs(10), rs(10), rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,10), rf(0,100), rDNA (30)),
    pytest.mark.xfail((rs(10), rs(10), -1, ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), -1, ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), ri(1,100), -1, ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), ri(1,100), ri(0,100), -1, ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), ri(1,100), ri(0,100), ri(0,100), -1, ri(0,100), ri(0,100), ri(0,100), rf(0,100), rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), -1, ri(0,100), ri(0,100), rf(0,100), rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), -1, ri(0,100), rf(0,100), rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), -1, rf(0,100), rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), -1, rf(0,100), rDNA (30))),
    pytest.mark.xfail((rs(10), rs(10), rf(0,100), ri(1,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), ri(0,100), rf(0,100), -1, rDNA (30)))
    ])

def test_BlastHit(q_id, s_id, identity, length, mis,gap, q_start, q_end, s_start, s_end, evalue, bscore, q_seq):
    """Test BlastHit success and with failure with various parameters"""
    from BlastHit import BlastHit
    BlastHit(q_id, s_id, identity, length, mis,gap, q_start, q_end, s_start, s_end, evalue, bscore, q_seq)

# TESTS BLASTN ####################################################################################

class rand_subject_query_files (object):
    """Small helper class to generates fasta file on the fly"""

    DNA_COMPLEMENT = {"A":"T","T":"A","C":"G","G":"C"}

    def __init__(self, len_subject, len_query, random_query=False):
        """"""
        print ("Create reference and query fasta files")
        self.s_seq = rDNA(len_subject)

        # Generate either a random query sequence in forward or reverse orientation from the subject sequence
        if random_query:
            self.s_orient, self.s_start, self.s_end = True, 0, 0
            self.q_seq = rDNA(len_query)
        elif ri(0,1):
            self.s_orient = True
            self.s_start = ri(0, len_subject-len_query)
            self.s_end = self.s_start + len_query
            self.q_seq = self.s_seq [self.s_start:self.s_end]
        else:
            self.s_orient = False
            self.s_end = ri(0, len_subject-len_query)
            self.s_start = self.s_end + len_query
            self.q_seq = "".join ([self.DNA_COMPLEMENT[base] for base in self.s_seq[self.s_start-1:self.s_end-1:-1]])

        # Write files in temporary directories
        self.s_path = tempfile.mkstemp()[1]
        with open (self.s_path, "w") as fp:
            fp.write (">ref\n{}\n".format(self.s_seq))

        self.q_path = tempfile.mkstemp()[1]
        with open (self.q_path, "w") as fp:
            fp.write (">query\n{}\n".format(self.q_seq))

    # Enter and exit are defined to use the with statement
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        """Destructor to remove the database and unziped fasta files if needed"""
        print ("Destroy fasta files")
        os.remove(self.s_path)
        os.remove(self.q_path)

def test_Blastn_hit():
    """Test Blastn class with simulated datasets when query are generated from the subject sequence"""

    # Import and create temp fasta
    from Blastn import Blastn

    for _ in range (100):
        with rand_subject_query_files(200, 20, random_query=False) as r:
            with Blastn(r.s_path) as blastn:
                hit_list = blastn(r.q_path)
                assert len(hit_list) == 1
                assert hit.s_orient == r.s_orient
                assert hit.s_start == r.s_start
                assert hit.s_end == r.s_end
                assert hit.q_seq == r.q_seq

def test_Blastn_nohit():
    """Test Blastn class with simulated datasets when query are randomly generated"""

    # Import and create temp fasta
    from Blastn import Blastn

    for _ in range (100):
        with rand_subject_query_files(200, 20, random_query=True) as r:
            with Blastn(r.s_path) as blastn:
                hit_list = blastn(r.q_path)
                assert hit_list == None
