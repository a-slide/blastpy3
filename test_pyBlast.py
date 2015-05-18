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
import sys, os, random, string

# Third party packages import
import pytest

# Import the current working dir in the python path to allow local package imports
sys.path.append(os.getcwd())

# Start tests
def test_MakeBlastDB():
    from MakeBlastDB import MakeBlastDB
    db_maker = MakeBlastDB()
    assert(db_maker (
            ref_path= "./test/Reference.fa" ,
            db_path= "../test/Reference"))

def test_BlastHit():
    from BlastHit import BlastHit
    for _ in range (100):
        assert(BlastHit(
                q_id = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(10)),
                s_id = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(10)),
                identity = random.uniform(0, 1),
                length = random.randint(1, 100),
                mis = random.randint(0, 100),
                gap = random.randint(0, 100),
                q_start = random.randint(0, 100),
                q_end = random.randint(0, 100),
                s_start = random.randint(0, 100),
                s_end = random.randint(0, 100),
                evalue = random.uniform(0, 1),
                bscore = random.uniform(0, 100),
                qseq = ''.join(random.choice(["A","T","C","G"]) for _ in range(20))))

def test_MakeBlastn():
    from MakeBlastn import MakeBlastn
    for _ in range (10):
        blast_maker = MakeBlastn(
            task = random.choice(['blastn', 'blastn-short', 'dc-megablast', 'megablast', 'rmblastn']),
            evalue = random.uniform(0, 1),
            best_per_query_seq = random.choice([True, False]))

        assert(blast_maker(
                query_path="./test/query_sample.fa",
                db_path="./test/Reference" ))
