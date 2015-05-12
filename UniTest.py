#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

import random
import string

if __name__ == '__main__':

    msg = ""

    ###############################################################################################

    print ("Testing MakeBlastDB")
    try:
        from MakeBlastDB import MakeBlastDB
        db_maker = MakeBlastDB()
        db_path = db_maker (ref_path="./test/Reference.fa" , db_path="./test/Reference")
        msg += "PASS\t MakeBlastDB test\n"
    except Exception as E:
        msg += "FAIL\t MakeBlastDB test\t {}\n".format(E)

    ###############################################################################################

    print ("Testing BlastHit")
    try:
        from BlastHit import BlastHit
        hit_list = []
        for i in range (100):
            hit_list .append(BlastHit (
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
                bscore = random.uniform(0, 1),
                qseq = ''.join(random.choice(["A","T","C","G"]) for _ in range(20))))
        assert len(hit_list) == 100, "Incorect number of hit generated"
        msg += "PASS\t BlastHit test\n"
    except Exception as E:
        msg += "FAIL\t BlastHit test\t {}\n".format(E)

    ###############################################################################################

    print ("Testing MakeBlastn")
    try:
        from MakeBlastn import MakeBlastn
        for _ in range (50):
            blast_maker = MakeBlastn (
                task = random.choice(['blastn', 'blastn-short', 'dc-megablast', 'megablast', 'rmblastn']),
                evalue = random.uniform(0, 1),
                best_per_query_seq = random.choice([True, False]))
            hits = blast_maker(query_path = "./test/query_sample.fa", db_path = db_path)

        msg += "PASS\t MakeBlastn test\n"
    except Exception as E:
        msg += "FAIL\t MakeBlastn test\t {}\n".format(E)

    print (msg)

    ###############################################################################################
