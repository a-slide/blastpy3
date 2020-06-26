# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
from collections import OrderedDict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BlastHit(object):
    """
    Object oriented class containing informations of one blast hit
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#
    ID_COUNT = 0

    #~~~~~~~CLASS METHODS~~~~~~~#

    @ classmethod
    def next_id (self):
        cur_id = self.ID_COUNT
        self.ID_COUNT +=1
        return cur_id

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self,
        q_id = "query",
        s_id = "subject",
        identity = 100,
        length = 10,
        mis = 0,
        gap = 0,
        q_start = 1,
        q_end = 10,
        s_start = 1,
        s_end = 10,
        evalue = 0,
        bscore = 0,
        q_seq = "NNNNNNNNNN"):
        """
        Create a BlastHit object which is automatically added to the class tracking instance list
        The object with the following parameters are required for object initialisation
        * q_id
            Query sequence name [STR]
        * s_id
            Subject sequence name [STR]
        * identity
            % of identity in the hit [FLOAT 0:100]
        * length
            length of the hit [INT >=1]
        * mis
            Number of mismatch in the hit [INT >=0]
        * gap
            Number of gap in the hit [INT >=0]
        * q_start
            Hit start position of the query (1 based) [INT >=1]
        * q_end
            Hit end position of the query (1 based) [INT >=1]
        * s_start
            Hit start position of the subject (1 based) [INT >=1]
        * s_end
            Hit end position of the subject (1 based) [INT >=1]
        * evalue
            E value of the alignement [FLOAT >=0]
        * bscore
            Bit score of the alignement [FLOAT >=0]
        * q_seq
            Sequence of the query aligned on the subject sequence [STR]
        """

        # Store parameters in self variables
        self.id = self.next_id()
        self.q_id = q_id.decode() if type(q_id) == bytes else q_id
        self.s_id = s_id.decode() if type(s_id) == bytes else s_id
        self.identity = float(identity)
        self.length = int(length)
        self.mis = int(mis)
        self.gap = int(gap)
        self.evalue = float(evalue)
        self.bscore = float(bscore)
        self.q_seq = q_seq.decode() if type(q_seq) == bytes else q_seq

        # Correct coordinates of hit for python 0 based coordinates depending of the orientation
        if int(q_start) < int(q_end):
            self.q_orient = "+"
            self.q_start = int(q_start)-1
            self.q_end = int(q_end)
        else:
            self.q_orient = "-"
            self.q_start = int(q_start)
            self.q_end = int(q_end)-1

        if int(s_start) < int(s_end):
            self.s_orient = "+"
            self.s_start = int(s_start)-1
            self.s_end = int(s_end)
        else:
            self.s_orient = "-"
            self.s_start = int(s_start)
            self.s_end = int(s_end)-1

        # Verify the hit validity
        self._test_arg()

    def __str__(self):
        msg = "HIT {}".format(self.id)
        msg += "\tQuery\t{}:{}-{}({})\n".format(self.q_id, self.q_start, self.q_end, self.q_orient)
        msg += "\tSubject\t{}:{}-{}({})\n".format(self.s_id, self.s_start, self.s_end, self.s_orient)
        msg += "\tLenght : {}\tIdentity : {}%\tEvalue : {}\tBit score : {}\n".format(self.length, self.identity, self.evalue, self.bscore)
        msg += "\tAligned query seq : {}\n".format(self.q_seq)
        return (msg)

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def get_report (self, full=False):
        """
        Generate a report under the form of an Ordered dictionary
        """
        report = OrderedDict ()
        report["Query"] = "{}:{}-{}({})".format(
            self.q_id, self.q_start, self.q_end, self.q_orient)
        report["Subject"] = "{}:{}-{}({})".format(
            self.s_id, self.s_start, self.s_end, self.s_orient)

        if full:
            report["Identity"] = self.identity
            report["Evalue"] = self.evalue
            report["Bit Score"] = self.bscore
            report["Hit length"] = self.length
            report["Number of gap"] = self.gap
            report["Number of mismatch"] = self.mis

        return report

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _test_arg(self):
        assert 0 <= self.identity <= 100, "Identity value out of range [0:100]"
        assert self.length >= 1, "length value out of range [>= 1]"
        assert self.mis >= 0, "mis value out of range [>= 0]"
        assert self.gap >= 0, "gap value out of range [>= 0]"
        assert self.q_start >= 0, "q_start value out of range [>= 1]"
        assert self.q_end >= 0, "q_end value out of range [>= 1]"
        assert self.s_start >= 0, "s_start value out of range [>= 1]"
        assert self.s_end >= 0, "s_end value out of range [>= 1]"
        assert self.evalue >= 0, "evalue value out of range [>= 0]"
        assert self.bscore >= 0, "bscore value out of range [>= 0]"
