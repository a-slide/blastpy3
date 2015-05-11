# -*- coding: utf-8 -*-

"""
@package    pyBlast
@brief      Class to perform a blast alignment through Shell Command line
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library packages import
from os import remove, path
from multiprocessing import cpu_count
from subprocess import Popen, PIPE

# Local library packages
from BlastHit import BlastHit

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class MakeBlastn(object):
    """
    @class  MakeBlast
    @brief  Perform a blastn of a query against a blast database. If hits are found, a list of
    BlastHit objects is returned.
    Blast+ 2.8+ needs to be install and eventually added to the path.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, blast_exec="blastn", blastn_opt="", task="megablast", evalue=10,
        best_per_query_seq = False):
        """
        Initialize the object and index the reference genome if necessary
        @param blast_exec Path ot the blast executable. By Default blastn will be used
        @param blastn_opt Blastn command line options as a string
        @param task Type of blast to be performed ('blastn' 'blastn-short' 'dc-megablast'
        'megablast' 'rmblastn')
        @param evalue E Value cuttoff to retain alignments
        """
        # Creating object variables
        self.blast_exec = blast_exec
        self.blastn_opt = "{} -num_threads {} -task {} -evalue {} -outfmt 6 -dust no".format(
            blastn_opt, cpu_count(), task, evalue)
        self.best_per_query_seq = best_per_query_seq
        self.evalue = evalue

    def __str__(self):
        msg = "BLAST WRAPPER\n"
        msg += "Blastn path : {}\n".format(self.blast_exec)
        msg += "Options : {}\n".format(self.blastn_opt)
        msg += "Output best hit per query : {}\n".format("T" if self.best_per_query_seq else "F")
        return msg

    def __repr__(self):
        return "<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__ (self, query_path, db_path):
        """
        Blast query against a subject database and return a list of BlastHit object
        @param  query Path to a fasta file containing the query sequences
        @param blastdb Blast database object NewDB or ExistingDB
        @return A list of BlastHit objects if at least one hit was found
        """

        # Build the command line string
        cmd = "{} {} -query {} -db {}".format(self.blast_exec, self.blastn_opt, query_path, db_path)
        print ("CMD : " + cmd)

        # Execute the command line in the default shell
        proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate()

        # Verify the return code
        if proc.returncode == 1:
            msg = "An error occured during execution of following command :\n"
            msg += "COMMAND : {}\n".format(cmd)
            msg += "STDERR : {}\n".format(stderr)
            raise Exception (msg)

        # If no hit were found
        if not stdout:
            return None

        # Return a simple list if the unicity of hit per query is not required
        if not self.best_per_query_seq:
            hit_list = []
            for line in stdout.splitlines():
                h = line.split()
                hit_list.append(BlastHit(h[0], h[1] , h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11]))

            print ("\t{} hits found".format(len(hit_list)))
            return hit_list

        # The most complicated situation where only the best hit per query is returned
        # Create an intermediate dictionary to store all hits per query
        hit_dict = {}
        # Parse each result lines and create a list of BlastHit objects per query
        for i, line in enumerate(stdout.splitlines()):
            h = line.split()
            if h[0] in hit_dict:
                hit_dict[h[0]].append(BlastHit(h[0], h[1] , h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11]))
            else:
                hit_dict[h[0]] = [BlastHit(h[0], h[1] , h[2], h[3], h[4], h[5], h[6], h[7], h[8], h[9], h[10], h[11])]

        print ("\t{} hits found from {} query".format(i, len(hit_dict)))

        # Flatten the dictionary in a list keeping only the best alignment per query
        hit_list = []
        for query_hits in hit_dict.values():
            best_evalue = self.evalue
            best_hit = None
            for hit in query_hits:
                if hit.evalue <= best_evalue:
                    best_evalue = hit.evalue
                    best_hit = hit
            hit_list.append(best_hit)

        print ("\t{} hits retained".format(len(hit_list)))
        return hit_list
