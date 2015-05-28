# -*- coding: utf-8 -*-

"""
@package    pyBlast
@brief      Class to create a blast database through Shell Command line
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2014
* <adrien.leger@gmail.com> <adrien.leger@inserm.fr> <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

# Standard library packages import
from os import path
from shutil import rmtree
from subprocess import Popen, PIPE
from tempfile import mkdtemp
from multiprocessing import cpu_count

# Local library packages
from BlastHit import BlastHit

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Blastn(object):
    """
    @class Blastn
    @brief Wrapper for Blastn. Create a subject database from a fasta file at initialisation and
    perform a blastn alignments by calling the object.
    Blast+ 2.8+ needs to be install and correctly added to the path.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, ref_path, makeblastdb_exec="", makeblastdb_opt="",
        dbtype="nucl", input_type="fasta"):
        """
        Create a blastdb from a reference fastq file
        @param ref_path Path to the reference fasta file (not gzipped). Mandatory
        @param makeblastdb_exec Path of the makeblastdb executable. Default = "makeblastdb"
        @param makeblastdb_opt makeblastdb command line options as a string. Default = ""
        @param dbtype Molecule type of target db ('nucl', 'prot'). Default = "nucl"
        @param input_type Type of the data specified in input_file ('asn1_bin', 'asn1_txt',
        'blastdb', 'fasta'). Default = "fasta"
        """
        # Creating object variables
        self.ref_path = ref_path
        self.makeblastdb_exec = makeblastdb_exec if makeblastdb_exec else "makeblastdb"
        self.makeblastdb_opt = makeblastdb_opt
        self.dbtype = dbtype
        self.input_type = input_type
        self.db_dir = mkdtemp()
        self.db_basename = self.ref_path.rpartition('/')[2].partition('.')[0]
        self.db_path = path.join(self.db_dir, self.db_basename)

        cmd = "{} {} -dbtype {} -input_type {} -in {} -out {}".format(
            self.makeblastdb_exec, self.makeblastdb_opt, self.dbtype, self.input_type,
            self.ref_path, self.db_path)

        #~ print ("CREATE DATABASE: {}\n".format(cmd))
        # Run the command line without stdin and asking both stdout and stderr
        try:
            # Execute the command line in the default shell
            proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            stdout, stderr = proc.communicate()

            # Verify the return code
            if proc.returncode == 1:
                raise Exception ("COMMAND: {}\nSTDERR: {}\n".format(cmd, stderr.strip()))

            # Verify the output
            if not stdout:
                raise Exception ("Error, no data received from standard output\n{}\n".format(stderr.strip()))

        except Exception as E:
            print (E)
            self.rm_db()
            self.db_dir = self.db_path = None

    # Enter and exit are defined to use the context manager "with"
    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        """Destructor to remove the database and unziped fasta files if needed"""
        if self.db_dir:
            self.rm_db()

    # Typical string methods
    def __str__(self):
        msg = "MAKEBLASTDB CLASS\tParameters list\n"
        # list all values in object dict in alphabetical order
        keylist = [key for key in self.__dict__.keys()]
        keylist.sort()
        for key in keylist:
            msg+="\t{}\t{}\n".format(key, self.__dict__[key])
        return (msg)

    def __repr__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__ (self, query_path, blastn_exec="", blastn_opt="", task="dc-megablast",
        evalue=1, best_query_hit=False):
        """
        Blast query against a subject database and return a list of BlastHit object
        @param  query_path Path to a fasta file containing the query sequences (not gzipped). Mandatory
        @param blastn_exec Path of the blast executable. By Default blastn will be used. Default = "blastn"
        @param blastn_opt Blastn command line options as a string. Default = ""
        @param task Type of blast to be performed ('blastn' 'blastn-short' 'dc-megablast'
        'megablast' 'rmblastn'). Default = "dc-megablast"
        @param evalue E Value cuttoff to retain alignments. Default = 1
        @param best_query_hit find and return only the best hit per query. Default = False
        @return A list of BlastHit objects if at least one hit was found
        """

        blastn_exec = blastn_exec if blastn_exec else "blastn"

        cmd = "{} {} -num_threads {} -task {} -evalue {} -outfmt \"6 std qseq\" -dust no -query {} -db {}".format(
            blastn_exec, blastn_opt, cpu_count(), task, evalue, query_path, self.db_path)

        #~ print ("MAKE BLAST: {}\n".format(cmd))
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
        if not best_query_hit:
            hit_list = []
            for line in stdout.splitlines():
                # split the line in a list and send the expanded list to blasthit
                hit_split = line.split()
                hit_list.append(BlastHit(*hit_split))

            #~ print ("\t{} hits found".format(len(hit_list)))
            return hit_list

        # The most complicated situation where only the best hit per query is returned
        # Create an intermediate dictionary to store all hits per query
        hit_dict = {}
        # Parse each result lines and create a list of BlastHit objects per query
        for i, line in enumerate (stdout.splitlines(), 1):
            # split the line in a list and send the expanded list to blasthit
            hit_split = line.split()
            query_id = hit_split[0]
            if query_id in hit_dict:
                hit_dict[query_id].append(BlastHit(*hit_split))
            else:
                hit_dict[query_id] = [BlastHit(*hit_split)]

        #~ print ("\t{} hits found from {} query".format(i, len(hit_dict)))

        # Flatten the dictionary in a list keeping only the best alignment per query
        hit_list = []
        for query_hits in hit_dict.values():
            best_evalue = 100 # Start with a very high evalue
            best_hit = None # and no hit
            for hit in query_hits:
                if hit.evalue <= best_evalue:
                    best_evalue = hit.evalue
                    best_hit = hit
            hit_list.append(best_hit)

        #~ print ("\t{} hits retained".format(len(hit_list)))
        return hit_list

    def rm_db(self):
        print (" * Cleaning up blast DB files for \"{}\"\n".format(self.db_basename))
        rmtree(self.db_dir)
