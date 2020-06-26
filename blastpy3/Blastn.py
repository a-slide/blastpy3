# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Std lib imports
from os import path, remove
from shutil import rmtree
from subprocess import Popen, PIPE
from tempfile import mkdtemp, mkstemp
from collections import OrderedDict

# Local lib packages
from blastpy3.BlastHit import BlastHit

# Third party lib imports
import pyfaidx

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Blastn(object):
    """
    Wrapper for Blastn. Create a subject database from a fasta file at initialisation and
    perform a blastn alignments by calling the object.
    Blast+ 2.8+ needs to be install and correctly added to the path.
    """
    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self,
        ref_path,
        makeblastdb_exec="",
        makeblastdb_opt="",
        verbose=False):
        """
        Create a blastdb from a reference fastq file
        * ref_path
            Path to the reference fasta file (non-gzipped) .
        * makeblastdb_exec
            Path of the makeblastdb executable. Default = "makeblastdb"
        * makeblastdb_opt
            makeblastdb command line options as a string. Default = ""
        """
        # Creating object variables
        self.ref_path = ref_path
        self.makeblastdb_exec = makeblastdb_exec if makeblastdb_exec else "makeblastdb"
        self.makeblastdb_opt = makeblastdb_opt

        self.db_dir = mkdtemp()
        self.db_path = path.join(self.db_dir, path.split(self.ref_path)[-1].partition('.')[0])
        self.verbose = verbose

        cmd = "{} {} -dbtype nucl -input_type fasta -in {} -out {}".format(
            self.makeblastdb_exec, self.makeblastdb_opt, self.ref_path, self.db_path)

        if self.verbose:
            print ("CREATE DATABASE: {}\n".format(cmd))

        # Run the command line without stdin and asking both stdout and stderr
        try:
            # Execute the command line in the default shell
            proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
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
        """Destructor to remove the database"""
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

    #~~~~~~~PROPERTIES~~~~~~~#

    @property
    def seq_len (self):
        seq_len = OrderedDict()
        with pyfaidx.Fasta(self.ref_path) as fa:
            for seq in fa:
                seq_len[seq.name] = len(seq)
        return seq_len

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def get_ref_seq (self, chrom, start=0, end=-1):
        try:
            with pyfaidx.Fasta(self.ref_path) as fa:
                return str(fa[chrom][start:end])
        except:
            return ""

    def align (self,
        query_path="",
        query_seq="",
        query_name="seq",
        blastn_exec="",
        blastn_opt="",
        task="blastn",
        evalue=1,
        max_hits=10,
        threads=2,
        sort_by="evalue"):
        """
        Blast query against a subject database and return a list of BlastHit object
        * query_path
            Path to a fasta file containing the query sequences (non-gzipped). Mandatory if query_seq is not set
        * query_seq
            Sequence string. Mandatory if query_path is not set
        * blastn_exec
            Path of the blast executable. By Default blastn will be used. Default = "blastn"
        * blastn_opt
            Blastn command line options as a string. Default = ""
        * task
            Type of blast to be performed among the following values: 'blastn', 'blastn-short', 'dc-megablast', 'megablast'. Default = "blastn"
        * evalue
            E Value cuttoff to retain alignments. Default = 1
        * max_hits
            Maximum number of alignments to keep per query sequence. Default = 10
        * sort_by
            key to sort hits on among the following values: 'identity', 'bscore', 'evalue'
        """

        blastn_exec = blastn_exec if blastn_exec else "blastn"

        # Requires at least a query sequence or a fasta query file
        if not query_path and not query_seq:
            raise Exception ("Please provide a query file or a query sequence")

        # If query sequence, write to disk as a fasta file
        elif query_seq:
            fd, query_path =  mkstemp(suffix=".fa")
            with open(query_path, "w") as fa_fp:
                fa_fp.write(">{}\n{}\n".format(query_name, query_seq))
            if self.verbose:
                print ("Write sequences to temporary fasta file: {}\n".format(query_path))


        cmd = "{} {} -num_threads {} -task {} -evalue {} -outfmt '6 std qseq' -dust no -query {} -db {}".format(
            blastn_exec, blastn_opt, threads, task, evalue, query_path, self.db_path)

        if self.verbose:
            print ("BLAST QUERY: {}\n".format(cmd))

        # Execute the command line in the default shell
        proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True, text=True)
        stdout, stderr = proc.communicate()

        # Remove temp fasta file if required
        if query_seq:
            remove(query_path)

        # Verify the return code
        if proc.returncode == 1:
            raise Exception ("COMMAND: {}\nSTDERR: {}\n".format(cmd, stderr.strip()))

        # Collect hits as Blasthits object
        hit_list = []
        for line in stdout.splitlines():
            hit_split = line.split()
            hit_list.append(BlastHit(*hit_split))

        if self.verbose:
            print ("\t{} hits found".format(len(hit_list)))

        if sort_by == "identity":
            hit_list.sort(key=lambda x: x.identity, reverse=True)
        elif sort_by == "bscore":
            hit_list.sort(key=lambda x: x.bscore, reverse=True)
        else:
            hit_list.sort(key=lambda x: x.evalue)

        # Remove extra hist if too many
        if max_hits and len(hit_list) > max_hits:
            hit_list = hit_list[0:max_hits]

        return hit_list

    def rm_db(self):
        if self.verbose:
            print ("Cleaning up blast DB files")
        rmtree(self.db_dir)
