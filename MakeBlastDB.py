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
from os import remove, path
from subprocess import Popen, PIPE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class MakeBlastDB(object):
    """
    @class MakeBlastDB
    @brief Wrapper for makeblastdb. Create a subject database from a fasta file
    Blast+ 2.8+ needs to be install and correctly added to the path.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, makeblastdb_exec="makeblastdb", makeblastdb_opt="", dbtype="nucl", input_type="fasta"):
        """
        Create a blastdb from a reference fastq file
        @param makeblastdb_exec Path of the makeblastdb executable by default "makeblastdb"
        @param makeblastdb_opt makeblastdb command line options as a string
        @param dbtype Molecule type of target db ('nucl', 'prot')
        @param input_type Type of the data specified in input_file ('asn1_bin', 'asn1_txt', 'blastdb', 'fasta')
        """
        # Creating object variables
        self.makeblastdb_exec = makeblastdb_exec
        self.makeblastdb_opt = "{} -dbtype {} -input_type {}".format(makeblastdb_opt, dbtype, input_type)

    def __str__(self):
        msg = "MAKEBLASTDB CLASS\n"
        msg += "Makeblastdb path : {}\n".format(self.makeblastdb_exec)
        msg += "Options : {}\n".format(self.makeblastdb_opt)
        return msg

    def __repr__(self):
        return "\n<Instance of {} from {} >\n".format(self.__class__.__name__, self.__module__)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self, ref_path, db_path="./out"):
        """
        Create a blastn database from ref_path using makeblastdb
        @param ref_path Path of the fasta file containing the reference sequence. Cannot be gziped
        @param db_path Outname for the blast db files basename.
        @return The absolute path of the database basename
        """

        # Build the command line
        cmd = "{} {} -in {} -out {}".format(self.makeblastdb_exec, self.makeblastdb_opt, ref_path, db_path)
        print ("CMD : " + cmd)

        # Run the command line without stdin and asking both stdout and stderr
        try:

            # Execute the command line in the default shell
            proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
            stdout, stderr = proc.communicate()

            # Verify the return code
            if proc.returncode == 1:
                msg = "An error occured during execution of following command :\n"
                msg += "COMMAND : {}\n".format(cmd)
                msg += "STDERR : {}\n".format(stderr)
                raise Exception (msg)

            # Verify the output
            if not stdout:
                raise Exception ("Error, no data received from standard output\n"+stderr)

            return path.abspath(db_path)

        # In case of exception during DB building remove the created files and reraise Exception
        except Exception as E:
            print("Remove database files")

            for ext in ["00.nhr", "nhr", "00.nin", "nin", "00.nsq", "nsq"]:
                f = "{}.{}".format(db_path, ext)
                if path.isfile (f):
                    remove (f)

            raise Exception (E.message+"Impossible to generate a valid database from the reference sequence")
