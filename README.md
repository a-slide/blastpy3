# pyBlast 0.1

See [GitHub page]( http://a-slide.github.io/pyBlast)

**Simple and lightweight Python 2.7 wrapper module for NCBI BLAST+**

**Creation : 2015/05/18**

**Last update : 2015/05/22**

## BlastHit

Python object representing a hit found by blastn. The object contains the following public fields:

* id: Auto incremented unique identifier [INT]
* q_id: Query sequence name [STR]
* s_id: Subject sequence name [STR]
* identity: % of identity in the hit [FLOAT 0:100]
* length: length of the hit [INT >=0]
* mis: Number of mismatch in the hit [INT >=0]
* gap: Number of gap in the hit [INT >=0]
* q_start: Hit start position of the query sequence [INT >=0]
* q_end: Hit end position of the query sequence [INT >=0]
* s_start: Hit start position of the subject sequence [INT >=0]
* s_end: Hit end position of the subject sequence [INT >=0]
* evalue: E value of the alignment [FLOAT >=0]
* bscore: Bit score of the alignment[FLOAT >=0]
* q_seq: Sequence of the query aligned on the subject sequence [STR]
* q_orient: Orientation of the query sequence [+ or -]
* s_orient: Orientation of the subject sequence [+ or -]

The validity of numeric value is checked upon instantiation. Invalid values will raise assertion errors.

## Blastn

This class contain the wrapper for Blastn and require the installation of ncbi Blast+ 2.2.28+.

### Setup Blastn object: Create subject database

Upon instantiation, a database is created from the user-provided subject sequence. Database files are created in a temporary directory.
The following parameters can be customized at Blastn objects instantiation
* ref_path: Path to the reference fasta file (not gzipped). Mandatory
* makeblastdb_exec: Path of the makeblastdb executable. Default = "makeblastdb"
* makeblastdb_opt: makeblastdb command line options as a string. Default = ""
* dbtype: Molecule type of target db ('nucl', 'prot'). Default = "nucl"
* input_type: type of the data specified in input_file ('asn1_bin', 'asn1_txt', 'blastdb', 'fasta'). Default = "fasta"

To ensure a proper database files deletion at the end of the execution it is possible to call the object using the `with` statement.
Alternatively you can call the `rm_db` method at the end of the Blastn usage.

**Code**
```
with Blastn(ref_path="./subject.fa") as blastn:
    print (blastn)
```
**Output**
```
CREATE DATABASE: makeblastdb  -dbtype nucl -input_type fasta -in subject.fa -out /tmp/tmpihszgZ/subject

MAKEBLASTDB CLASS	Parameters list
	db_basename	subject
	db_dir	/tmp/tmpihszgZ
	db_path	/tmp/tmpihszgZ/subject
	dbtype	nucl
	input_type	fasta
	makeblastdb_exec	makeblastdb
	makeblastdb_opt	
	ref_path	subject.fa

Cleaning up blast DB files for "subject"
```

### Calling Blastn object: Perform Blastn and return a list of hits

The Blastn object can then be directly called, as many time as desired, with query fasta files, that can contain several sequences.
The method will automatically call blastn in a multiprocessing fashion, using as many threads as possible.
The following parameters can be customized at Blastn objects calling:

* query_path: Path to a fasta file containing the query sequences (not gzipped). Mandatory
* blast_exec: Path of the blast executable. By Default blastn will be used. Default = "blastn"
* blastn_opt: Blastn command line options as a string. Default = ""
* task: Type of blast to be performed ('blastn' 'blastn-short' 'dc-megablast' 'megablast' 'rmblastn'). Default = "dc-megablast"
* evalue: E Value cuttoff to retain alignments. Default = 1
* best_query_hit: find and return only the best hit per query. Default = False

A list containing 1 BlastHit object for each query hit found in the subject will be returned, except if not hit were found in which situation 'None' will be returned. 
If the best_query_hit flag was set to True, Only the best hit per query sequence from the query file will be returned.

**Code**
```
with Blastn(ref_path="./subject.fa") as blastn:
    hit_list = blastn(query_path="./query.fa")
    for hit in hit_list:
        print (hit)
```
**Output**
```
CREATE DATABASE: makeblastdb  -dbtype nucl -input_type fasta -in ./subject.fa -out /tmp/tmp1ZBlfT/subject

MAKE BLAST: blastn  -num_threads 4 -task dc-megablast -evalue 1 -outfmt "6 std qseq" -dust no -query ./query.fa -db /tmp/tmp1ZBlfT/subject

	2 hits found
HIT 0	Query	query1:0-48(+)
	Subject	subject:19-67(+)
	Lenght : 48	Identity : 100.0%	Evalue : 2e-23	Bit score : 87.8
	Aligned query seq : GCATGCTCGATCAGTAGCTCTCAGTACGCATACGCTAGCATCACGACT

HIT 1	Query	query2:0-48(+)
	Subject	subject:89-137(+)
	Lenght : 48	Identity : 100.0%	Evalue : 2e-23	Bit score : 87.8
	Aligned query seq : CGCATCGACTCGATCTGATCAGCTCACAGTCAGCATCAGCTACGATCA

Cleaning up blast DB files for "subject"
```


## Testing pyBlast module

The module can be easily tested thanks to pytest

* Install pytest with pip `pip instal pytest`
* Run test with py.test-2.7  -v

Example of output if successful. Please note than some tests might fail due to the random sampling of DNA sequences, and uncertainties of Blastn algorithm.
```
========================================== test session starts ===========================================
platform linux2 -- Python 2.7.5 -- py-1.4.27 -- pytest-2.7.0 -- /usr/bin/python
rootdir: /home/adrien/Programming/Python/pyBlast, inifile: 
collected 21 items 

test_pyBlast.py::test_BlastHit[4.16866907958-57-98-69-88-12-100-43-1.40452897105-47.3666242716] PASSED
test_pyBlast.py::test_BlastHit[-1-7-10-20-73-54-25-45-98.7921480151-45.2397166228] xfail
test_pyBlast.py::test_BlastHit[8.92741377413--1-100-36-34-33-14-71-18.8547135761-97.6604693294] xfail
test_pyBlast.py::test_BlastHit[10.5987790458-46--1-45-78-81-86-86-73.8740266727-56.887410005] xfail
test_pyBlast.py::test_BlastHit[66.8213911219-62-48--1-91-10-60-20-88.7850139735-81.7901609219] xfail
test_pyBlast.py::test_BlastHit[86.6626174287-29-83-34--1-53-57-68-17.9799756069-7.83036609495] xfail
test_pyBlast.py::test_BlastHit[5.23985331666-43-85-33-7--1-14-3-74.2130782704-88.9289495285] xfail
test_pyBlast.py::test_BlastHit[75.6935977321-8-78-68-10-39--1-74-44.1447867052-22.5203082483] xfail
test_pyBlast.py::test_BlastHit[39.8692596061-60-5-49-77-9-31--1-2.59963139531-46.3133849683] xfail
test_pyBlast.py::test_BlastHit[15.7192632366-24-92-1-64-82-83-90--1-75.5540618409] xfail
test_pyBlast.py::test_BlastHit[18.6627439886-34-57-60-5-45-26-40-77.7840842678--1] xfail
test_pyBlast.py::test_Blastn[blastn-Queries from Subject] PASSED
test_pyBlast.py::test_Blastn[blastn-Random queries] xfail
test_pyBlast.py::test_Blastn[blastn-short-Queries from Subject] PASSED
test_pyBlast.py::test_Blastn[blastn-short-Random queries] xfail
test_pyBlast.py::test_Blastn[dc-megablast-Queries from Subject] PASSED
test_pyBlast.py::test_Blastn[dc-megablast-Random queries] xfail
test_pyBlast.py::test_Blastn[megablast-Queries from Subject] PASSED
test_pyBlast.py::test_Blastn[megablast-Random queries] xfail
test_pyBlast.py::test_Blastn[rmblastn-Queries from Subject] PASSED
test_pyBlast.py::test_Blastn[rmblastn-Random queries] xfail

================================== 6 passed, 15 xfailed in 5.91 seconds ==================================
```

## Dependencies

* [ncbi Blast+ 2.2.28+](http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [python package pytest](http://pytest.org/latest/): `pip instal pytest`

## Authors and Contact

Adrien Leger - 2015

* <adrien.leger@gmail.com> - <adrien.leger@inserm.fr> - <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089](http://www.atlantic-gene-therapies.fr/)
