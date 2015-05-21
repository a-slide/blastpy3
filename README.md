# pyBlast 0.1

**Simple Python 2.7/3.3 wrapper for BLAST+**

**Creation : 2015/05/18**

**Last update : 2015/05/18**

## MakeBlastDB

...

## MakeBlastn

...

## BlastHit

...

## test_pyBlast

* Install pytest
* Run test with py.test-2.7  -v

Example of output if succesfull. Please note than some tests might fail due to the random sampling of DNA sequences.
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
ncbi blast+
python package pytest for tests

## Authors and Contact

Adrien Leger - 2015

* <adrien.leger@gmail.com> - <adrien.leger@inserm.fr> - <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089](http://www.atlantic-gene-therapies.fr/)