
README File on acdFIMD.py

System requirements:

CentOS 7.6.1810
Python-3.7.1
Python package(latest version): pysam, numpy, scipy.stats
No specific non-standard hardware needed

Installation instruction:

1. Download and copy acdFMID.py to a local directory
2. chmod 700 acdFMID.py

Typical installation time:
Within 5 seconds

Usage:

python3 acdFMID.py -h
(Usage of the software can be expected if the installation process is a success)

usage: acdFMID.py [-h] [-v INVCF] [-b INBAM] [-l LIBID] [-q LIB2QC] [-i INSSM]
                  [-c RINSC] [-d OUTDIR]

=============================================================================================================

Version: 1.1.20201106   Date: Dec 11 CST 2020
Author: Jing Ruilin     Contact: see data share disclaimer.

Description

   ACD & FMID Filter On SNV Only VCF


=============================================================================================================

optional arguments:
  -h, --help            show this help message and exit
  -v INVCF, --invcf INVCF
                        Input VCF File
  -b INBAM, --inbam INBAM
                        Input BAM File
  -l LIBID, --libid LIBID
                        Sample Library ID (In accord with VCF FORMAT)
  -q LIB2QC, --lib2qc LIB2QC
                        LibID to QC Info File
  -i INSSM, --inssm INSSM
                        InsertSize Smoothing Instead of Filtering
  -c RINSC, --rinsc RINSC
                        Ref InsertSize Correction By Alt
  -d OUTDIR, --outdir OUTDIR
                        General Output Directory

Demo run:

Demo command:
(Demo input data have be modified as chr4:1799442-1806167)

/XXXX/software/python/python-3.7.1/bin/python3.7 /XXXX/projects/acd_fmid/acdFMID.py 
-l F300001661_L01_612_2030579P1NP
-v /XXXX/projects/acd_fmid/testdata/F300001661_L01_612_2030579P1NP.snv.vcf.gz
-b /XXXX/projects/acd_fmid/testdata/F300001661_L01_612_2030579P1NP.bam
-q /XXXX/projects/acd_fmid/testdata/libid2qc.xls
-d /XXXX/projects/acd_fmid/demoResults

Expected output on screen (STDOUT):

Feature Activation Status:

# InsertSize Smoothing: on
# Ref InsertSize Correction By Alt: on

Typical run time:
Within 1 minutes

Expected output files:
/XXXX/projects/acd_fmid/demoResults/F300001661_L01_612_2030579P1NP.acd_fmid.xls
/XXXX/projects/acd_fmid/demoResults/F300001661_L01_612_2030579P1NP.addVariantInsInfo.bam (used for debug)
