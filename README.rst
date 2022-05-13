Quick Start
===========


The objective of this tutorial is to reproduce part of the results presented in Baurens et al (2019) and Ahmed et al (2019), using respectively VCFHunter and TraceAncestor.

The outputs of these programs can then be used in the GeMo webapps.

Installation requirements
~~~~~~~~~~~~~~~~~~~~~~~~~

This tutorial is developed to run on Linux or Apple (MAC OS X) operating systems. There are no versions planned for Windows.

Software requirements:

-  Perl 5 for TraceAncestor
-  Python 3 for VCFHunter


Testing your Perl installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test that Perl 5 is installed, enter on the command line

::

    perl -version

Testing your Python installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To test that Python 3 is installed, enter on the command line

::

    python3 --version

Now, you can clone the repository, create a virtualenv and install several additional package using pip, as describe below.

::

   git clone https://github.com/gdroc/GeMo_tutorials.git
   python3 -m venv $PWD/venv
   source venv/bin/activate
   pip install numpy
   pip install matplotlib
   pip install scipy


Download Dataset
~~~~~~~~~~~~~~~~

For this tutorial, dataset that will be used by TraceAncestor or by VCFHunter have been deposited on Zenodo https://doi.org/10.5281/zenodo.6539270

To download this, you only need to launch the script download_dataset.pl without any parameter

::

   perl download_dataset.pl

This script create a new directory data

::

   data/
   ├── Ahmed_et_al_2019_color.txt
   ├── Ahmed_et_al_2019_individuals.txt
   ├── Ahmed_et_al_2019_origin.txt
   ├── Ahmed_et_al_2019.vcf
   ├── Baurens_et_al_2019_color.txt
   ├── Baurens_et_al_2019_individuals.txt
   ├── Baurens_et_al_2019_origin.txt
   ├── Baurens_et_al_2019_chromosome.txt
   └── Baurens_et_al_2019.vcf

These files are require for this tutorial to run VCFHunter or TraceAncestor


Input
~~~~~

- Baurens_et_al_2019_origin.txt : A two column file with individuals in the first column and group tag (i.e. origin) in the second column

=========== ======
individuals origin
=========== ======
P2      AA
T01     BB
T02     BB
T03     AA
T04     AA
T05     AA
T06     AA
T07     AA
T08     BB
=========== ======

- Baurens_et_al_2019.vcf : A vcf file with ancestral and admixed individuals

::

   grep #CHROM data/Baurens_et_al_2019.vcf
   #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ACC48-FPG	ACC48-FPN	ACC48-P_Ceylan	ACC48-Red_Yade	DYN163-Kunnan	DYN275-Pelipita	DYN359-Safet_Velchi	GP1	GP2	P1	P2	T01	T02	T03	T04	T05	T06	T07	T08	T10	T11

- Baurens_et_al_2019_individuals.txt : List of individuals to scan origin from vcf, as defined in the VCF headerline

- Baurens_et_al_2019_color.txt : A color file with 4 columns: col1=group and the three last column corresponded to RGB code.

===== ========== === === =
group name       r   g   b
===== ========== === === =
AA    acuminata  0   255 0
BB    balbisiana 255 0   0
===== ========== === === =

Run workflow using create_gemo_input.pl
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

   perl create_gemo_input.pl --help
   Parameters :
       -v, --vcf         A vcf file [required]
       -o, --origin       A two column file with individuals in the first column and group tag (i.e. origin) in the second column [Required]
       -i, --individuals  List of individuals to scan from vcf, as defined in the VCF headerline [Required]
       -m, --method      Permissible values: vcfhunter traceancestor (String). Default vcfhunter
       -c, --color       A color file with 4 columns: col1=group and the three last column corresponded to RGB code.
       -t, --threads     Number of threads
       -d, --dirout      Path to the output directory (Default method option name)
       -h, --help        display this help

For VCFHunter, you must use the dataset prefixed with Baurens_et_al.
::

    perl create_gemo_input.pl --vcf data/Baurens_et_al_2019.vcf --origin data/Baurens_et_al_2019_origin.txt --individuals data/Baurens_et_al_2019_individuals.txt --method vcfhunter --color data/Baurens_et_al_2019_color.txt --threads 4

For TraceAncestor, you must use the dataset prefixed with with Ahmed_et_al.
::

    perl create_gemo_input.pl --vcf data/Ahmed_et_al_2019.vcf --origin data/Ahmed_et_al_2019_origin.txt --individuals data/Ahmed_et_al_2019_individuals.txt --method traceancestor --color data/Ahmed_et_al_2019_color.tx




References
~~~~~~~~~~

-  `Summo, Marilyne. (2022). GeMo : a web-based platform for the visualization and curation of mosaic genomes [Data set]. Zenodo. <https://doi.org/10.5281/zenodo.6539270>`__
-  `Baurens,F.-C. et al.(2019) Recombination and Large Structural
   Variations Shape Interspecific Edible Bananas Genomes. Mol Biol Evol,
   36, 97–111. <https://doi.org/10.1093/molbev/msy199>`__
-  `Martin et al., 2020a. Martin G, Cardi C, Sarah G, Ricci S, Jenny C,
   Fondi E, Perrier X, Glaszmann J-C, D’Hont A, Yahiaoui N. 2020. Genome
   ancestry mosaics reveal multiple and cryptic contributors to
   cultivated banana. Plant J.
   102:1008–1025. <https://doi.org/10.1111/tpj.14683>`__
-  `Ahmed,D. et al. (2019) Genotyping by sequencing can reveal the
   complex mosaic genomes in gene pools resulting from reticulate
   evolution: a case study in diploid and polyploid citrus. Annals of
   Botany, 123, 1231–1251. <https://doi.org/10.1093/aob/mcz029>`__
   
