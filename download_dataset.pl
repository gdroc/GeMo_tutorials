#!/usr/bin/perl
 

use warnings;
use strict;  

my $citrus_dataset = "https://zenodo.org/record/6542870/files/Ahmed_et_al_2019.zip";
my $banana_dataset = "https://zenodo.org/record/6542870/files/Baurens_et_al_2019.zip";

unless (-e "data/Baurens_et_al_2019.vcf") {
    system("mkdir data");
    system("curl --silent $banana_dataset -o data/Baurens_et_al_2019.zip;cd data; unzip Baurens_et_al_2019.zip");
    system("curl --silent $citrus_dataset -o data/Ahmed_et_al_2019.zip;cd data;unzip Ahmed_et_al_2019.zip");
}
