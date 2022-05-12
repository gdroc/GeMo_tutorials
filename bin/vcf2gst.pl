#!/usr/bin/perl

# Calculer les GST pour les individus dialéliques

use warnings;
use strict; 
use Getopt::Long; 
use FindBin;
use Pod::Usage;

my $vcf;
my $ancestors;
my $DM;
my $nbMark;
my $depth = 5;
my $output = "GSTmatrice.txt";
my $help = "";
my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help
Parameters :
    --vcf       vcf containing the ancestors [Required]
    --ancestor  ancestor file [Required]
    --depth     minimal depth for a snp to be used in the analysis (Default 5)
    --output    output file name (Default GSTmatrice.txt) 
    --help
/;
Getopt::Long::Configure ('bundling');
GetOptions (
	'v|vcf=s'      => \$vcf,
	'a|ancestor=s' => \$ancestors, 
	'd|depth=i'    => \$depth,
	'o|output=s'   => \$output,
 	'h|help!'      => \$help,
) or pod2usage(-message => $usage);

if ($help) { pod2usage(-message => $usage); }


## check to presence of necessary parameters
if (!$vcf){
    warn $usage;
	warn "\nWarn :: --vcf is empty. Please specify a vcf file\n\n";
    exit 0; 
}   
if (!$ancestors){
    warn $usage;
	warn "\nWarn :: --ancestor is empty. Please specify a ancestor file\n\n";
    exit 0; 
} 

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
my %ancetres; 
my @individusAncestraux;
my %list_ancestor; 

open(IN,$ancestors)  or die ("Error: ancestors wont open\n");  
while(<IN>){
    chomp;
    my ($accession,$ancestor) = (split(/\t/,$_));
    if ($accession ne "") {
        $ancetres{$accession} = $ancestor;
        push(@individusAncestraux, $accession);
        $list_ancestor{$ancestor} = 1; 
    }
}
close IN;
my @especesAncestrales = ( keys %list_ancestor ); 
 
my @valBlock;
my @individus;
my %heterozygotie; 
my $chr;
my $pos;
my %refalt;

my %Nref;
my $nbOfN;
my $nbOfInd;

my %FormatHash;
my @format;
my $cpt = 0;
open(F1,"$vcf") or die ("Error: vcf wont open\n"); #opening of the ancestor matrix
while (my $li = <F1>){ # for each marker
	chomp($li);
	if ($li =~ m/^#.+$/){
		if ($li =~ m/(#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT.+)$/){
			@individus = split(/\s+/,$li); # List of samples 
		}
	}
	else{
		$nbOfInd = 0;
		$nbOfN = 0;
		@valBlock = split(/\t/,$li);
		$chr = $valBlock[0]; # chromosome
		$pos = $valBlock[1]; # position
		$refalt{$chr}{$pos}{"REF"} = $valBlock[3]; # ref
		$refalt{$chr}{$pos}{"ALT"} = $valBlock[4]; # alt
		@format = split(":",$valBlock[8]);
		for (my $ind = 9; $ind < scalar(@individus); $ind++){ 
			my@formatIND = split(":", $valBlock[$ind]);
			for (my$f = 0; $f < scalar(@format); $f++){
                if ($ancetres{$individus[$ind]}){
                    $FormatHash{$individus[$ind]}{$format[$f]} = $formatIND[$f]; 
                }
				
			}
		}
		for my $formInd(sort keys %FormatHash){ 
			$nbOfInd++;
			my $het;
			if ($FormatHash{$formInd}{"GT"} =~ m/^(1|0)\/(1|0)/){
                if($FormatHash{$formInd}{"DP"} >= $depth){
                    if ($FormatHash{$formInd}{"GT"} eq "1/1") {$het = 1;}
                    if ($FormatHash{$formInd}{"GT"} eq "1/0") {$het = 0.5;}
                    if ($FormatHash{$formInd}{"GT"} eq "0/1") {$het = 0.5;}
                    if ($FormatHash{$formInd}{"GT"} eq "0/0") {$het = 0;} 
				}
				else{
					$het = "NA";
					$nbOfN++;
				}
			}
			else{
				$het = "NA";
				$nbOfN++;
			}  
            push @{$heterozygotie{$ancetres{$formInd}}{$chr}{$pos}} , $het;
        }
		if($nbOfInd == 0){
			print"Error: Check your ancestors file. No match found.\n";
			exit;
		}
		else{
			$Nref{$chr}{$pos} = $nbOfN / $nbOfInd;
		}
	}
}
close F1;
 
 
### calcul Fréquences pour chaque ancètre

my %freqParAncetre; #FM
my $valueofhet;
my $nbofindNotNA;
for my $ancetre (sort keys %heterozygotie){
	for my $chr(sort keys %{$heterozygotie{$ancetre}}){
		for my $pos(sort {$a<=>$b} keys %{$heterozygotie{$ancetre}{$chr}}){
			$valueofhet = 0;
			$nbofindNotNA = 0;
            foreach my $value (@{$heterozygotie{$ancetre}{$chr}{$pos}}) {
                if ($value ne "NA"){
                    $valueofhet = $valueofhet + $value;
                    $nbofindNotNA++;
                }
            }
			if($nbofindNotNA != 0){
				$freqParAncetre{$chr}{$pos}{$ancetre} = $valueofhet / $nbofindNotNA;
			}
			else{
				$freqParAncetre{$chr}{$pos}{$ancetre} = "NA";
			}
		}		
	}
}


### calcul Fréquences de tous les autres ancètres que celui auquel on s'interresse

my %freqSfAncetre; #F-M
my $calculSfAnc;
for my $chr (sort keys %freqParAncetre){
	for my $pos(sort {$a<=>$b} keys %{$freqParAncetre{$chr}}){	
		foreach my $anc(@especesAncestrales){
			my $booNa = 0;
			$calculSfAnc = 0;
			for my $a (sort keys %{$freqParAncetre{$chr}{$pos}}){
				if($anc ne $a){
					if($freqParAncetre{$chr}{$pos}{$a} ne "NA"){
						$calculSfAnc = $calculSfAnc + $freqParAncetre{$chr}{$pos}{$a};
					}
					else{
						$booNa = 1;
					}
				}
			}
			if($booNa == 0){
			 	$freqSfAncetre{$chr}{$pos}{$anc} = $calculSfAnc/(scalar(@especesAncestrales) - 1);
			}
			else{
				$freqSfAncetre{$chr}{$pos}{$anc} = "NA";
			}
		}
	}
}



### calcul des fréquences totales + H
# HM = 1 - (FM² + (1-FM²))
### calcul diversité de tous les autres ancètres que celui auquel on s'interresse H-M

my%freqTot; #FMtot
my%HSfAncetre; #H-M
my%H; #HM
for my$chr (sort keys %freqParAncetre){
	for my$pos(sort {$a<=>$b} keys %{$freqParAncetre{$chr}}){
		foreach my$anc(@especesAncestrales){
			if ($freqParAncetre{$chr}{$pos}{$anc} ne "NA" ){
				$H{$chr}{$pos}{$anc} = 1-(($freqParAncetre{$chr}{$pos}{$anc}*$freqParAncetre{$chr}{$pos}{$anc}) + (1-$freqParAncetre{$chr}{$pos}{$anc})*(1-$freqParAncetre{$chr}{$pos}{$anc}));
			}
			else{
				$H{$chr}{$pos}{$anc} = "NA";
			}
			if ($freqSfAncetre{$chr}{$pos}{$anc} ne "NA"){
				$HSfAncetre{$chr}{$pos}{$anc} = 1-(($freqSfAncetre{$chr}{$pos}{$anc}*$freqSfAncetre{$chr}{$pos}{$anc}) + (1-$freqSfAncetre{$chr}{$pos}{$anc})*(1-$freqSfAncetre{$chr}{$pos}{$anc}));
			}
			else{
				$HSfAncetre{$chr}{$pos}{$anc} = "NA";
			}
			if ($freqParAncetre{$chr}{$pos}{$anc} ne "NA" and $freqSfAncetre{$chr}{$pos}{$anc} ne "NA"){
				$freqTot{$chr}{$pos}{$anc} = ($freqParAncetre{$chr}{$pos}{$anc} + $freqSfAncetre{$chr}{$pos}{$anc})/2;
			}
			else{
				$freqTot{$chr}{$pos}{$anc} = "NA";
			}
		}
	}
}

## Calcul diersité totale HMTot et GST
# GSTM = (HMtot - (HM + H-M)/2) / HMtot
my%Htot;
my%GST;
for my$chr (sort keys %H){
	for my$pos(sort {$a<=>$b} keys %{$H{$chr}}){	
		foreach my$anc(@especesAncestrales){
			if (($H{$chr}{$pos}{$anc} ne "NA") and ($HSfAncetre{$chr}{$pos}{$anc} ne "NA")){
				$Htot{$chr}{$pos}{$anc} = 1 - (($freqTot{$chr}{$pos}{$anc}*$freqTot{$chr}{$pos}{$anc}) + (1-$freqTot{$chr}{$pos}{$anc})*(1-$freqTot{$chr}{$pos}{$anc}));
				if ($Htot{$chr}{$pos}{$anc} != 0){
					$GST{$chr}{$pos}{$anc} = ($Htot{$chr}{$pos}{$anc} - ($H{$chr}{$pos}{$anc} + $HSfAncetre{$chr}{$pos}{$anc})/2) / $Htot{$chr}{$pos}{$anc};

				}
				else{
					$Htot{$chr}{$pos}{$anc} = "NA";
					$GST{$chr}{$pos}{$anc} = "NA";
				}
			}
			else{
				$Htot{$chr}{$pos}{$anc} = "NA";
				$GST{$chr}{$pos}{$anc} = "NA";
			}
		}
	}
}


##Ecriture de la matrice des GST
open(F2,">$output") or die ("Error: output wont open\n");
print F2 join("\t","#CHROM","POS","REF","ALT","\%Nref");
for(my$i = 0; $i < scalar(@especesAncestrales); $i++){
	my $anc = $especesAncestrales[$i];
	print F2 "\tGST$anc";
}
for(my$i = 0; $i < scalar(@especesAncestrales); $i++){
	my$anc = $especesAncestrales[$i];
	print F2 "\tF$anc";
}
print F2 "\n";
for my $chr (sort keys %GST){
	for my $pos(sort {$a<=>$b} keys %{$GST{$chr}}){
		my $boo = 0;
		for(my$i = 0; $i < scalar(@especesAncestrales); $i++){
			my$anc = $especesAncestrales[$i];
			if($GST{$chr}{$pos}{$anc} eq "NA" || $freqParAncetre{$chr}{$pos}{$anc} eq "NA"){
				$boo = 1; 
			}
		}
		if($boo == 0){
			print F2 join("\t",$chr,$pos,$refalt{$chr}{$pos}{'REF'},$refalt{$chr}{$pos}{'ALT'},$Nref{$chr}{$pos});
			for(my$i = 0; $i < scalar(@especesAncestrales); $i++){
				my $anc = $especesAncestrales[$i];
				print F2 "\t$GST{$chr}{$pos}{$anc}";
			}
			for(my$i = 0; $i < scalar(@especesAncestrales); $i++){
				my $anc = $especesAncestrales[$i];
				print F2 "\t$freqParAncetre{$chr}{$pos}{$anc}";
			}
			print F2 "\n";
		}
	}
}

close F2;

