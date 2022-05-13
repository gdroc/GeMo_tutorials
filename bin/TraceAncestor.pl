#!/usr/bin/perl

# Date : 22/10/21
# Author : Aurore Comte

use warnings;
use strict;

use Getopt::Long; 
use FindBin;
use Pod::Usage;

# Initialization of files and variables that can be chosen by the user to launch the script.
my $matrix;
my $vcf;  
my $window_size = 10;
my $window_cut = 100;
my $LOD = 3;
my $freq = 0.99;
my $individuals;
my $help = ""; 
my $merge; 
my $dirout = "result";
my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help
Parameters :
    --matrix     Diagnosis matrix [Required]
    --vcf       vcf of the hybrid population 
    --individuals    A two column file with individuals to scan for origin (same as defined in the VCF headerline) in the first column and the ploidy in the second column [Required]
    --window    number of markers by window (Default 10)
    --lod       LOD value to conclude for one hypothesis (Default 3)
    --freq      theoretical frequency used to calcul the LOD (Default 0.99)
    --cut       number of K bases in one window (Default 100) 
    --dirout    Directory output (Default result)
    --help      display this help
/;
Getopt::Long::Configure ('bundling');
GetOptions (
	'm|matrix=s'     => \$matrix, ## reference matrix
	'v|vcf=s'       => \$vcf, ## vcf of hybrids. Hybrids should not have space in their names and no special characters. 
	'i|individuals=s' => \$individuals, 
	'w|window=i'    =>\$window_size, ## number of ancestral markers by windows on a chromosome
	'l|lod=i'       =>\$LOD, ## LOD value to conclude for one hypothesis or an other
	't|threshold=f' =>\$freq, ## Theoretical frequency used to calcul the LOD
	'k|cut=i'       => \$window_cut, ## number of K bases in one window 
    'd|dirout=s'      => \$dirout,
	'h|help!'       => \$help
)
or pod2usage(-message => $usage);

if ($help) { pod2usage(-message => $usage); }
if (!$matrix){
    warn $usage;
	warn "\nWarn :: --matrix is empty. Please specify a Diagnisis matrix\n\n";
    exit 0; 
}  
 
if (!$vcf){
    warn $usage;
	warn "\nWarn :: --vcf is empty. Please specify a vcf file\n\n";
    exit 0; 
}
if (!$individuals){
    warn $usage;
	warn "\nWarn :: --individuals is empty. Please specify a file with individuals to scan for origin\n\n";
    exit 0; 
}  
#------------------------------------------------------------------------------------------------------------------------------------------------------
# STEP 0 : If we work on WGS --> we redo an ancestor matrix with only snp that are also in the vcf file.

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------
my %hybrid;  

open(IN,$individuals)  or die ("Error: ancestors wont open\n");  
while(<IN>){
    chomp;
    my ($accession,$ploidy) = (split(/\t/,$_));
    if ($accession ne "") {
        $hybrid{$accession} = $ploidy; 
    }
}
close IN;  
 
system("mkdir $dirout") unless -e $dirout;

#group	name	hex
#C	Citron	#fff34d
#P	Pummello	#2c6ee8
#M	Mandarin	#ff0000
#Mic	Papeda	#2aa72e
#undefined	undefined	#24a6dd
 

#------------------------------------------------------------------------------------------------------------------------------------------------------

# STEP 1: creation of the not overlapping windows of $window_size (parameter chose by the user) markers

my%bloc_position; # key 1 = ancestor / key 2 = chromosome / key 3 = position of the end of the window => value = every marker positions in one window separated by a -
my%endChro; # key = chromosome => value = position of the last marker of every chromosome. We will assume this is the length of the chromosome.
my@ref; # list of the name of each ancestral genome.
my%allele_refalt; # key 1 = ancestor / key 2 = chromosome / key3 = position of the marker => value = base of the ancestral allele
my%catPos; # key 1 = ancestor / key 2 = chromosome => value = concatenation of every position of a windows separated by a -
my%hashRef; # marker counter by window
my@chromosomList; # list of chromosomes in the ancestor matrix
open(MATRIX,$matrix) or die ("Error: matrix of ancestors wont open\n"); #opening of the ancestor matrix

while (my $line = <MATRIX>){ # for each marker
	chomp($line);
	# if line is not the header
	if ($line !~ m/^ancestor/){
		my ($ancetre, $chromosome, $position, $allele) = (split(/\t/, $line)); 
		# list of last position for every chromosome
		$endChro{$chromosome} = $position;
		if (!grep { $_ eq $ancetre } @ref){
			push (@ref, $ancetre);
		}
		# list of chromosomes
		if (!grep { $_ eq $chromosome } @chromosomList){
			push (@chromosomList, $chromosome);
		}
		# We make windows of 10 consecutives markers -> %bloc_position
		if (!defined($catPos{$ancetre}{$chromosome})){
			$catPos{$ancetre}{$chromosome} = $position; #initialization of concatenation of positions at the begining of each chromosome.
			$hashRef{$ancetre}{$chromosome} = 1; #initialization of the marker counter
		}
		if ($hashRef{$ancetre}{$chromosome} != $window_size-1){ # if the counter is less than the wanted number of markers by windows
			if ($catPos{$ancetre}{$chromosome} eq ""){
				$catPos{$ancetre}{$chromosome} = $position;
				$hashRef{$ancetre}{$chromosome} = 1;
			}
			elsif ($catPos{$ancetre}{$chromosome} ne $position){
				$catPos{$ancetre}{$chromosome} = $catPos{$ancetre}{$chromosome}."-".$position;
				$hashRef{$ancetre}{$chromosome}++;
			}
		}
		elsif($hashRef{$ancetre}{$chromosome} == $window_size-1){ # if the counter equal the wanted number of markers by windows
			$bloc_position{$ancetre}{$chromosome}{$position} = $catPos{$ancetre}{$chromosome}."-".$position;
			$hashRef{$ancetre}{$chromosome} = 0; # reinitialization
			$catPos{$ancetre}{$chromosome} = "";
		}
		$allele_refalt{$ancetre}{$chromosome}{$position} = $allele;
	}
}
close MATRIX;
# We are numbering each windows in a hash %blocs.			
my %blocs;
my %blocspos;
my$end;
for my $parent (sort keys %bloc_position) {
    for my $chro (sort keys %{$bloc_position{$parent}}) {
        my $boolean = 0; #if boo==0, it's the begining of a chromosome
		my $num = 0;
		for my $pos (sort {$a<=>$b} keys %{$bloc_position{$parent}{$chro}}) {
			$num++;
			$blocs{$parent}{$chro}{$num} = $bloc_position{$parent}{$chro}{$pos};
			my@marqueur = split("-",$blocs{$parent}{$chro}{$num});
			if ($boolean == 0){
				$blocspos{$parent}{$chro}{$num}= "1-".$marqueur[scalar(@marqueur)-1];
			}
			else{
				$blocspos{$parent}{$chro}{$num}=$end."-".$marqueur[scalar(@marqueur)-1];
			}
			$end = $marqueur[scalar(@marqueur)-1];
			$boolean = 1;
		}
	}
}

#--------------------------------------------------------------------------------------------------------

# STEP 2: Calculation of the frequency of specifics reads for each ancestors by window.

open(VCF,$vcf) or die ("Error: can't open the vcf file\n");
my@individus; # list of samples. Samples names begin at $individu[9]
my@valBlock; # On line of the vcf
my%nbAncestralOrNot; # key 1: ancestor / key 2 : hybrid / key 3 : chromosome / key 4: window number / key 5 : "NOT" or "ANCESTRAL" / value : If "NOT" = sums of non-ancestral reads. If "ANCESTRAL" = sum of ancestral reads
my$chr;
my$pos="";
my$numbloc; 
my$AlleleRef; # Ancestral allele for a given position and a given ancestor
my$Ancestral = "ANCESTRAL";
my$Not = "NOT";
my$base1; # first base for one position in the vcf
my$base2; # second base for one position in the vcf (can be multiple)
my%frequences; # Frequency of ancestral alleles by marker
my%NM; # hash counting the number of marker really present in the vcf for each window (could be usefull later).
my%realEndChro; #real length of each chromosome. We will use these values if they are available.
my$booEndChro = 0; # if = 0 we use %endChro else we use %realEndChro as end of chromosome

my%FormatHash;
my@format;

while (my $line = <VCF>){
	chomp($line);
	if ($line =~ m/^#.+$/){
		if ($line =~ m/(#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT.+)$/){
			@individus = split(/\s/,$line); # List of samples
		}
		if ($line =~ m/##contig=<ID=(\S+)\,length=(\d+)\>/){
			$realEndChro{$1} = $2; # length ($2) for every chromosomes ($1)
			$booEndChro = 1;
		}
	}
	else{ # for each position / line in the vcf
		@valBlock = split("\t",$line);
		$chr = $valBlock[0]; # chromosome
		if(!$chr){
			next;
		}
		if(!grep { $_ eq $chr } @chromosomList){ # we work only on the chromosomes that are in the matrix of ancestors.
			next;
		};
		$pos = $valBlock[1]; # position
		$base1 = $valBlock[3];
		$base2 = $valBlock[4];
		@format = split(":",$valBlock[8]);
		for (my$ind = 9; $ind < scalar(@individus); $ind++){ 
			my@formatIND = split(":", $valBlock[$ind]);
			for (my$f = 0; $f < scalar(@format); $f++){ 
				$FormatHash{$individus[$ind]}{$format[$f]} = $formatIND[$f];
			}
		}
		my@bases2; # list for tue multiple alternative base in the vcf alternative base
		if ($base2 =~ /[\w|\*\,]+/){ # if base2 is multiple, we consider the first one (the more frequent)
			@bases2 = split(",", $base2);
		}
		for my $parent (sort keys %blocs){
			foreach my $num(sort {$a<=>$b} keys(%{$blocs{$parent}{$chr}})){
				my @allpos = split("-",$blocs{$parent}{$chr}{$num}); # list of the markers in this window
				if (grep { $_ == $pos } @allpos){ # We check that the snp is in the reference marker list.
					$NM{$parent}{$chr}{$num}++; # incrementation of the number of existing markers by windows.
					$AlleleRef = $allele_refalt{$parent}{$chr}{$pos};
					
					for my$formInd(sort keys %FormatHash){ 
						if($FormatHash{$formInd}{"DP"} ne "."){
							if ($FormatHash{$formInd}{"AD"} =~ m/^(\-\d+|\d+)\,([\d+\,\-]+)/){
								my $refer = $1; # number of reads for the base1
								my @alters = split(",",$2);
								my $sumreads = $FormatHash{$formInd}{"DP"}; # number of reads for the for a marqueur (base1 + base2)
								if($sumreads > 0){ # if base1 + base2 != 0
									if ($AlleleRef eq $base1){ # if base1 if the ancestral allele
										$frequences{$parent}{$formInd}{$chr}{$num}{$pos} = $refer / $sumreads; # frequency of the ancestral allele for this marker
										if (defined($nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Ancestral})){
											$nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Ancestral} = $nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Ancestral} + $refer;
											$nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Not} = $nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Not} + ($sumreads-$refer);
										}
										else{
											$nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Ancestral} = $refer;
											$nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Not} = ($sumreads-$refer);
										}
									}
									else{
										for(my$a=0; $a < scalar(@bases2); $a++){
											if ($AlleleRef eq $bases2[$a]){ # if one of the base2 if the ancestral allele
												$frequences{$parent}{$formInd}{$chr}{$num}{$pos} = $alters[$a] / $sumreads; # frequency of the ancestral allele for this marker
												if (defined($nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Ancestral})){
													$nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Ancestral} = $nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Ancestral} + $alters[$a];
													$nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Not} = $nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Not} + ($sumreads-$alters[$a]);
												}
												else{
													$nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Ancestral} = $alters[$a];
													$nbAncestralOrNot{$parent}{$formInd}{$chr}{$num}{$Not} = ($sumreads-$alters[$a]);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
close VCF;


my%percentage; # frequency of ancestral allele in a given window (ancestral reads / sum of ancestral and non ancestral reads). If there is no marker of a given widow in the vcf --> window is NA (will become NE later).
for my $parent(sort keys(%blocs)){
	for my $chr (sort keys(%{$blocs{$parent}})) {
		for my $num (sort {$a <=> $b} keys(%{$blocs{$parent}{$chr}})) {
            foreach my $accession (keys %hybrid){
                if (!exists($nbAncestralOrNot{$parent}{$accession}{$chr}{$num}{$Ancestral})) {
                	$percentage{$parent}{$accession}{$chr}{$num} = "NA";
                }
                else {
                	my $sum = $nbAncestralOrNot{$parent}{$accession}{$chr}{$num}{$Ancestral} + $nbAncestralOrNot{$parent}{$accession}{$chr}{$num}{$Not};
                	if ($sum == 0) {
                		$percentage{$parent}{$accession}{$chr}{$num} = "NA";
                	}
                	else {
                		$percentage{$parent}{$accession}{$chr}{$num} = $nbAncestralOrNot{$parent}{$accession}{$chr}{$num}{$Ancestral} / $sum;
                	}
                }
            }   
        }
    }
}
my @list_parent;
for my $parent (sort keys(%blocspos)) {
	push(@list_parent,$parent)
}
my %gemo;
my %coordinate;

my %data;
foreach my $accession (keys %hybrid){
    my $file_curve = $dirout."/". $accession."_curve.txt";
    my $file_ideogram = $dirout."/". $accession."_ideo.txt";
    my $file_color = $dirout."/". $accession."_color.txt";
    my $file_chromosome = $dirout."/". $accession."_chrom.txt";
    my $file_ancestor = $dirout."/". $accession."_ancestor.txt";
    my $ploidy = $hybrid{$accession};
    open (Ffreq, ">$file_ancestor") or die ("cant write $file_ancestor\n"); # output => frequency of ancestral alleles along chromosomes
    print Ffreq join("\t","Hybrid","Ancestry","Chromosome","Position_Start","Position_End","Frequence"),"\n"; # header of ancestorfreq
    for my $parent (sort keys(%blocspos)){
    	for my $chro (sort keys(%{$blocspos{$parent}})){
    		for my $num (sort {$a<=>$b} keys(%{$blocspos{$parent}{$chro}})){
    			my @position = split("-",$blocspos{$parent}{$chro}{$num});
    			my $start = $position[0];
    			# Correction of the begining of each bloc. -> +1 postition.
    			my $start_cor = 1;
    			if ($start != 1){
    				$start_cor = $start + 1;
    			}
    			my $end = $position[1];
    			print Ffreq join("\t",$accession , $parent,$chro ,$start_cor,$end ,$percentage{$parent}{$accession}{$chro}{$num}),"\n"; 
                push @{$data{$chro}} , {
                    start => $start_cor,
                    end   => $end,
                    hybrid => $parent,
                    percent => $percentage{$parent}{$accession}{$chro}{$num} 
                }; 
                push @{$coordinate{$chro}} , $start_cor;
            }
    	}
    }
    close Ffreq; 


    open(CURVE, ">$file_curve") or die("cant write $file_curve\n"); # output => frequency of ancestral alleles along chromosomes
    print CURVE join("\t","chr","start","end",join("\t",@list_parent)),"\n"; # header
    
    my %matrix;
    foreach my $chromosome (sort keys %coordinate) {
        my $last;
        foreach my $start (sort {$a <=> $b}  @{$coordinate{$chromosome}}) {
            if ($last) {
                next if $last == $start;
                my $end = $start - 1;
                my %percent;
                my @value;
                foreach my $data ( @{$data{$chromosome}}) {
                    
                    if ($last >= $data->{start} && $end <= $data->{end}) {   
                        push @{$matrix{$chromosome}{$last}{$end}}, {
                            hybrid => $data->{hybrid},
                            percent => $data->{percent}
                        }; 
                        $percent{$data->{hybrid}} =  $data->{percent};
                    }
                }
            
                for my $anc (@list_parent) {
                    if ($percent{$anc} ){
                        push @value,$percent{$anc};
                    }
                    else {
                        push @value,"0.0";
                    }
                }
            
                print CURVE join("\t",$chromosome, $last,$end,join("\t",@value)) ,"\n";
                $last = $start;
            }
            else {
                $last = $start;          
            } 
        }
    }
    close CURVE;   
    #--------------------------------------------------------------------------------------------------------------------
    # STEP 3 : Estimate allelic dosage of each ancestors by windows. We use the Maximum of likelihood compared to theoretical frequency (LOD). 
    
    
    
    my @dosage; # estimated frequencies.
    if($ploidy == 2){
        $dosage[0] = 1-$freq;
        $dosage[1] = 0.50;
        $dosage[2] = $freq;
    }
    if($ploidy == 3){
        $dosage[0] = 1-$freq;
        $dosage[1] = 0.33;
        $dosage[2] = 0.66;
        $dosage[3] = $freq;
    }
    if($ploidy == 4){
        $dosage[0] = 1-$freq;
        $dosage[1] = 0.25;
        $dosage[2] = 0.50;
        $dosage[3] = 0.75;
        $dosage[4] = $freq;
    }
    
    # calculation of LOD
    my%dosageAllelique; # key 1 = hybrid / key 2 = parent / key 3 = chromosome / key 4 = window number => value = allelic dosage of an ancestor for a window.
    my$dos; 
    my$ref; # Number of ancestral reads
    my$alt; # Number of reads non ancestral 
    for my $parent (sort keys %percentage) {
        for my $ind (sort keys %{$percentage{$parent}}) {
            for my $chr (sort keys %{$percentage{$parent}{$ind}}) {
                for my $num (sort {$a<=>$b} keys %{$percentage{$parent}{$ind}{$chr}}) {
                    if ($percentage{$parent}{$ind}{$chr}{$num} eq "NA"){
                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NE";
                    }
                    else{
                        $ref = $nbAncestralOrNot{$parent}{$ind}{$chr}{$num}{$Ancestral}; # ancestral reads
                        $alt = $nbAncestralOrNot{$parent}{$ind}{$chr}{$num}{$Not}; # non ancestral reads
        
                        if($ploidy == 2){ # simplified version
                            # hyp 0 vs hyp 1 :	
                            my$lod0vs1 = $ref*log10($dosage[0]/$dosage[1]) + $alt*log10($dosage[2]/$dosage[1]);
                            # If lod > 3 --> hyp 0 kept
                            if ($lod0vs1 > $LOD){ # hyp 0 > hyp 1
                                $dos = 0;
                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                            }
                            elsif($lod0vs1 < -$LOD){ # hyp 1 > hyp 0
                                # hyp 1 vs hyp 2 :
                                my$lod1vs2 = $ref*log10($dosage[1]/$dosage[2]) + $alt*log10($dosage[1]/$dosage[0]);
                                # If lod 3 --> hyp 1 kept
                                if ($lod1vs2 > $LOD){ # hyp 1 > hyp 2
                                    $dos = 1;
                                    $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                }
                                elsif($lod1vs2 < -$LOD){ # hyp 1 < hyp 2
                                    $dos = 2;
                                    $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                }
                                else{ # indetermination between hyp1 and hyp2
                                    $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
                                }
                            }
                            else{ # indetermination between hyp1 and hyp0
                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
                            }
                        }
                        elsif($ploidy == 4){ # simplified version
                            # hyp 0 vs hyp 1 :
                            my$lod0vs1 = $ref*log10($dosage[0]/$dosage[1]) + $alt*log10($dosage[4]/$dosage[3]);
                            # If lod > 3 --> hyp 0 kept
                            if ($lod0vs1 > $LOD){ # hyp 0 > hyp 1
                                $dos = 0;
                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                            }
                            elsif($lod0vs1 < -$LOD){ # hyp 1 > hyp 0
                                # hyp 1 vs hyp 2 :
                                my$lod1vs2 = $ref*log10($dosage[1]/$dosage[2]) + $alt*log10($dosage[3]/$dosage[2]);
                                # si lod > 3 --> hyp 1 kept
                                if ($lod1vs2 > $LOD){ # hyp 1 > hyp 2
                                    $dos = 1;
                                    $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                }
                                elsif($lod1vs2 < -$LOD){ # hyp 2 > hyp 1
                                    # hyp 2 vs hyp 3 :
                                    my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[2]/$dosage[1]);
                                    if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
                                        $dos = 2;
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                    }
                                    elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
                                        # hyp 3 vs hyp 4 :
                                        my$lod3vs4 = $ref*log10($dosage[3]/$dosage[4]) + $alt*log10($dosage[1]/$dosage[0]);
                                        if ($lod3vs4 > $LOD){ # hyp 3 > hyp 4
                                            $dos = 3;
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                        }
                                        elsif ($lod3vs4 < -$LOD){ # hyp 4 > hyp 3
                                            $dos = 4;
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                        }
                                        else{ #indetermination between hyp4 and hyp3
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
                                        }
                                    }
                                    else{ #indetermination between hyp3 and hyp2
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
                                    }
                                }
                                else{ #indetermination between hyp1 and hyp2
                                    $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
                                }
                            }
                            else{ #indetermination between hyp1 et hyp0
                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; # Indetermination 
                            }
        
        
                        }
                        elsif($ploidy == 3){ # normal version
                        # hyp 0 vs hyp 1 :	
                            my$lod0vs1 = $ref*log10($dosage[0]/$dosage[1]) + $alt*log10($dosage[3]/$dosage[2]);
                            # If lod > 3 --> hyp 0 kept
                            if ($lod0vs1 > $LOD){ # hyp 0 > hyp 1
                                $dos = 0;
                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                            }
                            elsif($lod0vs1 < -$LOD){ # hyp 1 > hyp 0
                                # hyp 1 vs hyp 2 :
                                my$lod1vs2 = $ref*log10($dosage[1]/$dosage[2]) + $alt*log10($dosage[2]/$dosage[1]);
                                # If lod > 3 --> hyp 1 kept
                                if ($lod1vs2 > $LOD){ # hyp 1 > hyp 2
                                    $dos = 1;
                                    $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                }
                                # else, hyp 2 vs hyp 3 :
                                elsif($lod1vs2 < -$LOD){ # hyp 2 < hyp 1
                                    my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
                                    # si lod > 3 --> hyp 2 kept
                                    if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
                                        $dos = 2;
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                    }
                                    # si lod < -3  --> hyp 3 kept
                                    elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
                                        $dos = 3;
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                    }
                                    else{ 
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                    }
                                }
                                else{ #indetermination between hyp1 and hyp2
                                    # hyp1 vs hyp3
                                    my$lod1vs3 = $ref*log10($dosage[1]/$dosage[3]) + $alt*log10($dosage[2]/$dosage[0]);
                                    if ($lod1vs3 > $LOD){ # hyp 1 > hyp 3;
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                    }
                                    else { # hyp 3 > hyp 1 ou indetermination
                                        # hyp 2 vs hyp 3
                                        my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
                                        # si lod > 3 --> hyp2 kept
                                        if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                        }
                                        # si lod < -3  --> hyp3 kept
                                        elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
                                            $dos = 3;
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                        }
                                        else{ 
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                        }
                                    }
                                }
                            }
                            else {  #indetermination between hyp0 and hyp1
                                my$lod0vs2 = $ref*log10($dosage[0]/$dosage[2]) + $alt*log10($dosage[3]/$dosage[1]);
                            # hyp 0 vs hyp 2 :
                                if ($lod0vs2 > $LOD){ # hyp 0 > hyp 2
                                    $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                }
                                # lod < 3 --> hyp 2 kept but we have to test others hypothesis.
                                elsif($lod0vs2 < -$LOD){ # hyp 2 > hyp 0
                                    # we test hyp 1 vs hyp2 :
                                    my$lod1vs2 = $ref*log10($dosage[1]/$dosage[2]) + $alt*log10($dosage[2]/$dosage[1]);
                                    # if lod > 3 --> hyp 1 kept
                                    if ($lod1vs2 > $LOD){ # hyp 1 > hyp 2
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                    }
                                    # else, hyp 2 vs hyp 3 :
                                    elsif($lod1vs2 < -$LOD){ # hyp 2 < hyp 1
                                        my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
                                    
                                        # if lod > 3 --> hyp 2 kept
                                        if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
                                            $dos = 2;
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                        }
                                        # si lod < -3  --> hyp 3 kept
                                        elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
                                            $dos = 3;
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                        }
                                        else{ 
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                        }
                                    }
                                    else{ #indetermination between hyp1 and hyp2
                                        # then hyp1 vs hyp3
                                        my$lod1vs3 = $ref*log10($dosage[1]/$dosage[3]) + $alt*log10($dosage[2]/$dosage[0]);
                                        if ($lod1vs3 > $LOD){ # hyp 1 > hyp 3
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                        }
                                        elsif ($lod1vs3 < -$LOD){
                                            # hyp 2 and hyp 3
                                            my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
                                            # if lod > 3 --> hyp 2 kept
                                            if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
                                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                            }
                                            # if lod < -3  --> hyp 3 kept
                                            elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
                                                $dos = 3;
                                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                            }
                                            else{ 
                                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                            }
                                        }
                                        else { # indetermination 1v3
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA";	
                                        }
                                    }
                                }
                                else{ # indetermination between 0 and 2
                                    # hyp 0 vs hyp 3
                                    my$lod0vs3 = $ref*log10($dosage[0]/$dosage[3]) + $alt*log10($dosage[3]/$dosage[0]);
                                    if ($lod0vs3 > 3){ # hyp 0 > hyp 3
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                    }
                                    elsif ($lod0vs3 < -$LOD){ # hyp 3 > hyp 0 --> we have to test other hypothesis
                                        # hyp 1 vs hyp 3 :
                                        my$lod1vs3 = $ref*log10($dosage[1]/$dosage[3]) + $alt*log10($dosage[2]/$dosage[0]);
                                        # if lod > 3 --> ID
                                        if ($lod1vs3 > $LOD){ # hyp 1 > hyp 3
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                        }
                                        # else, hyp 2 vs hyp 3 :
                                        elsif($lod1vs3 < -$LOD){ # hyp 2 < hyp 1
                                            my$lod2vs3 = $ref*log10($dosage[2]/$dosage[3]) + $alt*log10($dosage[1]/$dosage[0]);
                                            # if lod > 3 --> hyp 2 kept
                                            if ($lod2vs3 > $LOD){ # hyp 2 > hyp 3
                                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                            }
                                            # if lod < -3  --> hyp 3 kept
                                            elsif ($lod2vs3 < -$LOD){ # hyp 3 > hyp 2
                                                $dos = 3;
                                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = $dos;
                                            }
                                            else{ 
                                                $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                            }
                                        }
                                        else{ #indetermination between hyp1 and hyp3
                                            $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                        }
                                    }
                                    else{ # indétermintation between 0 and 3
                                        $dosageAllelique{$ind}{$parent}{$chr}{$num} = "NA"; #Indetermination 
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
     
    
    # --------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    # Lenght of chromosome determination.
    my%theEndChr; 
    my%theNumChr; 
    if ($booEndChro == 0){
        my$n = 0;
        for my$x(sort keys %endChro){
            if (grep { $_ eq $x } @chromosomList){
                $n++;
                $theEndChr{$x} = $endChro{$x};
                $theNumChr{$x} = $n;
            }
        }
        %endChro = %theEndChr;
    }
    elsif ($booEndChro == 1){
        my$n = 0;
        for my$x(sort keys %realEndChro){
            if (grep { $_ eq $x } @chromosomList){
                $n++;
                $theEndChr{$x} = $realEndChro{$x};
                $theNumChr{$x} = $n;
            }
        }
        %endChro = %theEndChr;
    }
    
    # step of keys inversion (parent/chrom/num) of %blocs to (chrom/parent/num) %blocsInv 
    my%blocsInv;
    for my$parent(sort keys %blocspos){
        for my $CHROMOSOMES(sort keys %{$blocspos{$parent}}){
            for my $num (sort {$a<=>$b} keys %{$blocspos{$parent}{$CHROMOSOMES}}){
                $blocsInv{$CHROMOSOMES}{$parent}{$num} = $blocspos{$parent}{$CHROMOSOMES}{$num};
            }
        }
    }
    
    # fill little windows with allelic dosages. The sum of every ancestors for one little window should be equal to the ploidy. Else it's an indetermination -> NA 
    
    my$cut=$window_cut*1000;
    my$startslice=1;
    my$endslice=$cut;
    my%dosbywindows;
    my$sum;
    #open(W,">windows");
    #print W "Sample	Chr	Window	Ancestor	Dose\n";
    for my$ind (sort keys %dosageAllelique){
        for my$CHROMOSOMES(sort keys %blocsInv){
            $startslice=1;
            $endslice=$cut;
            while ($startslice < $endChro{$CHROMOSOMES}){
                for my $parent(sort keys %{$blocsInv{$CHROMOSOMES}}){
                    $sum = 0;
                    for my$num (sort {$a<=>$b} keys %{$blocsInv{$CHROMOSOMES}{$parent}}){
                        my@markers = split("-",$blocsInv{$CHROMOSOMES}{$parent}{$num});
                        my$start = $markers[0];
                        my$end = $markers[1];
                        if($startslice >= $start && $startslice < $end){
                            $dosbywindows{$ind}{$CHROMOSOMES}{$startslice}{$parent} = $dosageAllelique{$ind}{$parent}{$CHROMOSOMES}{$num};
    #						print W "$ind   $CHROMOSOMES    $startslice     $parent $dosbywindows{$ind}{$CHROMOSOMES}{$startslice}{$parent}\n";
                        }
                    }
                }
                $startslice = $endslice;
                $endslice = $startslice + $cut;
            }
        }
    }
    #close(W);
    
    for my$ind (sort keys %dosbywindows){
        for my$CHROMOSOMES(sort keys %{$dosbywindows{$ind}}){
            for my$startslice(sort {$a<=>$b} keys %{$dosbywindows{$ind}{$CHROMOSOMES}}){
                foreach my$parent(@ref){
                    if (!exists($dosbywindows{$ind}{$CHROMOSOMES}{$startslice}{$parent})){
                        $dosbywindows{$ind}{$CHROMOSOMES}{$startslice}{$parent} = "NE"; 
                    } 
                }
            }
        }
    }
    
    my %hasher;
 #   my $ind = $accession;
    for my $CHROMOSOMES (sort keys %{$dosbywindows{$accession}}){
        for my $slice (sort {$a<=>$b} keys %{$dosbywindows{$accession}{$CHROMOSOMES}}){
            my$booNa = 0; # boolean if NA for at least one ancestor in the slice
            my%dosparent; # allelic dosage by parent for this slice
            my$dosBySlice = 0; # sum of allelic dosage for a slice
            for my $parent (sort keys %{$dosbywindows{$accession}{$CHROMOSOMES}{$slice}}){
                if ($dosbywindows{$accession}{$CHROMOSOMES}{$slice}{$parent} eq "NE"){ # NE = blocs with no corresponding markers in vcf => Dosage = 0
                    $dosbywindows{$accession}{$CHROMOSOMES}{$slice}{$parent} = 0;
    
                    if ($booNa == 0){
                        $dosparent{$parent} = $dosbywindows{$accession}{$CHROMOSOMES}{$slice}{$parent};
                        $dosBySlice = $dosBySlice + $dosparent{$parent};
                    }
                    elsif($booNa == 1){
                        $dosBySlice = "NA";
                    }
                }
                elsif($dosbywindows{$accession}{$CHROMOSOMES}{$slice}{$parent} eq "NA"){ # NA = bloc with LOD indetermintation -> NA
                    $dosparent{$parent} = "NA";
                    $booNa = 1;
                    $dosBySlice = "NA";
                }
                else{ # real dosage
                    if ($booNa == 0){
                        $dosparent{$parent} = $dosbywindows{$accession}{$CHROMOSOMES}{$slice}{$parent};
                        $dosBySlice = $dosBySlice + $dosparent{$parent};
                    }
                    elsif($booNa == 1){
                        $dosBySlice = "NA";
                    }
                }
            }
            # fill $hasher with color for each slice
            # If NA of != ploidy --> NA. Can't conclude anything.
            if ($dosBySlice eq "NA"){
                for (my$p = 0; $p < $ploidy; $p++){
                    my$slice2 = $slice + $cut;
                    $hasher{$CHROMOSOMES}{$p}{$slice} = "NA";
                }
            }
            elsif ($dosBySlice != $ploidy){
                for (my$p=0; $p < $ploidy; $p++){
                    my$slice2 = $slice + $cut;
                    $hasher{$CHROMOSOMES}{$p}{$slice} = "NA";
                }
            }
            # Else, if = ploidy. We can paint.
            elsif ($dosBySlice == $ploidy){
                my$slice2 = $slice + $cut;
                my$p = 0; # incrementation of p until the poidy each time we paint an haplotype
                for my$parent (sort keys %dosparent){
                    my$x = $dosparent{$parent}; # We paint as much haplotype than allelic dosage for this parent
                    if ($x != 0){
                        for(my$i = 1; $i <= $x; $i++){
                            if ($p == 0){
                                $hasher{$CHROMOSOMES}{$p}{$slice} = $parent;
                            }
                            elsif ($p == 1){
                                $hasher{$CHROMOSOMES}{$p}{$slice} = $parent;
                            }
                            elsif ($p == 2){
                                $hasher{$CHROMOSOMES}{$p}{$slice} = $parent;
                            }
                            elsif ($p == 3){
                                $hasher{$CHROMOSOMES}{$p}{$slice} = $parent;
                            }
                            $p++;
                        }
                    }
                }
    
            }
        }
    }
    open(F11, ">$file_ideogram") or die ("Error: can't open output for $file_ideogram\n"); # for a particular hybrid. Output = ideogram format.
    print F11 join("\t","chr","haplotype","start","end","ancestral_group"),"\n";
    for my $CHROMOSOMES (sort keys %hasher){
        my $endC = $endChro{$CHROMOSOMES};
        my $numchr = $theNumChr{$CHROMOSOMES}; 
        for my$p (sort {$a<=>$b} keys %{$hasher{$CHROMOSOMES}}){
            my$color = "";
            my$begin = 1;
            my$endS = 0;
            my$sli;
            for my$slice (sort {$a<=>$b} keys %{$hasher{$CHROMOSOMES}{$p}}){
                if ($color eq ""){
                    $color = $hasher{$CHROMOSOMES}{$p}{$slice};
                }
                # we join slice of the same color together before we print
                if ($color ne $hasher{$CHROMOSOMES}{$p}{$slice}){
                    $endS = $slice;
                    if ($endS <= $endC){
                        if($color ne "NA") {
                            print F11 join("\t",$numchr ,$p ,$begin ,$endS ,$color),"\n";
                        }
                    }
                    else{
                        if($color ne "NA"){
                            print F11 join("\t",$numchr, $p ,$begin ,$endC ,$color),"\n";
                        }
                    }
                    $begin = $slice + 1;
                }
                $color = $hasher{$CHROMOSOMES}{$p}{$slice};
                $sli = $slice;
            }
            if ($endS < $endC){
                if ($sli + $cut > $endS){
                    if ($sli + $cut >= $endC){ # If the chromosome longer than the end of the last slice -> grey 'til the end
                        $endS =	$endC;
                        if($color ne "NA") {
                            print F11 join("\t",$numchr, $p ,$begin,$endS, $color),"\n"; #color of the last slices
                        }
                    }
                    else{
                        if($color ne "NA"){
                            $endS =	$sli + $cut;
                            if($color ne "NA"){
                                print F11 join("\t",$numchr ,$p ,$begin ,$endS, $color),"\n"; #color of the last slices
                            }
                            $begin = $endS + 1;
                            $color = "NA";
                            if($color ne "NA") {
                                print F11 join("\t", $numchr , $p, $begin, $endC, $color),"\n"; #end of chro grey
                            }
                        }
                        else{ #if color = color of NA
                            if($color ne "NA") {
                                print F11 join("\t",$numchr, $p, $begin, $endC, $color),"\n";
                            }
                        }
                    }
                }
                else {
                    $begin = $endS + 1;
                    $color = "NA";
                    if($color ne "NA") {
                        print F11 join("\t",$numchr, $p, $begin, $endC, $color),"\n"; #end of chro grey
                    }
                }
            }
    
        }
    }
    close F11;
    # Print the length file.
    # chr length_chromosome haplotypes_names
    open(CHROMOSOME, ">$file_chromosome") or die ("Error : can't open $file_chromosome\n");
    print CHROMOSOME join("\t","chr","len","centromereInf","centromereSup","label"),"\n";
    for my $chromosome (sort keys %endChro){ 
    
        if ($ploidy == 2){
            print CHROMOSOME "$chromosome\t$endChro{$chromosome}\t\t\t01\n";
        }
        if ($ploidy == 3){
            print CHROMOSOME "$chromosome\t$endChro{$chromosome}\t\t\t012\n";
        }
        if ($ploidy == 4){
            print CHROMOSOME "$chromosome\t$endChro{$chromosome}\t\t\t0123\n";
        }
    }
    close CHROMOSOME; 

}

# fonction of log10
sub log10 {
     my $n = shift;
     return log($n)/log(10);
}