#!/usr/bin/perl
 

use warnings;
use strict;
use lib 'lib';
use Parallel::ForkManager;
use Getopt::Long; 
use FindBin;
use Pod::Usage;
 
my $origin;
my @vcf;
my $individuals;
my $method = "vcfhunter"; 
my $threads = 1; 
my $color = "";
my $help = "";  
my $clean;
my $dirout = "";
my $usage = $FindBin::Bin ."/". $FindBin::Script.q/ --help
Parameters :
    -v, --vcf         A vcf file [required]
    -o, --origin       A two column file with individuals in the first column and group tag (i.e. origin) in the second column [Required]
    -i, --individuals  List of individuals to scan from vcf, as defined in the VCF headerline [Required]
    -m, --method      Permissible values: vcfhunter traceancestor (String). Default vcfhunter
    -c, --color       A color file with 4 columns: col1=group and the three last column corresponded to RGB code.
    -t, --threads     Number of threads
    -d, --dirout      Path to the output directory (Default method option name)
    -h, --help        display this help
/;
Getopt::Long::Configure ('bundling');
GetOptions (
	'o|origin=s'     => \$origin, ## ancestor file
    'v|vcf=s'       => \@vcf, ## vcf of hybrids. Hybrids should not have space in their names and no special characters. 
	'i|individuals=s' => \$individuals, 
	'c|color=s'     => \$color, 
	'm|method=s'    => \$method,
	't|threads=i'    => \$threads,  
    'd|dirout=s'      => \$dirout,
	'h|help!'       => \$help
)
or pod2usage(-message => $usage);



if ($help) { pod2usage(-message => $usage); }
if (!$origin){
    warn $usage;
	warn "\nWarn :: --origin is empty. Please specify a group file\n\n";
    exit 0; 
}  
 
if (!@vcf){
    warn $usage;
	warn "\nWarn :: --vcf is empty. Please specify a vcf file\n\n";
    exit 0; 
}
if (!$individuals){
    warn $usage;
	warn "\nWarn :: --individuals is empty. Please specify an accesion file\n\n";
    exit 0; 
}
 
my $pm = Parallel::ForkManager->new($threads);
$dirout = $method if $dirout eq "";
my $chromosome = &vcf_header($vcf[0],"chromosome$$.txt");

if ($method eq "vcfhunter") {
    &exec("mkdir -p temp/step1 temp/step2 temp/step3 temp/step4 $dirout");
    open(CONF,">list_vcf.txt");
    print CONF join("\n",@vcf) ,"\n";
    close CONF;
    my $cmd_step1 = "bin/IdentPrivateAllele.py -c list_vcf.txt -g ".$origin ." -o temp/step1 -a y  -m y -t $threads";
    &exec($cmd_step1);
    my $cmd_step2 = "bin/allele_ratio_group.py -g ". $origin. " -p _ratio.tab.gz -o temp/step2 -i temp/step1";
    &exec($cmd_step2);
    open(ACCESSION,$individuals);
    LOOP:
    while(<ACCESSION>){
        chomp;
        my $pid = $pm->start and next LOOP;
        my ($accession,$ploidy) = (split(/\t/,$_));
        my $cmd_step3 = "bin/allele_ratio_per_acc.py -c list_vcf.txt -g ".$origin. " -i temp/step2 -o temp/step3 -a ". $accession;
        &exec($cmd_step3);
        
        my $cmd_step4 = "bin/PaintArp.py -a ". $accession ." -r temp/step3/". $accession."_ratio.tab.gz -c ".$color ."  -o temp/step4/".$accession ."  -w 12 -O 0 -s ".$chromosome ." -p $ploidy";
        &exec($cmd_step4);
        
        my $cmd_step5 = "bin/convertForIdeo.py --name ".$accession ." --dir temp/step4/ --col ".$color ." --size ".$chromosome ." --prefix ".$dirout. "/".$accession ." --plo ".$ploidy;
        &exec($cmd_step5);
        $pm->finish;
    }
    $pm->wait_all_children; 
    &exec("rm -Rf temp list_vcf.txt $chromosome");


}
else {
    &exec("mkdir -p $dirout");
    my $cmd_step1 = "bin/vcf2gst.pl --vcf ". $vcf[0] ." --ancestor ". $origin ." --output ". $dirout. "/GSTMatrice.txt";
    &exec($cmd_step1);
    my $cmd_step2 = "bin/prefilter.pl --input ". $dirout. "/GSTMatrice.txt --output ". $dirout. "/Diagnosis_matrix.txt";
    &exec($cmd_step2);
    my $cmd_step3 = "bin/TraceAncestor.pl --vcf ". $vcf[0] ." --dirout ". $dirout ."  --input ". $dirout. "/Diagnosis_matrix.txt --hybrid ".$individuals;
    &exec($cmd_step3);
}



sub exec {
    my $cmd = shift; 
    system($cmd);
}

sub vcf_header{
    my ($in , $out) = @_;
    open(OUT,">$out");

    my @header = `grep "^##" $in`;

    foreach my $line(@header){
        chomp($line);
        if ($line =~ /^#.*ID=(.*),length=(\d+).*/){
            print OUT join("\t",$1,$2),"\n";
        }
    }
    close OUT;
    return $out;
}
