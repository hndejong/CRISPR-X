#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

#prend en entrée un fichier vcf des positions à tester et la liste des bam sur lesquels effectuer l'étude
#en sortie tous les reads qui overlappent le candidat, la position au sein du read, l'allele vu sur le read, l'allele porté par l'adn, l'allele rdd 

my @a_bam_list;
my @snp;

########################################OPTIONS############################################
# Options definitions (name,format,default,required)
my @Options = ( 
    ['h|help',undef,undef,0],
    ['m|man',undef,undef,0],
    ['v|vcf','s',undef,1],
    );

# defaults values
my %Values=();
my %Required=();
map {$Values{$_->[0]}=$_->[2];$Required{$_->[0]}=$_->[3];} @Options;

# build options list
my %Options=();
map {$Options{defined($_->[1])?$_->[0].'='.$_->[1]:$_->[0]}=\$Values{$_->[0]};} @Options;

# retrieve options
GetOptions(%Options) || pod2usage(2);
pod2usage(1) if ($Values{'help'});
pod2usage(-exitstatus => 0, -verbose => 2) if ($Values{'m|man'});

# check required
map {pod2usage("$_ is required") if ($Required{$_}&&!defined($Values{$_}));} keys %Values;

# BAM provided as arguments
if (scalar(@ARGV)!=0)  
{
    @a_bam_list = @ARGV;
}
else {
    print STDERR "\n\nCANNOT PROCEED\nat least 1 BAM file must be provided\n\n";
}
##########################################################################################

open(VCF,$Values{'v|vcf'})||die "Can't open file $Values{'v|vcf'}";
my %coord=();


while (my $line=<VCF>){
    chomp $line;
    next if ($line =~ /#/);
    my @a_vcf      = split ("\t",$line);
    
    my $region       = $a_vcf[0];
    my $snp_position = $a_vcf[1];
    my $DNA_allele = $a_vcf[3];
    my $alt_allele=$a_vcf[4];
    my $coord        = $region.":".$snp_position."-".$snp_position;    
    
    # Process each vcf position for all provided BAMs
    foreach my $bam (@a_bam_list){
	my $reads;
	my @reads;
	$reads     = `samtools view -F 4 $bam $coord`;
	@reads     = split ("\n",$reads);
	
	foreach my $read (@reads){	
	    next if (!defined $read);
	    $bam =~m/(\d+_\w+)_L\S+/;
	    my $bam_name   = $1;
	    my @a_read     = split ("\t",$read);
	    
	    # get read sequence & build an array with bases
	    my $sequence   = $a_read[9];
	    my @bases	   = split ("", $sequence);

	    my $snp_name   = $coord;
	    my $readName   = $a_read[0];
	    my $read_coord = $a_read[3];
	    my $CIGAR      = $a_read[5];
	    my $RG         = $a_read[12];
	    my $direction;	    
	    my $flag       = $a_read[1];
	    if ($flag & 16){
		$direction = "reverse";
	    }
	    else{
		$direction = "forward";
	    }    		      
	    
	    my $snp_distance_from_read_start = $snp_position - $read_coord + 1;			    
	    my $distance_from_snp  = $snp_distance_from_read_start;
	    my $position_on_read   = 0;
	    my $letter;
	    my $count;
	    
	    # CIGAR LINE Management
	    $CIGAR =~ s/(\d+\w)/$1\t/g;
	    my @F_cigar = split( /\t/, $CIGAR );
	    
	    foreach my $couple (@F_cigar) {
		$couple    =~ /(\d+)(\w)/;		
		my $digit     = $1;
		my $letter    = $2;
		#print join ("\t",$digit,$letter,$distance_from_snp,"\n");
		if ( $letter eq 'D' || $letter eq 'N'){		    
		    #print "D|N\n";
		    if ($distance_from_snp <= $digit){
			$position_on_read = 0;
			last;
		    }else{
			$distance_from_snp -= $digit;
		    }
		}
		if ( $letter eq 'S' ){		   
		    #print "S\n";
		    if ($distance_from_snp <= $digit){
			$position_on_read = 0;
			last;
		    }else{
			$distance_from_snp -= $digit;
			$position_on_read  += $digit;
		    }
		}
		if ( $letter eq 'I' || $letter eq 'M'){
		    #print "I|M\n";
		    if ($distance_from_snp <= $digit){
			#print $distance_from_snp."\t".$digit."\n";
			$position_on_read  += $distance_from_snp; 
			last;
		    }else{
			$distance_from_snp -= $digit;
			$position_on_read  += $digit;
		    }
		}
	    }
	    my $allele = $bases[$position_on_read - 1];
	    $CIGAR =~s/\t//g;
	    print join ("\t",$snp_name,$readName,$read_coord,$direction, $CIGAR,$position_on_read,$allele,$DNA_allele,$alt_allele,"\n") if ($position_on_read ne '0');    
	}
    }
}
close VCF;
#=======================================================================

=pod
    
    =head1 NAME

    snp_position_on_reads.pl
    
=head1 SYNOPSIS

snp_positions_on_reads.pl -v <VCF FILE> <BAM1> ...

=head2 OPTIONS :

    COMMAND              TYPE     DETAILS

    -v | --vcf            <FILE>  VCF file
    -h | --help           -       View this documentation
    --man                 -       View a more exaustive documentation 

=head1 DESCRIPTION :

gives snp coordinates, read name, read start, direction, mapping CIGAR line & snp position on read for each position given by the vcf input file 

=head1 AUTHORS
 
 DJARI Anis
 KLOPP Christophe
 Questions can be posted to the sigenae mailing list:
 sigenasupport@jouy.inra.fr

=head1 VERSION

1

=head1 DATE

06/2013

=cut
