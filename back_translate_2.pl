#!/usr/bin/perl -w
use strict;


my $usage = "\nusage: perl back_translate_2.pl <coding sequence input file (fasta format)> <protein alignment input file (aligned fasta format)> <codon alignment output file (aligned fasta format)> <log output file>\n\n".
            "This script produces codon alignments using protein alignments and the DNA sequences that code for them.The back-translation is checked for error.  Inconsistencies will produce 'Undefined value' warnings.  They are handled by printing 'xxx' in place of the codon, and noting the exception in a log file.  The script assumes the input files to be FASTA formatted, with headers and sequences in single, alternating lines. Headers are searched for regular expressions found in lines 134 and 186 of the script.\n\nWritten by Christopher J. Grassa at the University of Florida, under the instruction of Dr. Gordon Burleigh.\nDecember 3rd, 2009\n\n";
            


my $cdsfile = shift or die $usage;
my $alignment_in_file = shift or die $usage;
my $alignment_out_file = shift or die $usage;
my $log_file = shift or die $usage;

my(%geneticcode) =
    (
    'TCA' => 'S', # serine
    'TCC' => 'S', # serine
    'TCG' => 'S', # serine
    'TCT' => 'S', # serine
    'TCN' => 'S', # serine
    
    
    'TTC' => 'F', # phenylalanine
    'TTT' => 'F', # phenylalanine
    'TTA' => 'L', # leucine
    'TTG' => 'L', # leucine
    'TAC' => 'Y', # tyrosine
    'TAT' => 'Y', # tyrosine
    'TAA' => '*', # stop
    'TAG' => '*', # stop
    'TGC' => 'C', # cysteine
    'TGT' => 'C', # cysteine
    'TGA' => '*', # stop
    'TGG' => 'W', # tryptophan
    
    'CTA' => 'L', # leucine
    'CTC' => 'L', # leucine
    'CTG' => 'L', # leucine
    'CTT' => 'L', # leucine
    'CTN' => 'L', # leucine
    
    
    'CCA' => 'P', # proline
    'CCC' => 'P', # proline
    'CCG' => 'P', # proline
    'CCT' => 'P', # proline
    'CCN' => 'P', # proline
    
    'CAC' => 'H', # histidine
    'CAT' => 'H', # histidine
    'CAA' => 'Q', # glutamine
    'CAG' => 'Q', # glutamine
    
    'CGA' => 'R', # arginine
    'CGC' => 'R', # arginine
    'CGG' => 'R', # arginine
    'CGT' => 'R', # arginine
    'CGN' => 'R', # arginine
    
    
    'ATA' => 'I', # isoleucine
    'ATC' => 'I', # isoleucine
    'ATT' => 'I', # isoleucine
    'ATG' => 'M', # methionine
    
    'ACA' => 'T', # threonine
    'ACC' => 'T', # threonine
    'ACG' => 'T', # threonine
    'ACT' => 'T', # threonine
    'ACN' => 'T', # threonine
    
    
    'AAC' => 'N', # asparagine
    'AAT' => 'N', # asparagine
    'AAA' => 'K', # lysine
    'AAG' => 'K', # lysine
    'AGC' => 'S', # serine
    'AGT' => 'S', # serine
    'AGA' => 'R', # arginine
    'AGG' => 'R', # arginine
    
    'GTA' => 'V', # valine
    'GTC' => 'V', # valine
    'GTG' => 'V', # valine
    'GTT' => 'V', # valine
    'GTN' => 'V', # valine
    
    
    'GCA' => 'A', # alanine
    'GCC' => 'A', # alanine
    'GCG' => 'A', # alanine
    'GCT' => 'A', # alanine
    'GCN' => 'A', # alanine
    
    'GAC' => 'D', # aspartic acid
    'GAT' => 'D', # aspartic acid
    'GAA' => 'E', # glutamic acid
    'GAG' => 'E', # glutamic acid
    
    'GGA' => 'G', # glycine
    'GGC' => 'G', # glycine
    'GGG' => 'G', # glycine
    'GGT' => 'G', # glycine
    'GGN' => 'G', # glycine
    );



###################################################
##                                               ## 
##    SECTION 1: Storing the Coding Sequences    ##
##                                               ##
###################################################

# open the file containing the coding sequences for all species
# build a giant hash with 'Transcript ID' => 'coding sequence'
# the file needs to be formatted as alternating lines of header, sequence


open CDSFILE, "<", $cdsfile or die "could not locate coding sequences\n";;

my $cds_transcript_id = "";
my $coding_sequence = "";
my %coding_transcripts = ();

my $hashcounter = 0;
    
while ($cdsfile = <CDSFILE>)
{
    my $fileline = $cdsfile;
    chomp $fileline;
    # match a header line; capture the Transcript ID
    if ($fileline =~ m/^>ENS[A-Z]{3}G\d{11}\|(ENS[A-Z]{3}T\d{11})/)
    {
        ++$hashcounter;
        print "$hashcounter nucleotide sequences loaded so far\n"; 
        $cds_transcript_id = "$1";    
    }
    
    #if not a header line, must be a sequence line
    else
    {
        $coding_sequence = $fileline;
        
        #load the hash
        $coding_transcripts{"$cds_transcript_id"} = "$coding_sequence";
        
        #unnecessarily clear variables, for peace of mind
        $cds_transcript_id = "";  
        $coding_sequence = "";
    }
}



#################################################################
##                                                             ## 
##    SECTION 2: Back-translating the amino acid alignmnets    ##
##                                                             ##
#################################################################


open INALIGN, "<", $alignment_in_file or die "could not locate protein alignmnets\n";


open OUTALIGN, ">", $alignment_out_file or die "could not create codon alignmnet file\n";


open LOG, ">", $log_file or die "could not create log alignmnet file\n";

my $transcript_id;
my $header;
my $set_number;
my $amino_acid_sequence;
my @amino_acid_sequence;
my @codons;
my $counter = 0;
my $translation;
while ($alignment_in_file = <INALIGN>)
{
    my $fileline = $alignment_in_file;
    chomp $fileline;
    
    # match a header line; capture the Orthologous Set Number & Transcript ID
    if ($fileline =~ m/^>(\d+)\|ENS[A-Z]{3}G\d{11}\|(ENS[A-Z]{3}T\d{11})/)
    {
        ++$counter;
        $set_number = "$1";
        $transcript_id = "$2";
        $header = $fileline;
        print OUTALIGN "\n$header\n";
        print "back-translating sequence number $counter, found in set number $set_number\n";
    }
    
    #if not a header line, must be a sequence line
    else
    {
        $amino_acid_sequence = $fileline;
        #split amino acids on nothing
        @amino_acid_sequence = ($amino_acid_sequence =~ m/./g);
        
        #fetch the cds
        $coding_sequence = $coding_transcripts{"$transcript_id"};
        #split into codons
        @codons = ($coding_sequence =~ m/.../g);
        
        foreach (@amino_acid_sequence)
        {
            my $aminoacid = $_;
            
            # if an amino acid indel, print three nucleotide indels
            if ($aminoacid eq "-")
            {
                print OUTALIGN "---";
            }
            
            # otherwise, print its codon
            else
            {
                my $codon = shift @codons;
                #(upper case)
                $codon = uc $codon;
                #print "$codon";
                #fetch codon translation
                
                $translation = $geneticcode{"$codon"};
                #print "$translation\n";
                
                #verify that it codes for the amino acid
                if ($translation eq $aminoacid)
                {
                    print OUTALIGN "$codon";
                    
                }
                #if not, print xxx, and log it
                else
                {
                    print OUTALIGN "xxx";
                    print LOG "codon does not translate to amino acid in sequence $header of set $set_number\n";
                }

            
            }
        }
        
        #clear variables, just to be sure
        $transcript_id="";
        $header="";
        $set_number="";
        $amino_acid_sequence="";
        @amino_acid_sequence="";
        @codons="";
    }



}