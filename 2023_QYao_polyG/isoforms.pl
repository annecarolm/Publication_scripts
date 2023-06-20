#!/usr/bin/perl
## ACMdS, Illinois Institute of Technology, 2021

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd qw(cwd);

my $name = 'isoforms.pl';
my $version = '0.7';
my $update = '2023-06-20';

my $usage=<<"USAGE";

NAME		${name}
VERSION		${version}
UPDATED		${update}

SYNOPSYS	This script was built for ROSMAP datasets from the AD knowledge portal and BaxD2.
		    It parses BAM alignment files to find insertion/deletions in reads mapped against Bax/
		    BaxD2 using the CIGAR string (6th column), then summarizes the findings into a tab- /
		    delimited table and counts polyG in reads per dataset. 
        It can also generate images of the aligned regions using BAMSNAP. 

REFERENCES	AD Knowledge Portal: https://adknowledgeportal.synapse.org/
		    SAM/BAM format: https://samtools.github.io/hts-specs/SAMv1.pdf

EXAMPLE		${name} \\
		  -bam *.bam \\
		  -b biospecimen.csv \\
		  -c clinical.csv \\
		  -o out.tsv

OPTIONS

# Required input
-bam                        BAM alignment file(s)
-c (--clinical)             ROSMAP clinical metadata in csv format
-b (--biospecimen)          ROSMAP biospecimen metadata in csv format

# Specified output
-o (--out)                  Table file name (must add the .tsv file extension)
-n (--num)                  Start numbering patients at specified value [Default = 1]

# Samtools view options
-ref (--reference_name)     Name of reference gene used for mapping [Default: BaxD2]
-r1 (--range1)              Start nucleotide position for SAM window [Default: 20]
-r2 (--range2)              End nucleotide position for SAM window [Default: 50]
-q (--MapQ)                 Filter bam file by uniqueness of alignment [Default: 10]
-F (--excl_flag)            Exclude alignments with specified samtools flag [Default: 256]   

# Bamsnap (https://github.com/parklab/bamsnap)
-bamsnap                    Generate alignment images with BAMSNAP
-fa (--fasta)               Required gene sequence in fasta format
-pos (--position_view)      Gene nucleotide position to view [Default: 30]
-m (--margin)               Genomic margin size [Default: 15]

USAGE
die "$usage\n" unless @ARGV;

## Defining variables

## Required input
my @bam;                    # List of input BAM files
my $clinical;               # ROSMAP clinical file
my $biospecimen;            # ROSMAP biospecimen file
my $table_out;              # output table name in tsv format
my $start = 1;              # Start numbering patients at
my $ref_name = 'BaxD2';     # Ref gene name

## Samtools specific options
my $range1 = 20;            # Start nucleotide position to filter SAM and BAM
my $range2 = 50;            # End nucleotide position to filter SAM and BAM
my $quality_score = 10;     # Quality score to filter reads in BAM file
my $flag = 256;             # Define samtools flag to exclude alignments from filtered files

## Only for bamsnap analysis
my $bamsnap;                # Flag to use bamsnap analysis
my $ref_fasta;              # Fasta of reference gene used for previous mapping
my $bamsnap_pos = 30;       # Gene nucleotide position to view
my $margin = 15;            # Genomic margin size

GetOptions(
    'bam=s@{1,}' => \@bam,
    'c|clinical=s' => \$clinical,
    'b|biospecimen=s' => \$biospecimen,
    
    'o|out=s' => \$table_out,
    'n|num=i' => \$start,
    
    'ref|reference_name=s' => \$ref_name,
    'r1|range1=i' => \$range1,
    'r2|range2=i' => \$range2,
    'q|MapQ=i' => \$quality_score,
    'F|excl_flag=i' => \$flag,
    
    'bamsnap' => \$bamsnap,
    'fa|fasta=s' => \$ref_fasta,
    'pos|position_view=i' => \$bamsnap_pos
);

## Generating SAM and BAM files for the specified reference range
## Setting up directory name for new SAM and BAM files
my $sam_dir = "${ref_name}.r${range1}-${range2}.SAM";
my $bam_dir = "${ref_name}.r${range1}-${range2}.BAM";

my @filt_bam;      # List of filtered BAM files
my @filt_sam;      # List of filtered SAM file name

my $path = cwd;    # Keeping absolute path to working directory
my $filtered_bam;  # Undef filtered BAM file
my $filtered_sam;  # Undef filtered SAM file
my $index_bam;     # Undef filtered BAM files to use with -bamsnap


#######################################################
### From BAM files
#######################################################

system "echo Generating filtered SAM and BAM files...";

while (my $bam_file = shift@bam){

    system "echo Working on $bam_file\n";

    ## Generating filtered SAM and BAM output file names
    my ($sam_name) = fileparse($bam_file);
    $sam_name =~ s/.bam/.filtered.sam/;
    my ($bam_name) = fileparse($bam_file);
    $bam_name =~ s/.bam/.filtered.bam/;
    
    ## Populating SAM and filtered BAM arrays for later use
    push(@filt_sam, $sam_name);
    push(@filt_bam, $bam_name);
    
    ## Creating directory to store SAM files
    unless (-d $sam_dir){
        mkdir ($sam_dir, 0755) or die "Cannot create $sam_dir: $!\n";
    }

    ## Generating the SAM files mapped to a specific range of the ref gene
    $filtered_sam = "${path}/${sam_dir}/${sam_name}";
    
    # BAM files must be indexed
    unless (-f "$bam_file.csi"){
        print "Indexing $bam_file\n";
        system (
            "samtools \\
                index \\
                -c \\
                $bam_file"
        ); 
    }

    system (
        "samtools view \\
            -h \\
            -q $quality_score \\
            -F $flag \\
            $bam_file \\
            ${ref_name}:${range1}-${range2} \\
            > $filtered_sam"
    );
    
    ## Creating directory to store filtered BAM files
    unless (-d $bam_dir){
        mkdir ($bam_dir, 0755) or die "Cannot create $bam_dir: $!\n";
    }
    
    ## Generating the BAM files mapped to a specific range of the ref gene
    $filtered_bam = "${path}/${bam_dir}/${bam_name}";
        
    system (
        "samtools view \\
            -h \\
            -q $quality_score \\
            -F $flag \\
            $bam_file \\
            ${ref_name}:${range1}-${range2} \\
            -bo $filtered_bam"
    );

    if($bamsnap){
        $index_bam = "${path}/${bam_dir}/${bam_name}.bai";
        
        system (
            "samtools index \\
                $filtered_bam \\
                $index_bam"
        );
    }

}

#######################################################
### Biospecimen file
#######################################################

## Creating biospecimen db # grabbing the individual ID from the biospecimen file based on sample id
system "echo Working on biospecimen file...";
system "dos2unix $biospecimen";     # Converting file from windows to unix

my %biospecimen;

open BIOS, "<", "$biospecimen" or die "$!\n";

my $assay_num;
while (my $line = <BIOS>){

    chomp $line;

    ## Getting rid of the header

    if ($line =~ /^individualID/){ 
        my @header = split(',', $line);
        for my $num (0..$#header){
            if ($header[$num] =~ /assay/){
                $assay_num = $num;
            }
        }

    }
    
    ## Retrieving individual ID based on biospecimen ID 
    else {
        my @columns = split (',', $line);
        my $individual_id = $columns[0];
        my $biospecimen_id = $columns[1];
        my $assay = $columns[$assay_num];

        if ($assay =~ /rnaSeq/i){ $biospecimen{$biospecimen_id} = $individual_id; }
    }
}

#######################################################
### Clinical file
#######################################################

## Creating database of clinical data based on individual ID
system "echo Working on clinical file...";
system "dos2unix $clinical"; # Converting file from windows to unix

my %clinic;

open CLI, "<", "$clinical" or die "$!\n";

## Working on clinical csv file
while (my $line = <CLI>){
    chomp $line;
   
    ## Getting rid of quotes ""
    $line =~ s/\"//g;

    ## Grabbing clinical information
    my @columns = split (',', $line);
    my $id = $columns[17];
    $clinic{$id}{'msex'} = $columns[2];
    $clinic{$id}{'race'} = $columns[4];
    $clinic{$id}{'spanish'} = $columns[5];
    $clinic{$id}{'apoe_genotype'} = $columns[6];
    $clinic{$id}{'braaksc'} = $columns[13];
    $clinic{$id}{'ceradsc'} = $columns[14];
    $clinic{$id}{'cogdx'} = $columns[15];
    $clinic{$id}{'dcfdx_lv'} = $columns[16];

}

#######################################################
#######################################################

## Creating the output table with parsed info
open TSV, ">", "$table_out" or die "Can't create $table_out: $!\n";

## Printing table header
my @out_columns = (
    'msex',
    'race',
    'spanish',
    'apoe_genotype',
    'braaksc',
    'ceradsc',
    'cogdx',
    'dcfdx_lv',
    'BaxD2',
    'iso1',
    'iso2',
    'iso3',
    'G_count'
);

print TSV "\#Patient";
foreach (@out_columns){ print TSV "\t$_"; }
print TSV "\n";

## Deidentifying patient data name, start patient numbering from defines in command line 
my $patient_count = $start - 1;

## Creating directory for deidentified SAM files
my $sam_deid_dir = 'deidentified_sam';

unless (-d $sam_deid_dir){ 
    mkdir ($sam_deid_dir, 0755) or die "Cannot create $sam_deid_dir: $!\n";
}


###########################################################
### Working on filtered SAM, deidentifying data and BAMSNAP
###########################################################

my %isoforms;
my %Gs;

## Working on filtered SAM files
while (my $sam = shift@filt_sam){
    
    $patient_count++;
    $patient_count = sprintf("%05d", $patient_count);
    my $sam_deid = "patient_${patient_count}_deidentified.sam";

    ## Setting paths to open and write files to
    my $path_sam = "${path}/${sam_dir}/$sam";
    my $path_deid = "${path}/${sam_deid_dir}/$sam_deid";
    
    ## Getting the right bam file
    my $bam_path = "${path}/${bam_dir}";
    my $bam = $sam;
    $bam =~ s/sam$/bam/;
    $filtered_bam = "$bam_path/$bam";

    open SAM, "<", "$path_sam" or die "$!\n";     # opening filtered sam
    open DEID, ">", "$path_deid" or die "$!\n";   # writing deidentified sam 

    my $sID;                # specimen ID
    my $read_count = 0;
    my $sample_id;
    my $SM;

    my $read_id;
    my $start_align;
    my $end_align;
    my $cigar;
    my $sam_flag;
    my $seq;

    $isoforms{'BaxD2'} = 0;
    $isoforms{'iso1'} = 0;
    $isoforms{'iso2'} = 0;
    $isoforms{'iso3'} = 0;
    $Gs{'G_count'} = undef;

    while (my $line = <SAM>){

        my %db;
                
        chomp $line;
        
        ## Keeping SAM file header
        if ($line =~ /^\@HD/){ print DEID "$line\n"; }
        
        ## Keeping reference alignment sequence info
        elsif ($line =~ /^\@SQ/){ print DEID "$line\n"; }
        
        ## Getting the biospecimen id from the sam file and deidentifiying data name
        elsif ($line =~ /^\@RG\s+ID:(\S+)\s+(\S+)/){
            my $sample_id = $1;                                 # dataset name
            my $SM = $2;                                        # SM is reference mapped to
            my $base_name = fileparse($sample_id);              # grabbing only file name
            ($sID) = $base_name =~ /^(\w+)\./;                  # parsing data name to get biospecimen ID
            $sID =~ s/_S\d+_R[12]_\d+$//;                       # parsing file name to match biospecimen ID formatting
            print "$sID\n";
            my $id_len = "patient_${patient_count}.fastq.gz";   # deidentifying dataset name
            print DEID "\@RG ID:$id_len\t$SM\n";                # printing out to deidentified filtered SAM file
        }    
        
        ## Skipping command line details
        elsif ($line =~ /^\@PG/){ next; }
        
        ## Working on alignment data
        else {
            $read_count++;
            $read_count = sprintf("%02d", $read_count);

            ## Deidentifying lines and grabbin cigar info
            my @columns = split("\t", $line);
            ($read_id) = $columns[0] = "deidentified_read_$read_count";                 # Removing read ID
            $sam_flag = $columns[1];                                                    # Get flag to determine read mapping orientation
            $start_align = $columns[3];                                                 # start position of alignment
            $end_align = $columns[4];                                                   # end position of alignment
            $cigar = $columns[5];                                                       # Grabbing CIGAR info
            $seq = $columns[9];                                                         # Getting DNA sequence
            my $seq_length = 'N' x length($seq);                                        # Getting length of DNA seq and replacing by N
            $columns[9] = $seq_length;                                                  # Replacing DNA seq on file
            $columns[11] = "RG:Z:patient_${patient_count}.fastq.gz";                    # Deidentifying dataset name

            ## Printing deidentified SAM file
            for (0..$#columns - 1){ print DEID "$columns[$_]\t"; }
            print DEID "$columns[$#columns]\n";
            
            ## Trying to split cigar info
            my @split_cigar = split(/(\d+\S)/, $cigar);
                       
            ## Removing empty elements of cigar split array
            for my $item (@split_cigar){
                if ((defined $item) and !($item =~ /^$/)) {
                    push ( @{$db{$read_id}}, $item ); ## Hash of arrays
                }
            }
            ## Counting the number of Gs
            if ($seq =~ /(G{7,})/ ){
                my $G = $1;
                my $length = sprintf("%02d", length($G));
                if ($start_align <= 31){
                    $Gs{'G_count'}{$read_count} = $length;
                }
            }
        }

        my $match_pos;
        my $del_pos;
        
        ## Working on isoform identification from the read number
        for my $read (sort (keys %db)) {
            
            foreach (@{$db{ $read_id }} ) {
                my $cigar_align = $_;

                ## Skipping the soft clip reads
                if ($cigar_align =~ /(\d+S|\d+H)/){ next; }

                ## Grabbing the Match read info
                elsif ($cigar_align =~ /(\d+)M/){
                    my $match = $1;
                    $match_pos = ($start_align + $match) -1;    # correcting base position (-1)

                    ## Counting how many BaxD2 mapped, and keeping only the reads that matched from position 31 and later
                    if ($start_align <= 31 && $match_pos >= 50){ 
                        $isoforms{'BaxD2'} += 1; 
                    }
                    else { next; }
                }

                ## Calculating position of DEL
                elsif ($cigar_align =~ /(\d+)D/){
                    my $del = $1;
                    $del_pos = ($match_pos + $del) -1;      # correcting base position (-1)
                    
                    ## Separating double deletion isoforms 
                    if ($del == 2){

                        ## Counting how many reads are isoform 1
                        if ($del_pos == 36){ $isoforms{'iso1'} += 1; }

                        ## Counting how many reads are isoform 2
                        elsif ($del_pos == 31) { $isoforms{'iso2'} += 1; }
                    }
                    ## Grabing isoform 3
                    elsif ($del == 1){
                        
                        ## Counting how many reads are isoform 3 based on first deletion position 
                        if ($del_pos == 30){ $isoforms{'iso3'} += 1; }
                    }
                    last;
                }
            }  
        }
    }

    ## Getting clinical data for each patient based on biospecimen and individual IDs
    my $in_ID = $biospecimen{$sID};
    
    if (exists $clinic{$in_ID}){
        print TSV "$patient_count";
        foreach (@out_columns){ 

            if (exists $clinic{$in_ID}{$_}){
                print TSV "\t$clinic{$in_ID}{$_}";
            }
            
            if (exists $isoforms{$_}){
                print TSV "\t$isoforms{$_}";
            }
        } 
        print TSV "\t";
        
        ## Add the Gs 'G_count'
        if (defined $Gs{'G_count'}){
            foreach my $read_num (sort( keys %{ $Gs{'G_count'} })){
                my $g = $Gs{'G_count'}{$read_num};
                print TSV "read_${read_num}: ${g}Gs;";
            }
        }
        else {
            print TSV "NA";
        }
        print TSV "\n";
    }
    else {
        print "$in_ID not found in \$clinic\n";
    }

    ## Bamsnap analysis: deidentifyied outputs
    
    if ($bamsnap){

        # Creating bamsnap directory
        my $bamsnap_dir = 'bamsnap';

        unless (-d $bamsnap_dir){
            mkdir ($bamsnap_dir, 0755) or die "Cannot create $bamsnap_dir: $!\n";
        }
        
        # Generating bamsnap files
        my $bamsnap_out = "${path}/${bamsnap_dir}/patient_${patient_count}_${ref_name}_pos${bamsnap_pos}_m${margin}.png";

        system (
            "bamsnap \\
            -bam $filtered_bam \\
            -ref $ref_fasta \\
            -pos ${ref_name}:${bamsnap_pos} \\
            -out $bamsnap_out \\
            -title patient_${patient_count}_${ref_name}_pos${bamsnap_pos}_m${margin} \\
            -draw coordinates bamplot \\
            -read_group strand \\
            -margin $margin"
        );
    }
}
