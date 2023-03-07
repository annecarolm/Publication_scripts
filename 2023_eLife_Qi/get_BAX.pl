#!/usr/bin/perl
## ACMdS, Pombert Lab, Illinois Institute of Technology, 2023

use warnings;
use strict;
use Getopt::Long qw(GetOptions);
use File::Basename;
use Cwd qw(cwd);

my $script_name = 'get_BAX.pl';
my $version = '0.1a';
my $update = '2023-02-22';

my $usage=<<"USAGE";

NAME		  ${script_name}
VERSION		${version}
UPDATED		${update}

SYNOPSYS	This script was built for RNA-seq data from the ROSMAP study from the AD knowledge portal.
		    It parses BAM alignment files to determine the % of reads that maps to the Bax/BaxD2 juntions.
            It generates deidentified .SAM files of mapped region and summarizes the findings into a tab-delimited table. 
            It can also generate images of the aligned regions using BAMSNAP. 

REFERENCES	AD Knowledge Portal: https://adknowledgeportal.synapse.org/
		    SAM/BAM format: https://samtools.github.io/hts-specs/SAMv1.pdf

EXAMPLE		${script_name} \\
		  -bam *.bam \\
		  -b biospecimen.csv \\
		  -c clinical.csv \\
		  -o out.tsv

OPTIONS

# Required input
-bam                        BAM alignment file(s)
-c (--clinical)             ROSMAP clinical metadata in csv format
-b (--biospecimen)          ROSMAP biospecimen metadata in csv format
-n (--num)                  Start numbering patients at specified value [Default = 1]

# Output options
-o (--out)                  Table file name (must add the .tsv file extension) [Default: out.tsv]
-osam (--outdir_sam)        Directory name for deidentified SAM files [Default: deidentified_sam]
-obams (--outdir_bamsnap)   Directory name for bamsnap output [Default: bamsnap]

# Samtools view options
-ref (--reference_name)     Name of reference gene used for mapping [Default: BaxD2]
-r1 (--range1)              Start nucleotide position for SAM window [Default: 30]
-r2 (--range2)              End nucleotide position for SAM window [Default: 40]
-q (--MapQ)                 Filter bam file by uniqueness of alignment [Default: 10]
-F (--excl_flag)            Exclude alignments with specified samtools flag [Default: 256]   
-t (--threads)              Number of threads to use with samtools [Default: 16]

# Bamsnap (https://github.com/parklab/bamsnap)
-bamsnap                    Generate alignment images with BAMSNAP
-fa (--fasta)               Required gene sequence in fasta format
-pos (--position_view)      Gene nucleotide position to view [Default: 35]
-m (--margin)               Genomic margin size [Default: 15]

USAGE
die "$usage\n" unless @ARGV;

## Defining variables

## Required input
my @bam;                    # List of input BAM files
my $clinical;               # ROSMAP clinical file
my $biospecimen;            # ROSMAP biospecimen file
my @ref_name;               # Ref gene name
my $start = 1;              # Start numbering patients at

## Output options
my $table_out = 'out.tsv';  # output table name in tsv format
my $outdir_sam = 'deidentified_sam';
my $outdir_bamsnap = 'bamsnap';

## Samtools specific options
my $range1 = 30;            # Start nucleotide position to filter SAM and BAM
my $range2 = 40;            # End nucleotide position to filter SAM and BAM
my $quality_score = 10;     # Quality score to filter reads in BAM file
my $flag = 256;             # Define samtools flag to exclude alignments from filtered files
my $threads = 16;           # Number of threads for samtools

## Only for bamsnap analysis
my $bamsnap;                # Flag to use bamsnap analysis
my $ref_fasta;              # Fasta of reference gene used for previous mapping
my $bamsnap_pos = 35;       # Gene nucleotide position to view
my $margin = 15;            # Genomic margin size

GetOptions(
    'bam=s@{1,}' => \@bam,
    'c|clinical=s' => \$clinical,
    'b|biospecimen=s' => \$biospecimen,
    
    'n|num=i' => \$start,
    'o|out=s' => \$table_out,
    'osam|outdir_sam=s' => \$outdir_sam,
    'obams|outdir_bamsnap=s' => \$outdir_bamsnap,
    
    'ref|reference_name=s{1,}' => \@ref_name,
    'r1|range1=i' => \$range1,
    'r2|range2=i' => \$range2,
    'q|MapQ=i' => \$quality_score,
    'F|excl_flag=i' => \$flag,
    't|threads=i' => \$threads,
    
    'bamsnap' => \$bamsnap,
    'fa|fasta=s' => \$ref_fasta,
    'pos|position_view=i' => \$bamsnap_pos
);
die "References not specified!\n" unless @ref_name;

## Generating SAM and BAM files for the specified reference range
## Setting up directory name for new SAM and BAM files

my @patient_name;

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
        my $name = fileparse($bam_file);
    if ($name =~ /^(\S+?)\.fastq/){
        my $info = $1;
        push (@patient_name, $info);
    }

    ## Generating SAM and BAM output file names
    system "echo Working on $bam_file\n";
    
    foreach(@ref_name){
        my $ref = $_;
        
        my $sam_dir = "${ref}.r${range1}-${range2}.SAM";
        my $bam_dir = "${ref}.r${range1}-${range2}.BAM";
        
        my ($sam_name) = fileparse($bam_file);
        $sam_name =~ s/.bam/.${ref}.filtered.sam/;
        my ($bam_name) = fileparse($bam_file);
        $bam_name =~ s/.bam/.${ref}.filtered.bam/;

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
                    -@ $threads \\
                    $bam_file"
            ); 
        }

        system (
            "samtools view \\
                -h \\
                -q $quality_score \\
                -F $flag \\
                -h \\
                -@ $threads \\
                $bam_file \\
                ${ref}:${range1}-${range2} \\
                > $filtered_sam"
        );
        
        ## Creating directory to store filtered BAM files
        unless (-d $bam_dir){
            mkdir ($bam_dir, 0755) or die "Cannot create $bam_dir: $!\n";
        }
        
        ## Generating the BAM files mapped to a specific range of the ref gene
        $filtered_bam = "${path}/${bam_dir}/${bam_name}";

        # Generate a BAM file and BAI file (if bamsnap) for the specified gene range
        system (
            "samtools view \\
                -h \\
                -q $quality_score \\
                -F $flag \\
                -@ $threads \\
                $bam_file \\
                ${ref}:${range1}-${range2} \\
                -bo $filtered_bam"
        );

        if($bamsnap){
            $index_bam = "${path}/${bam_dir}/${bam_name}.bai";
            
            system (
                "samtools index \\
                    -@ $threads \\
                    $filtered_bam \\
                    $index_bam"
            );
        }
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
);
# Adding references to the title
push (@out_columns, @ref_name);
# Adding reference % column to the title
foreach(@ref_name){
    push(@out_columns, "$_".'%');
}
# Printing title in output
print TSV "\#Patient";
foreach (@out_columns){ print TSV "\t$_"; }
print TSV "\n";

## Deidentifying patient data name, start patient numbering from defined in command line 
my $patient_count = $start - 1;

## Creating directory for deidentified SAM files
my $sam_deid_dir = $outdir_sam;

unless (-d $sam_deid_dir){ 
    mkdir ($sam_deid_dir, 0755) or die "Cannot create $sam_deid_dir: $!\n";
}

###########################################################
### Working on filtered SAM, deidentifying data and BAMSNAP
###########################################################

my %db;

## Working on filtered SAM files
for my $data_name (@patient_name){
    $patient_count++;
    $patient_count = sprintf("%05d", $patient_count);
    print "####### patient_$patient_count\n";
    
    foreach(@filt_sam){
        my $sam = $_;
        if ($sam =~ /$data_name/){
            my ($gene) = $sam =~ /.*\.(\S+?).filtered.sam/;
            my $sam_deid = "patient_${patient_count}_${gene}_deidentified.sam";           
            my $sam_dir = "${gene}.r${range1}-${range2}.SAM";
            my $bam_dir = "${gene}.r${range1}-${range2}.BAM";

            my $path_sam = "${path}/${sam_dir}/$sam";
            my $path_deid = "${path}/${sam_deid_dir}/$sam_deid";
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
            my $sam_ref;
            my $start_align;
            my $end_align;
            my $cigar;
            my $sam_flag;
            my $seq;

            my $match_pos;

            while (my $line = <SAM>){

                chomp $line;

                ## Keeping SAM file header
                if ($line =~ /^\@HD/){ print DEID "$line\n"; }
                
                ## Keeping reference alignment sequence info
                elsif ($line =~ /^\@SQ/){ print DEID "$line\n"; }
                
                ## Getting the biospecimen id from the sam file and deidentifiying data name
                elsif ($line =~ /^\@RG\s+ID:(\S+)\s+(\S+)/){
                    my $id = $1;                                        # dataset name
                    $SM = $2;                                           # SM is reference mapped to
                    my $sample_id = fileparse($id);                     # grabbing only file name
                    ($sID) = $sample_id =~ /^(\w+)\./;                  # parsing data name to get biospecimen ID
                    $sID =~ s/_S\d+_R[12]_\d+$//;                       # parsing file name to match biospecimen ID formatting
                    my $id_len = "patient_${patient_count}.fastq.gz";   # deidentifying dataset name
                    print DEID "\@RG\tID:$id_len\t$SM\n";               # printing out to deidentified filtered SAM file
                }
                ## Skipping command line details
                elsif ($line =~ /^\@PG/){
                    
                    # Found some \t errors on the SAM header, fixing it
                    $line =~ s/\\t/\t/g;
                    
                    # replacing all occurances of patient data name for deidentified patient count
                    my $pos = 0;
                    while ($pos < length($line)) {
                        $pos = index($line, $sID, $pos);
                        if ($pos == -1){
                            last;
                        }
                        my $pat = "patient_${patient_count}";
                        substr($line, $pos, length($sID), $pat);
                        $pos = length($pat);
                    }
                    print DEID "$line\n";
                }
                else {
                    $read_count++;
                    $read_count = sprintf("%02d", $read_count);

                    ## Deidentifying lines and grabbin cigar info
                    my @columns = split("\t", $line);
                    ($read_id) = $columns[0] = "deidentified_read_$read_count";                 # Removing read ID
                    $sam_flag = $columns[1];                                                    # Get flag to determine read mapping orientation
                    $sam_ref = $columns[2];                                                     # Reference gene
                    $start_align = $columns[3];                                                 # start position of alignment
                    $end_align = $columns[4];                                                   # end position of alignment
                    $cigar = $columns[5];                                                       # Grabbing CIGAR info
                    $seq = $columns[9];                                                         # Getting DNA sequence
                    my $seq_length = 'N' x length($seq);                                        # Getting length of DNA seq and replacing by N
                    $columns[9] = $seq_length;                                                  # Replacing DNA seq on file
                    my $dataset = "RG:Z:patient_${patient_count}.fastq.gz";                     # Deidentifying dataset name
                    $columns[11] = $dataset;                                                    # Replacing dataset name
                    
                    # Printing deidentified SAM file
                    for (0..$#columns - 1){ print DEID "$columns[$_]\t"; }
                    print DEID "$columns[$#columns]\n";
                }
            }

            ## Getting clinical data for each patient based on biospecimen and individual IDs
            my $in_ID = $biospecimen{$sID};
            
            # Checking if the data is also present in clinical metadata
            if (exists $clinic{$in_ID}){
                $db{$patient_count}{$in_ID}{$gene} = $read_count;
            }
            
            ## Running bamsnap if required
            if ($bamsnap){
                
                # Checking for the required input fasta file
                die "ERROR: Cannot run bansnap without reference file! $!" unless defined $ref_fasta;
                
                # Creating bamsnap directory
                my $bamsnap_dir = $outdir_bamsnap;

                unless (-d $bamsnap_dir){
                    mkdir ($bamsnap_dir, 0755) or die "Cannot create $bamsnap_dir: $!\n";
                }
                            
                # Generating bamsnap files
                my $bamsnap_out = "${path}/${bamsnap_dir}/patient_${patient_count}_${gene}_pos${bamsnap_pos}_m${margin}.png";
                system (
                    "bamsnap \\
                    -bam $filtered_bam \\
                    -ref $ref_fasta \\
                    -pos ${gene}:${bamsnap_pos} \\
                    -out $bamsnap_out \\
                    -title patient_${patient_count}_${gene}_pos${bamsnap_pos}_m${margin} \\
                    -draw coordinates bamplot \\
                    -read_group strand \\
                    -margin $margin"
                );
            }
        }
    }
}


## Working on ouput result TSV table
for my $patient (sort (keys %db)){
    my $total = 0;
    for my $ids (keys %{ $db{$patient}}){
        print TSV "$patient\t";

        # Calculating total number of reads per patient and adding the % key to the database
        foreach (@ref_name){
            if (exists $db{$patient}{$ids}{$_}){
                my $value = $db{$patient}{$ids}{$_};
                $total += $value;
            }
            my $perc = "$_"."%";
            $db{$patient}{$ids}{"$perc"}=undef;
        }

        # Determining the percentages calculations
        for my $keys (keys %{$db{$patient}{$ids}}){
            if ($keys =~ /(\S+)%$/){
                my $gene_name = $1;
                my $count = $db{$patient}{$ids}{$gene_name};
                
                # For no reads mapping to ref, print values 0
                if ($count == 0){
                    $db{$patient}{$ids}{$keys} = '0';
                }
                # For values other than 0, calculate percentage
                else{
                    my $calc = ($count/$total)*100;
                    $calc = sprintf("%.2f", $calc);
                    $db{$patient}{$ids}{$keys} = $calc;
                }
            }
        }

        # Printing output table info 
        foreach(@out_columns){
            
            # Printing metadata
            if (exists $clinic{$ids}{$_}){
                print TSV "$clinic{$ids}{$_}\t";
            }

            # Printing reads count and calculations for reference genes
            if (exists $db{$patient}{$ids}{$_}){
                my $value = $db{$patient}{$ids}{$_};
                print TSV "$value\t";
            }
        }
        print TSV "\n";
    }
}