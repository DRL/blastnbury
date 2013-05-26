#!/usr/bin/perl
use strict;
use warnings;

##############################################################################
#
# File   :  	blastnbury.pl
# History:  	2013-04-23 (dominik) first implementation
# Coming soon:  [S]kip, [B]ack, [A]bort, [M]ore	
#
##############################################################################
#
#  Perl script for performing a remote BLAST of a query against nt, 
#  followed by user selection of sorted hits saved as CVS (OSIRIS import)
#
##############################################################################

=head1 

blast'n'bury - BLAST your sequences and bury them in a database

=head1 SYNOPSIS

B<blastnbury.pl> [B<-q>] [B<-noblast>] [B<-hit> B<-report> B<-save>]

=head1 DESCRIPTION

B<blastnbury> takes a fasta file of DNA sequences and BLASTs them against nt, 
followed by saving the each resulting BLAST report into a tab-separated file (.txt). 
The script then loops the files extracting relevant information, and presents them to 
the user. The user then decides which of the "best" hits are saved. The selected results 
of all sequences are then exported to a CSV file.

=head1 OPTIONS

=over 4

=item B<-q>
    
FASTA file of DNA sequence(s)

=item B<-noblast>
    
No remote BLAST is performed. For analysing existing BLAST reports.

=item B<-hit>

Number of BLAST hits per sequence that are saved in tab-delimited files (default: 100)

=item B<-report>

Number of "best" sorted Blast hits from which the user can choose (default: 15)

=item B<-save>

Number of "best" sorted Blast hits the user can choose for the CSV file (default: 5). Do not change w/o telling OSIRIS.

=item B<-man>
    
Prints the manual page and exits.

=back

=head1 DIAGNOSTICS

=over 4

=back

=head1 AUTHOR

Dominik R. Laetsch, dominik.laetsch@gmail.com

=cut

##############################################################################
#
#   								START
# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------
use Getopt::Long qw(:config pass_through no_ignore_case);
use Pod::Usage;
use Data::Dumper;
use Bio::SeqIO;
use Term::ANSIColor;
use POSIX;
# ----------------------------------------------------------------------------
#
##############################################################################
#				 Parsing BLAST results
##############################################################################
#
# For each of (%data => $seq_id) the respective *.CVS file is openend.
# Each line is split into @hit_data. The first three strings are $hit_desc, $hit_acc, and $hit_length, respectively (shifted from @hit_data).
#
# Since multiple HSPs per ACC are spread over several lines, population of fields is 
# guided by the existance of a defined ($data{$id}{'hits'}{$hit_acc}). If already defined, 
# then there are multiple HSPs.
#  
# Fields are populated generating the following structure:
#	
# 		%data 	=>	$seq_id	=>	'query'	=>	'seq'		=>	$query_object->seq
#										=>	'length'	=>	$query_object->length	
#							=>	'hits'	=>	$hit_acc	=>	'link'		=>	"http://www.ncbi.nlm.nih.gov/nuccore/".$hit_acc
#														=>	'length'	=>	$hit_length
#														=>	'desc'		=>	&TRUNCATE_DESC($hit_desc, $truncate_desc_to_length)
#														=> 	'main'		=>	'hsps'			=> $num_of_hsps + 1 (human numbering)
#														=> 				=> 	'alen'			=>	
#																		=> 	'mismatch'		=>	
#																		=> 	'bitscore'		=>
#																		=> 	'gaps'			=>		
#																		=> 	'nident'		=>	
#														=>	'hsps'		=>	$num_of_hsps	=> 'pident' 	=>
#																							=> 'alen' 		=>
#																							=> 'qstart' 	=>
#																							=> 'qend' 		=>
#																							=> 'sstart' 	=>
#																							=> 'send' 		=>
#																							=> 'frames' 	=>
#																							=> 'bitscore' 	=>
#																							=> 'gaps' 		=>
#																							=> 'nident' 	=>
#																							=> 'eval' 		=>
# 
##############################################################################


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

my ($query_file) = ("");
my ($blast_switch,$hit_threshold, $display_threshold, $report_threshold, $help, $man) = (1,100,15,5,"",""); 

# Get options
# ----------------------------------------------------------------------------
my @parameters=("-q","-blast","-noblast","-hit","-hits","-display","-save"); #&VALIDATE_ARGS uses this for validation of parameters

GetOptions (
	"q:s" => \$query_file,
 	"blast!" => \$blast_switch, # -blast || -noblast
  	"hit|hits:i" => \$hit_threshold,
  	"display:i" => \$display_threshold,
  	"save:i" => \$report_threshold,
  	"help" => \$help,
  	"man" => \$man,
) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if ($man);
# ============================================================================

# Sanitising user input
# ----------------------------------------------------------------------------
die pod2usage(2) unless (&VALIDATE_ARGS(\@ARGV, \@parameters));
die pod2usage(1) unless $query_file;
unless ($hit_threshold > $display_threshold && $display_threshold > $report_threshold) {
	die "\n\tRule: Number of -hits > Number of -display > Number of -save\n\t   => e.g.: blastnbury.pl -q test.fas -blast -hits 100 -display 15 -save 5 (default)\n\n";
}
# ============================================================================

# Welcome screen
# ----------------------------------------------------------------------------
&PRINT_LOGO('help');
# ----------------------------------------------------------------------------

my %data;
my $truncate_desc_to_length=50;
# Read file
# ----------------------------------------------------------------------------
my ($query_file_object, $query_object, $seq_id);
my (@query_objects, @query_ids, %query_seqs, %query_lengths);
my $datestamp;
my $start_of_read=time();
# The $query_file is read using Bio::Seq.
$query_file_object = Bio::SeqIO->new(
	'-format' => 'fasta', 
	'-file' => "<$query_file"
);
# The first two fields of %data are populated giving it the following structure:
#
# 		%data 	=>	$seq_id	=>	'query'	=>	'seq'	=>	$query_object->seq
#										=>	'length'=>	$query_object->length	
#
my $number_of_sequences=0;
while ($query_object = $query_file_object->next_seq) {
	if(defined $query_object){
		$number_of_sequences++;
		$seq_id= $query_object->id;
		#$data{$seq_id}{'obj'}=$query_object;
		$data{$seq_id}{'id'}=$query_object->id;
		$data{$seq_id}{'query'}{'seq'}=$query_object->seq;
		$data{$seq_id}{'query'}{'length'}=$query_object->length;
    }
}

print strftime("[%H:%M:%S]",localtime)." - Read ".$number_of_sequences." sequence(s) [".&TIME_DIFFERENCE($start_of_read)."sec]\n";
# ============================================================================

# BLAST each sequence unless -noblast
# ----------------------------------------------------------------------------
if ($blast_switch){
	# If the program was called with -blast a BLAST search is performed
	my ($query_seq, $start_of_blast, $number_of_hits);
	my $i=1;
	
	print strftime("[%H:%M:%S]",localtime)." - Starting BLAST search for ... ".$hit_threshold." hits\n";
	$start_of_blast=time();
	# For each of the $query_id the &BLAST subroutine is called by passing:
	#
	# 		(%data => $seq_id) and (%data => $seq_id =>	'query' => 'seq' => $query_seq)
	# 
	foreach my $query_id (sort keys %data){
		$query_seq = $data{$query_id}{'query'}{'seq'};
		&BLAST($query_id, $query_seq);
		print strftime("[%H:%M:%S]",localtime)." = ".$query_id." (".$i." of ".$number_of_sequences.") [".&TIME_DIFFERENCE($start_of_blast)."sec] \n"; 
		$start_of_blast=time();
		$i++;
	}
}
else{
	print strftime("[%H:%M:%S]",localtime)." - Skipping BLAST search. Trying to find '*_BLAST.txt' file(s)...\n"; 
}
# ============================================================================

# Parsing BLAST output
# ----------------------------------------------------------------------------
# Check for existance of file(s)
my $die_flag_file=0;
my $start_check=time();
foreach my $id (sort keys %data){
	my $result_file = $id."_BLAST.txt";
	if (&CHECK_FILE($result_file)){
		next;
	}
	else{
		print "   [ERROR] - ".$result_file." could not be found\n";
		$die_flag_file=1;
	}
} 
if ($die_flag_file){
	die "   [ERROR] - Please restart with [-blast]\n";
}
unless ($blast_switch){
	print strftime("[%H:%M:%S]",localtime)." - Found ".$number_of_sequences." file(s) [".&TIME_DIFFERENCE($start_check)."sec]\n";
}
print strftime("[%H:%M:%S]",localtime)." - Parsing BLAST results...\n";
my %hash_for_sorting;
my $sequence_counter = 0;

foreach my $id (sort keys %data){
	$sequence_counter++;
	my $hit_counter = 0; # 
	my $num_of_hsps = 0;
	my $start_of_get_desc = time();
	my @fields = ("bitscore", "qstart", "qend", "sstart", "send", "eval", "alen", "nident", "pident", "mismatch", "gaps", "frames"); # fields of the CSV file (determined by &BLAST)
 	my $blast_file = $id."_BLAST.txt";

 	open (BLAST_HITS, "<".$blast_file) or die "   [ERROR] - ".$blast_file." could not be opened\n";
 	while (<BLAST_HITS>){
 		chomp;
 		if (m/^#/){
 			$datestamp=$_;
 			$datestamp=~ s/#//;
 		}
 		else {
 			if ($hit_counter < $hit_threshold){
 				my @hit_data=split/\t/;
	# 			print Dumper(\@hit_data);
 				my $hit_desc=shift @hit_data;
 				my $hit_acc=shift @hit_data;
 				my $hit_length=shift @hit_data;
	    		# Checks whether an ACC has already been read before from the files (i.e. multiple HSPs per ACC)
 				if (defined $data{$id}{'hits'}{$hit_acc}){ 
 					$num_of_hsps++;
 					$data{$id}{'hits'}{$hit_acc}{'main'}{'hsps'}=$num_of_hsps+1;
 					$data{$id}{'hits'}{$hit_acc}{"main"}{'eval'}=1;
 					$data{$id}{'hits'}{$hit_acc}{"main"}{'qstart'}="N/A";
 					$data{$id}{'hits'}{$hit_acc}{"main"}{'qend'}="N/A";
 					$data{$id}{'hits'}{$hit_acc}{"main"}{'sstart'}="N/A";
 					$data{$id}{'hits'}{$hit_acc}{"main"}{'send'}="N/A";
					for (my $k=0; $k<scalar(@fields); $k++){
						if ($fields[$k] eq "bitscore" || $fields[$k] eq "alen" || $fields[$k] eq "nident" || $fields[$k] eq "mismatch" || $fields[$k] eq "gaps"){ # values in these fields can be added when multiple HSPs
							$data{$id}{'hits'}{$hit_acc}{'main'}{$fields[$k]}=$data{$id}{'hits'}{$hit_acc}{'main'}{$fields[$k]} + $hit_data[$k];
							$data{$id}{'hits'}{$hit_acc}{'hsps'}{$num_of_hsps}{$fields[$k]}=$hit_data[$k];
						}
 						if ($fields[$k] eq "sstart" || $fields[$k] eq "send" || $fields[$k] eq "qstart" || $fields[$k] eq "qend" || $fields[$k] eq "eval" || $fields[$k] eq "pident" || $fields[$k] eq "frames"){ # values in these fields unique for each HSPs
 							$data{$id}{'hits'}{$hit_acc}{'hsps'}{$num_of_hsps}{$fields[$k]}=$hit_data[$k];
 						}
 						if ($fields[$k] eq "pident"){
 							$data{$id}{'hits'}{$hit_acc}{'hsps'}{$num_of_hsps}{$fields[$k]}=$hit_data[$k];
 							$data{$id}{'hits'}{$hit_acc}{'main'}{$fields[$k]} = 100*$data{$id}{'hits'}{$hit_acc}{'main'}{'nident'}/$data{$id}{'hits'}{$hit_acc}{'main'}{'alen'};
 						}
					}
				}
				# For first hit for an ACC (first hit of multiple HSPs, or single HSP)
				else{
					$num_of_hsps = 0;
					$data{$id}{'hits'}{$hit_acc}{'acc'} = $hit_acc;
					$data{$id}{'hits'}{$hit_acc}{'main'}{'hsps'} = $num_of_hsps + 1;
					$data{$id}{'hits'}{$hit_acc}{'main'}{'qlen'} = $data{$id}{'query'}{'length'};
					$data{$id}{'hits'}{$hit_acc}{'link'}="http://www.ncbi.nlm.nih.gov/nuccore/".$hit_acc;
					$data{$id}{'hits'}{$hit_acc}{'length'}=$hit_length;
 					$data{$id}{'hits'}{$hit_acc}{'desc'}=&TRUNCATE_DESC($hit_desc, $truncate_desc_to_length);

 					$hit_counter++; # sorting counter
 					$hash_for_sorting{$id}{$hit_counter}=$hit_acc;

 					for (my $k=0; $k<scalar(@fields); $k++){
 						#print $fields[$k]." => ".$hit_data[$k]."\n"; 
 						if ($fields[$k] eq "pident"){
 							$data{$id}{'hits'}{$hit_acc}{'hsps'}{$num_of_hsps}{$fields[$k]}=$hit_data[$k];
 							$data{$id}{'hits'}{$hit_acc}{'main'}{$fields[$k]} = 100*$data{$id}{'hits'}{$hit_acc}{'main'}{'nident'}/$data{$id}{'hits'}{$hit_acc}{'main'}{'alen'};
 						}
						else{
							$data{$id}{'hits'}{$hit_acc}{"main"}{$fields[$k]}=$hit_data[$k];
							$data{$id}{'hits'}{$hit_acc}{"hsps"}{$num_of_hsps}{$fields[$k]}=$hit_data[$k];
						} 	
					}
				}
			$data{$id}{'hits'}{$hit_acc}{'main'}{'differences'} = $data{$id}{'hits'}{$hit_acc}{'main'}{'mismatch'} + $data{$id}{'hits'}{$hit_acc}{'main'}{'gaps'};	
			$data{$id}{'hits'}{$hit_acc}{'main'}{'cov'} = 100*($data{$id}{'hits'}{$hit_acc}{'main'}{'alen'}-$data{$id}{'hits'}{$hit_acc}{'main'}{'gaps'})/$data{$id}{'hits'}{$hit_acc}{'main'}{'qlen'};
			}		
		}
	}
	close (BLAST_HITS);
	# Last things
	#print Dumper(\%data);
 	print strftime("[%H:%M:%S]",localtime)." = ".$id." (".$sequence_counter." of ".$number_of_sequences.") [".&TIME_DIFFERENCE($start_of_get_desc)."sec] \n";
}

# ============================================================================
# Display output
# ----------------------------------------------------------------------------

# Header for CSV file
open (CSV_FILE, ">".$query_file."_results.csv") || die "Cannot write to ".$query_file."_results.csv\n";
print CSV_FILE "#BEA_CODE,length,seq,date";
for (my $j = 1; $j < $report_threshold+1; $j++){
	print CSV_FILE ",HIT_ACC_".$j.",HIT_DESC_".$j.",NUM_OF_HSPS_".$j.",EVAL_".$j.",ID_".$j.",COV_".$j.",SITES_".$j.",NUCDIFF_".$j;
}
print CSV_FILE "\n";
close (CSV_FILE);

foreach my $id ( sort keys %hash_for_sorting){
	&PRINT_RESULT_HEADER($id, $data{$id}{'query'}{'length'});
	my $display_counter=1;
	foreach my $hit_number (sort { $a <=> $b} keys %{$hash_for_sorting{$id}}){
		if ($display_counter <= $display_threshold){
			my $acc=$hash_for_sorting{$id}{$hit_number};
			my %hit_data=%{$data{$id}{'hits'}{$acc}};
			my $query_length=$data{$id}{'query'}{'length'};
			&PRINT_HIT($query_length, $hit_number, %hit_data); 	
		}
		$display_counter++;
	}
	my @user_input=&ASK_FOR_USER_INPUT($display_threshold, $report_threshold);
	my $sort_order;
	while (scalar(@user_input) == 1){
		my %hash_of_fields;
		if ($user_input[0] eq 'b'){
			foreach my $position (keys %{$hash_for_sorting{$id}}){
				my $current_acc = $hash_for_sorting{$id}{$position};
				$hash_of_fields{$current_acc}=$data{$id}{'hits'}{$current_acc}{'main'}{'bitscore'};
			}
			$sort_order="descending";
		}
		elsif ($user_input[0] eq 'p'){
			foreach my $position (keys %{$hash_for_sorting{$id}}){
				my $current_acc = $hash_for_sorting{$id}{$position};
				$hash_of_fields{$current_acc}=$data{$id}{'hits'}{$current_acc}{'main'}{'pident'};
			}
			$sort_order="descending"; 
		}
		elsif ($user_input[0] eq 'c'){
			foreach my $position (keys %{$hash_for_sorting{$id}}){
				my $current_acc = $hash_for_sorting{$id}{$position};
				$hash_of_fields{$current_acc}=$data{$id}{'hits'}{$current_acc}{'main'}{'cov'};
			}
			$sort_order="descending";
		}
		elsif ($user_input[0] eq 'd'){
			foreach my $position (keys %{$hash_for_sorting{$id}}){
				my $current_acc = $hash_for_sorting{$id}{$position};
				$hash_of_fields{$current_acc}=$data{$id}{'hits'}{$current_acc}{'main'}{'differences'};
			}
			$sort_order="ascending";
		}
		elsif ($user_input[0] eq 'a'){
			foreach my $position (keys %{$hash_for_sorting{$id}}){
				my $current_acc = $hash_for_sorting{$id}{$position};
				$hash_of_fields{$current_acc}=$data{$id}{'hits'}{$current_acc}{'main'}{'alen'};
			}
			$sort_order="descending";
		}
		elsif ($user_input[0] eq 'e'){
			foreach my $position (keys %{$hash_for_sorting{$id}}){
				my $current_acc = $hash_for_sorting{$id}{$position};
				$hash_of_fields{$current_acc}=$data{$id}{'hits'}{$current_acc}{'main'}{'eval'};
				$sort_order="ascending";
			}
		}
		my @array_of_values;
		foreach my $values (values %hash_of_fields){
			push @array_of_values, $values;
		}
		my @sorted_array_of_values=&SORT_BY(\$sort_order, \@array_of_values);
		my $rank=1;
		foreach my $value_of_sorted_array_of_values (@sorted_array_of_values){
			foreach my $acc_of_hash_of_fields (keys %hash_of_fields){
				if ($value_of_sorted_array_of_values ne "N/A"){
					if ($value_of_sorted_array_of_values == $hash_of_fields{$acc_of_hash_of_fields}){
						$hash_for_sorting{$id}{$rank}=$acc_of_hash_of_fields;
						delete($hash_of_fields{$acc_of_hash_of_fields});
						$rank++;
					}
				}
			}
		}
		my $display_counter=1;
		#&PRINT_LOGO;
		&PRINT_RESULT_HEADER($id, $data{$id}{'query'}{'length'});
		foreach my $sorted_hit_number (sort {$a <=> $b} keys %{$hash_for_sorting{$id}}){
			if ($display_counter <= $display_threshold){
				my $acc=$hash_for_sorting{$id}{$sorted_hit_number};
				my %hit_data=%{$data{$id}{'hits'}{$acc}};
				my $query_length=$data{$id}{'query'}{'length'};
				&PRINT_HIT($query_length, $sorted_hit_number, %hit_data); 
				$display_counter++;
			}	
		}
		@user_input=&ASK_FOR_USER_INPUT($display_threshold, $report_threshold);
	}

	if (scalar(@user_input) == $report_threshold){
		open (CSV_FILE, ">>".$query_file."_results.csv") || die "Cannot write to ".$query_file."_results.csv\n";
		print CSV_FILE &PRINT_SEQ_DATA_TO_FILE(\$datestamp, \%{$data{$id}});
		close (CSV_FILE);
		my %hash_of_selected_hits;
		for (my $m=0; $m<$report_threshold; $m++){
			my $number_of_selected_hit = $user_input[$m];
			my $accession_of_selected_hit = $hash_for_sorting{$id}{$number_of_selected_hit};
			open (CSV_FILE, ">>".$query_file."_results.csv") || die "Cannot write to ".$query_file."_results.csv\n";
			print CSV_FILE &PRINT_HIT_TO_FILE(\%{$data{$id}{'hits'}{$accession_of_selected_hit}});
			close (CSV_FILE);
		}
		open (CSV_FILE, ">>".$query_file."_results.csv") || die "Cannot write to ".$query_file."_results.csv\n";
		print CSV_FILE "\n";
		close (CSV_FILE);
	}
}


# ============================================================================
# Subroutines
# ============================================================================

# ----------------------------------------------------------------------------
# &PRINT_LOGO subroutine
# ----------------------------------------------------------------------------

sub PRINT_LOGO{
	system("clear");
	my $additional_string='';
	if ($_[0 eq 'help']){
		$additional_string = "\t\t###      (more under blastnbury.pl -help -man)   ###\n";
	}
	print ("\n\t\t####################################################\n");
	print ("\t\t###                                              ###\n");
	print ("\t\t###           blast'n'bury Version 0.5           ###\n");
	print ("\t\t###      Interface for remote BLASTn searches    ###\n");
	print $additional_string;
	print ("\t\t###                                              ###\n");
	print ("\t\t####################################################\n");
	print "\n";
}
# ----------------------------------------------------------------------------
# &VALIDATE_ARGS subroutine
# ----------------------------------------------------------------------------
sub VALIDATE_ARGS{
	my @arguments= @{$_[0]};
	my @parameters= @{$_[1]};
	my %parameters;
	my $return;
	foreach my $i (@parameters){
		$parameters{$i}='';

	}
	foreach my $argument (@arguments){
		$argument=lc($argument);
		unless (($argument =~ m/^-/ && defined $parameters{$argument}) || ($argument !~ m/^-/ && &CHECK_FILE($argument)) || ($argument !~ m/^-/ && $argument =~ /^[+-]?\d+$/)){
			return 0;
		}
	}
	return 1;
}
# ----------------------------------------------------------------------------
# &CHECK_FILE subroutine
# ----------------------------------------------------------------------------
sub CHECK_FILE{
	if (-e $_[0]) { #if exists return 1 else 0
 		return 1;
 	}
 	else{
 		return 0;
 	}
} 
# ----------------------------------------------------------------------------
# &BLAST subroutine
# ----------------------------------------------------------------------------
sub BLAST{
	# Building query
	my ($query_id, $query, $temp_file);
	$query_id=$_[0];
	$query=">".$_[0]."\n".$_[1];
	$temp_file=$_[0]."_BLAST.tmp";
	# Printing query to TEMP file
	open (TEMP, ">$temp_file") or die $!;
	print TEMP $query;
	close (TEMP);
	# System BLAST call
	my $outfile=$query_id."_BLAST.txt";
	$datestamp = strftime("%Y-%m-%d",localtime);
	my $blast_call = "blastn -query ".$temp_file." -outfmt '6 stitle sacc slen bitscore qstart qend sstart send evalue length nident pident mismatch gaps frames' -remote -db nt -max_target_seqs ".$hit_threshold." -out ".$outfile; 
	system($blast_call);
	### Prepend $datestamp to beginning of BLASTFILE
	open (BLASTFILE, "<".$outfile);
	my @outfile_data = <BLASTFILE>;	
	close (BLASTFILE);
	open (BLASTFILE, ">".$outfile);
	print BLASTFILE "#".$datestamp."\n";
	print BLASTFILE @outfile_data;
	close (BLASTFILE);
	#print $error;
 	 # Deleting TEMP file
    unlink($temp_file);
}

# ----------------------------------------------------------------------------
# &TRIM subroutine
# ----------------------------------------------------------------------------
sub TRIM{
	my $string = $_[0];
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

# ----------------------------------------------------------------------------
# &SORT_RESULTS subroutine
# ----------------------------------------------------------------------------

sub SORT_BY{	
	my $sub_sort_order=${shift(@_)}; 
	my @sub_array_of_values = @{shift(@_)};
	if ($sub_sort_order eq "ascending"){
		@sub_array_of_values = sort { $a <=> $b } @sub_array_of_values;
	}
	elsif($sub_sort_order eq "descending"){
		@sub_array_of_values = sort { $b <=> $a } @sub_array_of_values;
	}
	return @sub_array_of_values;
}
# ----------------------------------------------------------------------------
# &TIME_DIFFERENCE
# ----------------------------------------------------------------------------
sub TIME_DIFFERENCE{
	# curl call
	my $start=$_[0];
	my $end=time();
	my $diff=$end-$start;
	return $diff;
}

# ----------------------------------------------------------------------------
# &TRUNCATE_DESC subroutine
# ----------------------------------------------------------------------------
sub TRUNCATE_DESC{
	my $temp_desc = $_[0];
	my $truncate_desc_to_length = $_[1];
	# newlines
	chomp $temp_desc;
	# commas
	$temp_desc =~ s/[,]//g;
	#leading whitespace
	$temp_desc=~ s/^\s+//;
	#set string length to 75 characters
	if (length($temp_desc) < $truncate_desc_to_length){
		$temp_desc.=" "x($truncate_desc_to_length-length($temp_desc));
	}
	else{
		$temp_desc = substr($temp_desc, 0, $truncate_desc_to_length);
	}
	return $temp_desc;
}
# ----------------------------------------------------------------------------
# &PRINT_RESULT_HEADER subroutine
# ----------------------------------------------------------------------------

sub PRINT_RESULT_HEADER{
	my $print_id=$_[0];
	my $print_length=$_[1];
	
	print "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
    print color 'yellow';
	print "\t".$print_id." (".$print_length."nt)\n";
	print color 'reset';
	print "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	print "[ #] Description \t\t\t\t\tLength\tHSPs\tID\tCov\tEval\tqstart\tqend\tframe\t∆nt\tgaps\tURL/Comments\n";
	print "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
}

# ----------------------------------------------------------------------------
# &PRINT_HIT subroutine
# ----------------------------------------------------------------------------
sub PRINT_HIT{
	my($query_length, $hit_number, %hit_hash) = @_;
	my $hit_desc=$hit_hash{'desc'};
	my $hit_length=$hit_hash{'length'};
	my $hit_link=$hit_hash{'link'};
	my $num_of_hsps=$hit_hash{'main'}{'hsps'};
	my $hit_alen=$hit_hash{'main'}{'alen'};
	my $gaps=$hit_hash{'main'}{'gaps'};
	my $total_cov=sprintf("%.2f%%", 100*(($hit_alen-$gaps)/$query_length));
	#my $s_cov=sprintf("%.2f%%", 100*($hit_alen/$hit_length));
	my $total_nident=$hit_hash{'main'}{'nident'};
	my $total_diff="=> ∆nt = ".($hit_alen-$total_nident);
	my $pident=sprintf("%.2f%%",100*($total_nident/$hit_alen));
	if ($hit_number<10){
		print "[ ".$hit_number."] ";
	}
	else {
		print "[".$hit_number."] ";
	}
	print color 'magenta';
	print $hit_hash{'desc'}; #desc
	print color 'reset';
	print " ".$hit_length."nt \t";
	my $hsp_diff;
	foreach my $hsp (sort keys %{$hit_hash{'hsps'}}){
		###
		if ($hsp > 0){
			$num_of_hsps=' ';
			#$total_diff='';
			print " "x($truncate_desc_to_length)." \t\t";
		}
		my $hsp_pident=sprintf("%2s", sprintf("%.1f%%", $hit_hash{'hsps'}{$hsp}{'pident'}));
		my $hsp_gaps=$hit_hash{'hsps'}{$hsp}{'gaps'};
		my $hsp_cov=sprintf("%.1f%%", 100*(($hit_hash{'hsps'}{$hsp}{'alen'}-$hsp_gaps)/$query_length));
		my $hsp_eval=sprintf("%.0e", $hit_hash{'hsps'}{$hsp}{'eval'});
		my $hsp_qstart=$hit_hash{'hsps'}{$hsp}{'qstart'};
		my $hsp_qend=$hit_hash{'hsps'}{$hsp}{'qend'};
		my $hsp_sstart=$hit_hash{'hsps'}{$hsp}{'sstart'};
		my $hsp_send=$hit_hash{'hsps'}{$hsp}{'send'};
		my $hsp_frames = $hit_hash{'hsps'}{$hsp}{'frames'};
		$hsp_diff= $hit_hash{'hsps'}{$hsp}{'alen'} - $hit_hash{'hsps'}{$hsp}{'nident'};
		$hsp_frames=~ s/1/+/g;
		$hsp_frames=~ s/0/-/g;
		###
		printf("%4s\t%6s\t%6s\t%4s\t%4s\t%4s\t%3s\t%4s\t%4s\t", $num_of_hsps, $hsp_pident, $hsp_cov, $hsp_eval, $hsp_qstart, $hsp_qend, $hsp_frames, $hsp_diff, $hsp_gaps);
		if ($hsp+1 == $hit_hash{'main'}{'hsps'} && $num_of_hsps ne 1){
			print color 'bold';
			print $total_diff.", Gaps = ".$gaps." TotalCov: ".$total_cov."\n";
			print color 'reset'
		}
		if ($hsp+1 > 1 && $hsp+1 < $hit_hash{'main'}{'hsps'}){
			print "\n";
		}
		if ($hsp+1 == 1 && $num_of_hsps eq 1){
			print $hit_link."\n";
		}
		if ($hsp+1 == 1 && $num_of_hsps ne 1){
			print $hit_link."\n";
		}
	}
}
# ----------------------------------------------------------------------------
# &ASK_FOR_USER_INPUT subroutine
# ----------------------------------------------------------------------------
sub ASK_FOR_USER_INPUT{
	my $display_threshold= $_[0]; # not sure if starts with 0 or 1
	my $report_threshold=$_[1];
	my $user_input;
	my @input;
	my $input;
	my @answers = ('b', 'p', ,'c', 'd', 'a', 'e');
	my $red_flag = 1;
	my $message_one="\nSort (b)itscore, (p)ident, (c)ov, (d)ifferences, (a)len or (e)val ";
	my $message_two="or select ".$report_threshold." hits for the report : ";
	while ($red_flag){
		$red_flag = 0;
		print $message_one, $message_two;
		$user_input=<>;
		chomp($user_input);
		@input=split(/ /, $user_input);
		if (scalar(@input)==1){ #one element
			$input=shift @input;
			if ($input !~ m/\w{1}$/ || $input =~ m/\d{1}$/){
				$red_flag=1;
			}
		 	if (grep {/^$input$/} @answers){
		 		push @input, $input;
		 	}
		 	else {
		 		$red_flag=1;
		 	}

		}
		elsif (scalar(@input) == $report_threshold){ #array in length of report_threshold
			foreach my $number (@input) {
				if( $number <= 0 || $number > $display_threshold ) {
					print "\n   [ERROR] - ".$number." Not within range\n";
					$red_flag=1;
				}
			}
		}
		elsif (scalar(@input) < $report_threshold || scalar(@input) > $report_threshold ){ 
			$red_flag=1;
		}
	}
	return @input;
}

# ----------------------------------------------------------------------------
# &PRINT_HIT_TO_FILE subroutine
# ----------------------------------------------------------------------------
sub PRINT_SEQ_DATA_TO_FILE{
	my $date = ${shift(@_)}; 
	my %seq_data = %{shift(@_)};
	return $seq_data{'id'}.",".$seq_data{'query'}{'length'}.",",$seq_data{'query'}{'seq'}.",".$date;
}

# ----------------------------------------------------------------------------
# &PRINT_HIT_TO_FILE subroutine
# ----------------------------------------------------------------------------
sub PRINT_HIT_TO_FILE{
	my %temp_hit = %{shift()};
	my $temp_hit_acc = $temp_hit{'acc'};
	my $temp_hit_desc = $temp_hit{'desc'};
	my $temp_hit_num_hsps = $temp_hit{'main'}{'hsps'};
	my $temp_hit_eval;
	if ($temp_hit{'main'}{'eval'} != 1){
		$temp_hit_eval = $temp_hit{'main'}{'eval'};
	}
	else {
		$temp_hit_eval = "N/A";
	}
	my $temp_hit_ID = sprintf("%.2f",$temp_hit{'main'}{'pident'})."%";
	my $temp_hit_alen = $temp_hit{'main'}{'alen'};
	my $temp_hit_qlen = $temp_hit{'main'}{'qlen'};
	my $temp_hit_gaps = $temp_hit{'main'}{'gaps'};
	my $temp_hit_cov = sprintf("%.1f%%", 100*($temp_hit_alen - $temp_hit_gaps)/$temp_hit_qlen);
	my $temp_hit_sites = $temp_hit{'main'}{'nident'}."/".$temp_hit_alen; 
	my $temp_hit_nuc_diff = $temp_hit{'main'}{'differences'};
	return ",".$temp_hit_acc.",".$temp_hit_desc.",".$temp_hit_num_hsps.",".$temp_hit_eval.",".$temp_hit_ID.",".$temp_hit_cov.",".$temp_hit_sites.",".$temp_hit_nuc_diff;
}