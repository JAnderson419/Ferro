#!/usr/local/bin/perl

# Jackson Anderson
# Spring 2015 - Senior Design
# Microelectronic Engineering - Rochester Institute of Technology
# 
# For use parsing dynamic hysteresis measurements made on the TF Analyzer 1000
# Reads .dat files and outputs a folder named after the .dat file 
# with seperate .tsv files for each measurement table
#
# v2 Release Notes:
# Removes []+-/ from headers since matlab confuses them on tdfread
#
# v3 Release Notes:
# Karine: Added table number as unique measurement identifier
#
# v4
# Added initial support for fatigue & PUND measurements
# Removed $number - same thing as $table

use strict;
use warnings;
use Cwd;


my $fileName = '';
my $summaryTableRead = 0;
my $isDataTable = 0;
my $isDataTableHeader = 1;
my $firstline = 1;
my ($table, $sampleName, $frequency, $amplitude, $average, $cycles, $type);

my $workingDir = getcwd;
opendir(DIR, $workingDir) or die;
my @dir = readdir DIR;
 
foreach my $item (@dir) {
   print "$item\n";
}

closedir DIR;

print "Enter the .dat filename (with extension) you would like to parse\n";
$fileName = <STDIN>;
chomp $fileName;
unless (-e $fileName) {
	print "File not found. Exiting.\n";
	exit;
}

my $dirName = substr ($fileName,0,length($fileName)-4);

unless(-e $dirName or mkdir $dirName) {
        die "Unable to create $fileName\n";
    }

open DATA, "<", $fileName;

while (<DATA>) {
	chomp;
	if ($firstline){
		$type = 0 if /DynamicHysteresisResult/;
		$type = 1 if /Fatigue/;
		$type = 2 if /PulseResult/;
		$firstline = 0;
	}
	s/\t$//;
	#Skip measurement summary table
	$summaryTableRead = 1 if /^DynamicHysteresis$/ or /^Data Table/ or /^Pulse$/;
	next unless $summaryTableRead;



	$table = $1 if /Table (\[*\d*.*\d*\]*)/;
	$sampleName = $1 if /SampleName: (.*)/;
	if ($type != 2) {
		$frequency = $1 if /Hysteresis Frequency \[Hz\]: (\d*)/;
		$amplitude = $1 if /Hysteresis Amplitude \[V\]: (\d*)/;
	}
	elsif ($type == 2) {
		$frequency = $1 if /Pund Frequency \[Hz\]: (\d*)/;
		$amplitude = $1 if /Pund Amplitude \[V\]: (\d*)/;
	}
	#$amplitude=~ s/(\.\d{2})\d*/$1/;
	$average = $1 if /Averages: (.*)/;
	$cycles = $1 if /Total Cycles: (\d*)/;
	
	$isDataTable = 1 if /^Time \[s\]/;
	if ($isDataTable){
	    my $datapath;
		$datapath = "$dirName\\$sampleName $frequency".'Hz '."$amplitude".'V '."$average".'Average Table'."$table".'.tsv' if $type == 0;
	    $datapath = "$dirName\\$sampleName $frequency".'Hz '."$amplitude".'V '."$average".'Average Table'."$table $cycles".'Cycles.tsv' if $type == 1;
		$datapath = "$dirName\\$sampleName $frequency".'Hz '."$amplitude".'V Table'."$table".'.tsv' if $type == 2;
	
		unless(defined $datapath) {
			die "Measurement type not recognized";
		}
	
		$isDataTableHeader ? (open DATAFILE, ">" , $datapath) : (open DATAFILE, ">>",$datapath);
		if ($isDataTableHeader) {
			$isDataTableHeader = 0;
			s/\[//g;
			s/\]//g;
			s/\//_per_/g;
			s/\+/plus/g;
			s/-/minus/g;
		}
		print DATAFILE;
		print DATAFILE "\n";
		if (/^$/) {
			$isDataTable = 0;
			$isDataTableHeader = 1;
			close DATAFILE;
		}
	}
	

}

close DATA;

exit;


