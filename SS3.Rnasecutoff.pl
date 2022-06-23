#!/usr/bin/perl -w

#using strict;
# Identified RNase Footprinting sites

die "Usage: $0 BED, cutoff \n" if($#ARGV!=1);
my $inputfile=shift;
my $threshold=shift;


open(WIGGLEBED, $inputfile)||die "$!";

my @data=<WIGGLEBED>;
my %Rnasefoot=();
my $i=0;

foreach (@data)
{
	chomp;
	my @temp=split(/\t/,$_);	
	my $sites=$temp[2];
	my $cscore=$temp[3];
	if ($cscore <= $threshold)
	{
		$Rnasefoot{$i}.=$sites.";";			
	}
	else
	{
		$i++;
	}
}
# print "Footprinting hash construction\n";

foreach my $keytag(sort keys %Rnasefoot)
{

	my @sites=split(/;/,$Rnasefoot{$keytag});
	my $footprintslength=@sites; # length of array and footprinting length;
	my $start=$sites[0]-1;  # footprinting start site
	my $end=$sites[$footprintslength-1]; # footprinting end site
	print "chrM\t$start\t$end\t$footprintslength\n";	
}

close(WIGGLEBED); 
