#!/usr/bin/perl -w
use List::Util qw(max min);
#using strict;
# RNase accesibility to each RNA sites (footprinting Cscore)

die "Usage: $0 Wiggle \n" if($#ARGV!=0);
my $inputfile=shift;

open(WIGGLE, $inputfile)||die "$!";

my @data=<WIGGLE>;

foreach (@data)
{
	chomp;
	my @temp=split(/\t/,$_);
	my @rnase=();
	push (@rnase,$temp[2]);
	push (@rnase,$temp[3]);
	push (@rnase,$temp[4]);
	my $maxrnase = max(@rnase);
	my $cscore=log(($maxrnase+1)/($temp[1]+1))/log(10);

	print "$_\t";
    printf "%.8f\n",$cscore;
}

close(WIGGLE); 
