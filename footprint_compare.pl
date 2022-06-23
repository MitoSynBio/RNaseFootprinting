#!/usr/bin/perl -w

#use strict;
#use warnings;

# For a set of footprints from the experimental sample, output the Cscore and Fscore of the equivalent coordinates from the control sample

die "Usage: $0 KO footprints (BED), WT Cscore (BED), lower range, upper range, flanking length\n" if($#ARGV!=4);

my $KOfootprint=shift;
my $WTcscore=shift;

my $lowrange=shift;
my $uprange=shift;
my $flanking=shift;

open(KO,$KOfootprint) || die "$!";
open(WT,$WTcscore) || die "$!";


my @data=<WT>;
my %Cscore=();

# Split each WT Cscore input line by tab and store the components in @temp; define a unique $Cscore variable for each position using stored items in @temp
foreach (@data)
{
	chomp;
	my @temp=split(/\t/,$_);
	$Cscore{$temp[2]}=$temp[3];
}

my @footprints=<KO>;
my @mergefootprints=();
my $mergestart=-1;
my $mergeend=-1;
my $mergelength=0;
#my $merge=0;

foreach (@footprints)
{
	chomp;
	my @temp=split(/\t/,$_);
	if($temp[1] == 0 )
	{
		$mergeend=$temp[2];
		$mergelength+=$temp[3];
	}
	elsif($temp[2] == 16299)
	{
		$mergestart=$temp[1]; 
		$mergelength+=$temp[3];
	}
	else
	{
		push(@mergefootprints,$_);
	}

}


# Correct for circular chromosome

if ($mergeend != -1 && $mergestart != -1 )
{
	$merge="chrM"."\t".$mergestart."\t".$mergeend."\t".$mergelength;
	push(@mergefootprints,$merge);
}

elsif($mergeend != -1 && $mergestart == -1)
{
	$merge="chrM"."\t"."0"."\t".$mergeend."\t".$mergelength;
	push(@mergefootprints,$merge);
}
elsif($mergeend == -1 && $mergestart != -1)
{
	$merge="chrM"."\t".$mergestart."\t"."16299"."\t".$mergelength;
	push(@mergefootprints,$merge);
}
else
{

}


# Calculate WT average (central) Cscore and Fscore of footprints identified in KO

foreach (@mergefootprints)

{
	chomp;
	my @sites=split(/\t/,$_);
	my $footprintslength=$sites[3];
	my $start=$sites[1];
	my $end=$sites[2];
	
	if( $footprintslength >= $lowrange && $footprintslength <= $uprange )
	{	

		## average score of Cscore on footprinting region (Central/Core Cscore)
		my $totalscore=0;
		my $averagescore=0;
		if ($start > $end)
		{
			for (my $i=$start+1; $i<=16299;$i++)
			{
				$totalscore+=$Cscore{$i};
			}
			for (my $j=1;$j<=$end;$j++)
			{
				$totalscore+=$Cscore{$j};
			}
			$averagescore=$totalscore/$footprintslength;
		}
		else
		{
			for (my $k=$start+1;$k<=$end;$k++)
			{
				$totalscore+=$Cscore{$k};

			}
			$averagescore=$totalscore/$footprintslength;
		}


	###	footprinting flanking score	&& start 0 or end 16299

		my $leftscore=0;
		my $rightscore=0;
		my $flankingscore=0;
		my $averageleft=0;
		my $averageright=0;
		
		for(my $m=$start;$m>($start-$flanking);$m--)
		{
			if($m <=0 )
			{ 
				$leftscore+=$Cscore{$m+16299};
			}
			else	    
			{
				$leftscore+=$Cscore{$m};
			}
		}

		for (my $n=$end+1;$n<($end+$flanking+1);$n++)
		{
			if ($n > 16299)
			{	
				$rightscore+=$Cscore{$n-16299};
			}
			else
			{	
				$rightscore+=$Cscore{$n};
			}

		}


		$averageleft=$leftscore/$flanking;
		$averageright=$rightscore/$flanking;

			
		$flankingscore=10**$averagescore/10**$averageleft+10**$averagescore/10**$averageright;
		print "$_\t";
		printf "%.6f\t%.6f\t%.6f\t%.6f\n",$averagescore,$flankingscore,$averageleft,$averageright; # central Cscore, Fsocre, left Cscore, right Cscore
		
	}
}

close(WT);
close(KO);


