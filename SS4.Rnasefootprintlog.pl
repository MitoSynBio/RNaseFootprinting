#!/usr/bin/perl -w

#using strict;

# Identified RNase Footprinting sites

die "Usage: $0 output (sort -k2,2n) of Rnasecutoff, Cscore of each sites (Cscore bed), lowrange, uprange {5, 50}nt, length of flanking(>=1)\n" if($#ARGV!=4);

my $Rnasecutoff=shift;
my $inputfile=shift;

my $lowrange=shift;
my $uprange=shift;
my $flanking=shift;

open(RNASE,$Rnasecutoff) || die "$!";
open(WIGGLEBED, $inputfile)||die "$!";



my @data=<WIGGLEBED>;
my %Cscore=(); ## Cscore hash

foreach (@data)

{
	chomp;
	my @temp=split(/\t/,$_);	
	## Hash: Cscore for every site
	$Cscore{$temp[2]}=$temp[3];
}




my @footprints=<RNASE>;
my @mergefootprints=();
my $mergestart=-1; ## mark 
my $mergeend=-1; ## mark 
my $mergelength=0;

### Merge chrM start and end possible footprints region

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



## link circular chrM 

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
#	print "$merge\n";
	push(@mergefootprints,$merge);
}
else
{

}


## filtering size of footprinting sites 
## calculate average Cscore and flanking score

foreach(@mergefootprints)

{
	chomp;
	my @sites=split(/\t/,$_);
	my $footprintslength=$sites[3];
	my $start=$sites[1]; # footprinting start site
	my $end=$sites[2]; # footprinting end site

	if( $footprintslength >= $lowrange && $footprintslength <= $uprange )	
	{	

		## average score of Cscore on footprinting region
		my $totalscore=0;
		my $averagescore=0;
		if ($start > $end)
		{
			for(my $i=$start+1; $i<=16299;$i++)
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

		if ($averageleft > $averagescore && $averageright > $averagescore)
		{			
				$flankingscore=10**$averagescore/10**$averageleft+10**$averagescore/10**$averageright;
				print "$_\t";
				printf "%.6f\t%.6f\t%.6f\t%.6f\n",$averagescore,$flankingscore,$averageleft,$averageright; # central Cscore, Fsocre, left Cscore, right Cscore
		}
	}
}

close(WIGGLEBED); 
close(RNASE);




