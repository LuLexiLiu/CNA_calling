#! /usr/bin/perl -w
use strict;
my %goodzone;
my %goodzonenum;
my ($input_bam, $dynabins, $output, $binflag)=@ARGV;
my $i=0;
my $total_read=0;
my $total_readnum=0;
open In1,"< $dynabins" or die $!;

while(<In1>)
{
	my @array=split;
	if(exists($goodzonenum{$array[0]})){
	    $goodzonenum{$array[0]}++
	    }
    else{
        $goodzonenum{$array[0]}=1
        };
	${$goodzone{$array[0]}}[$goodzonenum{$array[0]}]=[$array[1],$array[3],0,$array[5],$array[2]];
	unless($array[0]=~/X/ ||$array[0]=~/Y/ || $array[0]=~/M/){$total_readnum+=$array[5];}
}
close In1;

open Input,"< $input_bam" or die $!;
open Error,"> $output.dyna_$binflag.error" or die $!;
open Output,"> $output.dyna_$binflag.bin" or die $!;

while(<Input>)
{
	$_=~/^@/ and next;
	$_=~/XS:i/ and next;
	$_=~/chrM/ and next;
	my @array=split;
	$array[1] & 4 and next;
    !exists($goodzonenum{$array[2]}) and next;
	my $sequence=$_;
	my $up=$goodzonenum{$array[2]};
	my $down=1;
	my $flag=0;
	my $mid=int(($up-$down)/2);
	while($up-$down>=0)
	{
		if($goodzone{$array[2]}->[$mid][1]<$array[3])
		{
			$down=$mid+1;
			$mid=int(($up-$down)/2)+$down;
		}
		elsif($goodzone{$array[2]}[$mid][0]>$array[3])
		{
			$up=$mid-1;
			$mid=int(($up-$down)/2)+$down;
		}
		elsif($goodzone{$array[2]}[$mid][0]<=$array[3] && $goodzone{$array[2]}[$mid][1] >=$array[3])
		{
			$goodzone{$array[2]}[$mid][2]++;
			$flag=1;
			last;
		}
	}
	if($flag==1 && !($array[2]=~/M/) ){$total_read++}
	if($flag==0){print Error $sequence;}
}

foreach (1..22,"X","Y")
{
	my $chr="chr$_";
	foreach my$num(1..$goodzonenum{$chr})
	{
		my $ratio=($goodzone{$chr}->[$num][2]/$goodzone{$chr}->[$num][3])/($total_read/$total_readnum);
		print Output "$chr\t$goodzone{$chr}->[$num][0]\t$goodzone{$chr}->[$num][4]\t$goodzone{$chr}->[$num][2]\t$ratio\n"
	}
}
