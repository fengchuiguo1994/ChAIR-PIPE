use strict;
use warnings;

#copy the CB tag value to read ID
#add a tag for cell barcode from readID
my @bc;
my @b;
my $zb;
my @z;

my %dashnum = ('1' => 'A', '2'=>'C', '3'=>'G', '4'=>'T', '5'=>'N');

while(<>){
	if($_ =~ /^@/){
		print $_;
		next;
	}
        chomp $_;
	my @bits = split(/\t/, $_);
        @bc = grep { /^CB:Z/ } @bits;
        @b = split(/:/, $bc[0]);
        #$zb=substr($b[2], 0, -2); #use if dash not needed
        #@z = split(/-/, $b[2]);
        #$zb = $z[0].'-'.$dashnum{$z[1]};
        #push @bits, sprintf("ZB:Z:%s", $zb);
  
	$bits[0] = $b[2].':'.$bits[0];
	print join("\t", @bits), "\n";
}

