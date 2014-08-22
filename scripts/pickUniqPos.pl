#!/share/bin/perl
use Bio::Seq;
use List::Util qw(sum);

die "perl $0 <sam>\n" if @ARGV<1;
open in,$ARGV[0];
my %pe;
while(<in>)
{
	chomp;
	my @f=split/\t/,$_,12;
	## read number 1 or 2
	my ($rnum)=$f[1]=~/(\d)$/;

	## XT:A:* 
	#my ($xt)=$f[11]=~/XT:A:(.)/;
	my ($xt)=$f[4];

	my $strand="+";
	## revcomp
	if($f[1]=~/r/ &&  $f[9]!~/\*/)
        {
                my $seq=Bio::Seq->new(-seq=>$f[9]);
                $f[9]=$seq->revcom->seq;
		$strand="-";
        }

	## parse CIGAR
	if($xt > 0)
        {
                # CIGAR
                my (@cigar_m)=$f[5]=~/(\d+)M/g;
                my (@cigar_d)=$f[5]=~/(\d+)D/g;
                my (@cigar_s)=$f[5]=~/(\d+)S/g;
                my (@cigar_i)=$f[5]=~/(\d+)I/g;
                my $aln_ln=sum(@cigar_m,@cigar_d);
		
		print $f[2],"\t",$f[3]-1,"\t",$f[3]-1+$aln_ln,"\t$f[0]/$rnum\t",$f[9],"\t",$strand,"\n";
	}
}
close in;

