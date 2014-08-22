#!/share/bin/perl
use Bio::Seq;
use List::Util qw(sum);

die "perl $0 <input.insertion.refined.bp> <fragment size>\n" if @ARGV<1;

my $title=$ARGV[0];
$title =~ s/insertion.refined.bp/sorted.bam/;
my $frag=$ARGV[1];

open (input, "<$ARGV[0]") or die "Can't open $ARGV[0] since $!\n";
print "Chr\tStart\tEnd\tTransposonName\tTransposonDirection\tClass\tVariantSupport\tFrequency\tJunction1\tJunction1Support\tJunction2\tJunction2Support\t5\'_Support\t3\'_Support\n";
while (my $line=<input>) {
    chomp($line);
    my @a=split(/\t/, $line);

    my @b=split(/\(/, $a[6]);
    my @c=split(/\(/, $a[7]);
    my $terminal="5\'";
    my $reverse=0;
    my $positive=$b[0];
    my $negative=$c[0];
    if ($b[0] > $c[0]) {
	$terminal="3\'";
	$reverse=1;
	my $swap=$b[0];
	$b[0]=$c[0];
	$c[0]=$swap;
    }
    my $lower=$b[0]-$frag;
    my $upper=$c[0]+$frag;
    system("samtools view -bu $title $a[0]\:$lower\-$upper > temp.bam");
    system("samtools view -Xf 0x2 temp.bam > temp.sam");
    
    open in,"temp.sam";
    my %ps=();
    my %me=();
    my $ref_sup=0;
    my $soft_clip=0;
    while(<in>)
    {
	chomp;
	my @f=split/\t/,$_,12;
	## read number 1 or 2
	my ($rnum)=$f[1]=~/(\d)$/;
	
	## XT:A:* 
	my ($xt)=$f[11]=~/XT:A:(.)/;
	
	## Coordinate
	if ($f[1]=~/r/)
	{
	    my (@cigar_m)=$f[5]=~/(\d+)M/g;
	    my (@cigar_d)=$f[5]=~/(\d+)D/g;
	    my (@cigar_s)=$f[5]=~/(\d+)S/g;
	    my (@cigar_i)=$f[5]=~/(\d+)I/g;
	    my $aln_ln=sum(@cigar_m,@cigar_d);
	    $me{$f[0]}=$f[3]+$aln_ln-1;
	}
	else
	{
	    $ps{$f[0]}=$f[3];
	}	    
	
#	${$pe{$f[0]}}[$rnum-1]=[$xt,$coor];
    }
    close in;

    
    foreach my $id (keys %ps)
    {
	if (defined $me{$id})
	{	    
	    if (((($ps{$id}+5)<=$positive)&&($me{$id}>$negative)) || ((($me{$id}-5)>=$negative)&&($ps{$id}<$positive))) {
#	    if (((($ps{$id}+5)<=$b[0])&&($me{$id}>$c[0])) || ((($me{$id}-5)>=$c[0])&&($ps{$id}<$b[0]))) {
		$ref_sup++;
#		print "$id\n";
	    }
	}
    }

    my $variant=$a[8]+$a[9]+$a[10]+$a[11];
    #my $ratio=sprintf("%.4f", $variant/($variant+$ref_sup));
    if (($a[0] =~ /^\d{1,2}$/) || ($a[0] eq "X") || ($a[0] eq "Y")) {$a[0]="chr$a[0]";}
    if ($reverse == 0) {
	print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$variant\tratio\t$b[0]\t$a[8]\t$c[0]\t$a[9]\t$a[10]\t$a[11]\n";
    }
    elsif ($reverse == 1) {
	print "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$variant\tratio\t$b[0]\t$a[9]\t$c[0]\t$a[8]\t$a[10]\t$a[11]\n";
    }
    system("rm temp.sam temp.bam");
}

