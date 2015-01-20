use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use Utility;
use ParseSam;


#### UNMAPPED###
#HWUSI-EAS300R:7:1:10:750#0      77      *       0       0       *       *       0       0       TCATAACATGAATAAATATGTTCGCAAGTTTGTAAACAACTGCATAACTTGTAGGACTTCCAAGTCGC
#CCTCTG      Xbb\ba`^b`Oaa\KbabbbbaOVb^X_PRbb_\aabaaX__K`VbX__\aababM_aV_^a\``ab_[FL`\b
#HWUSI-EAS300R:7:1:10:750#0      141     *       0       0       *       *       0       0       ATTTTAAGGGTGTGATACAGGTGTACAAGTTTCGTGAAGCCATCAATCTGAACGATGACATATTCTTT
#TAG a]`^\XG`^ODVJYGWGXT^[DWX_O_TFN`_``G[O_Z]SWO_XOOZGVUHS_OW\ZXPLLSZZXO^UXX
#HWUSI-EAS300R:7:1:10:463#0      83      2R      6819009 60      74M     =       6818842 -241    TAGTGGTGTATGTACATACATACATACTGCCTTACTTTTTTGGACCGTCCGCACTCATTAATTGCAAA
#CCGTTT      `YWPU_^V]_\W^_aWXX_]WSW_YP^a`YFX_YPaaaa`_a^\a_aa`ab_``_VY]\baab_U^a`a`bbba      XT:A:U  NM:i:1  SM:i:37 AM:i:37 X0:i:1  X1:i:0  XM:i:1  XO:i:0  XG:i:0  
#MD:Z:44A29


my $fastq1="";
my $fastq2="";
my $sam1="";
my $sam2="";
my $output="";
my $help=0;
my $debug=0;

GetOptions(
    "fq1=s"             =>\$fastq1,
    "fq2=s"             =>\$fastq2,
    "sam1=s"	        =>\$sam1,
    "sam2=s"            =>\$sam2,
    "output=s"          =>\$output,
    "debug"             =>\$debug,
    "help"	        =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);
pod2usage(-verbose=>2) if $help;

my $sr1=SamParser->new($sam1);
my $sr2=SamParser->new($sam2);
my $fqr1=Utility::get_fastq_reader($fastq1);
my $fqr2=Utility::get_fastq_reader($fastq2);

open my $ofh,">",$output or die "Could not open output file";

# first print the header
my @h1=@{$sr1->{header}};
my @h2=@{$sr2->{header}};
for my $i (0..(@h1-1))
{
    my $hl1=$h1[$i];
    my $hl2=$h2[$i];
    die "heaer in the two sam files are not identical" unless $hl1 eq $hl2;
    print $ofh $hl1."\n"; 
}

my $s1="";
my $s2="";
while(1)
{
    # first read some fastqfiles
    my $fq1=$fqr1->();
    my $fq2=$fqr2->();
    die "fastq files do not have equal length" if(not $fq1 and $fq2 or not $fq2 and $fq1);
    last unless $fq1;
    die "fastq files not 'in sync'; " unless $fq1->{readid} eq $fq2->{readid};
    
    my $readid=$fq1->{readid};
    $s1=$sr1->nextsam() unless $s1;
    $s2=$sr2->nextsam() unless $s2;
    
    if($debug){
        print "S1: $s1->{readid}\tS2: $s2->{readid}\tFQ: $readid\n";
    }
    
    my($a1,$a2);
    if($s1 and $s2 and $s1->{readid} eq $readid and $s2->{readid} eq $readid)
    {
        $a1=$s1;
        $a2=$s2;
        $s1=""; # both need to be reread in the next round
        $s2="";
    }
    elsif($s1 and $s1->{readid} eq $readid)
    {
        $a1=$s1;
        $s1=""; # only reread the first in the next round
        $a2=Utility::samsoniteFastq($fq2);
    }
    elsif($s2 and $s2->{readid} eq $readid)
    {
        $a1=Utility::samsoniteFastq($fq1);
        $a2=$s2;
        $s2=""; # only reread the second sam in the next round
    }
    else
    {
        # both are existing but none of them matches the fastqreadid
        $a1=Utility::samsoniteFastq($fq1);
        $a2=Utility::samsoniteFastq($fq2);
        # none of the two need to be reread
    }
    
    Utility::crosslinkPair($a1,$a2);
    print $ofh Utility::formatSam($a1)."\n";
    print $ofh Utility::formatSam($a2)."\n";
    my $bla=0;
}
if($s1)
{
    die "sam is not empty in the end"
}

close $ofh;
exit;

{
    package Utility;
    use strict;
    use warnings;
    sub samsoniteFastq
    {
        my $fastq=shift;
        
    #HWUSI-EAS300R:7:1:10:750#0      77      *       0       0       *       *       0       0       TCATAACATGAATAAATATGTTCGCAAGTTTGTAAACAACTGCATAACTTGTAGGACTTCCAAGTCGC
    #CCTCTG      Xbb\ba`^b`Oaa\KbabbbbaOVb^X_PRbb_\aabaaX__K`VbX__\aababM_aV_^a\``ab_[FL`\b
        # readid, flag, chr, pos, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
        my $sam=$fastq;
        $sam->{flag}=0;
        $sam->{flag}|=0x0004; # unmapped
        $sam->{chr}="*";
        $sam->{pos}=0;
        $sam->{mq}=0;
        $sam->{cigar}="*";
        $sam->{chrmate}="*";
        $sam->{posmate}=0;
        $sam->{distance}=0;
        $sam->{appendix}="";
        return $sam
    }

    sub crosslinkPair
    {

        my $r1=shift;
        my $r2=shift;
        
    
        $r1->{flag}|=0x0001; # paired in sequence
        $r2->{flag}|=0x0001; # paired in sequence
        $r1->{flag}|=0x0040; # first in a pair
        $r2->{flag}|=0x0080; # second in a pair
        
        crosslink($r1,$r2);
        crosslink($r2,$r1);
    }

    sub crosslink
    {
        my $a=shift; # the thing to update
        my $temp=shift;  # the mate
        
        
    
        if ($temp->{flag} & 0x0004)
        {
            $a->{flag}|=0x0008; # the mate is unmapped; set this flag even if the read itself is unmapped
            return undef;
        }
        return undef if $a->{flag} & 0x004;
        
        # FROM here on both are mapped
        # strand
        $a->{flag}|=0x0020 if $temp->{flag} & 0x0010;
        
        my $pp=1;
        
        # readid, flag, chr, pos, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
        $a->{distance}=0;
        if($a->{chr} eq $temp->{chr})
        {
            $a->{chrmate}="=";
            $a->{distance}=$temp->{pos}-$a->{pos};
        }
        else
        {
            $a->{chrmate}=$temp->{chr};
            $pp=0;
        }
        
        $a->{posmate}=$temp->{pos};
        
        $pp=0 if (($a->{flag} & 0x0010) == ($temp->{flag} & 0x0010));
        $a->{flag} |=0x0002 if $pp;
    }


    sub formatSam
    {
        my $p=shift;
        # readid, flag, chr, pos, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
        #HWUSI-EAS300R:7:1:1:422#0/1     0       X       7618779 26      1S73M   *       0       0       NTTACTTTAGTTTCACTGGTAGCATGGAGTTGTAAGTGAAAACTATTTATACTCCTTTTCAATATA
    #TGATTCNN      DPZZZY[[YY[[ZYY[Z[WVZ[[[[YZ[V[[[ZZ[[YYZ[[ZZ[Z[[[[[[[Z[Y[Y[[[[[[[[[[ZZBBBBB      AS:i:73 XS:i:0  XF:i:3  XE:i:1  XN:i:0
    my $str="$p->{readid}\t$p->{flag}\t$p->{chr}\t$p->{pos}\t$p->{mq}\t$p->{cigar}\t$p->{chrmate}\t$p->{posmate}\t$p->{distance}\t$p->{seq}\t$p->{qual}";
    $str.="\t$p->{appendix}" if $p->{appendix};
    return $str;
    }
    
    sub get_fastq_reader
    {
        my $inputfile=shift;
        open my $fh,"<",$inputfile or die "Could not open input file";
        #@HWUSI-EAS300R:7:1:0:697#0/1
        #NTCAAAAATGATTTTTCAATGTGTTTAAGGTNANNNNNNNNNCNNTNNNNNNNNNNNNNNNNNNNNNNNNNNNN
        #+HWUSI-EAS300R:7:1:0:697#0/1
        #DOOSV[[ZYQQUZXRXNXZTBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

        return sub
        {
            my $h1=<$fh>;
            my $seq=<$fh>;
            my $h2=<$fh>;
            my $qual=<$fh>;
            return undef unless $h1;
            chomp $h1; chomp $seq; chomp $h2; chomp $qual;
            $h1=~s/^@//;
            $h1=~s{/[12]$}{};
            
            return
            {
                readid  =>  $h1,
                seq     =>  $seq,
                qual    =>  $qual
            };
        }
    }

}



{
    package SamParser;
    use strict;
    use warnings;

    sub new
    {
        my $cls=shift;
        my $input=shift;

        open my $ifh,"<",$input or die "Could not open input file 1";
        my $self=bless
        {
            input   =>  $input,
            fh     =>  $ifh,
            buffer   =>  [],
            header=>[],
        }, __PACKAGE__;
        
        $self->_spoolFwd();
        return $self;
    }
    
    sub _spoolFwd
    {
        my $self=shift;

        
        while(my $l=$self->_nextline)
        {
            unless($l=~m/^@/)
            {
                $self->_bufferline($l);
                last;
            }
            chomp $l;
            push @{$self->{header}},$l;
        }
        
    }
    
    
    
    sub nextsam
    {
        my $self=shift;
        my $l=$self->_nextline;
        return undef unless $l;
        chomp $l;
        
        my $s=$self->_parseSam($l);
        my $readid=$s->{readid};
        
        #treat the fucking chimeras!!
        while(my $nl=$self->_nextline)
        {
            chomp $nl;
            my $ns=$self->_parseSam($nl);
            if($ns->{readid} ne $readid)
            {
                $self->_bufferline($nl);
                last;
            }
        }

		# This warning has been commented out as it cannot occur in mcclintock but still produces errors?
		#warn "if of read is supposed to end with /1 or /2 : $readid" unless $readid=~m{/[12]$};
        $readid=~s{/[12]$}{};
        $s->{readid}=$readid;
        return $s;
    }
    
    sub _nextline
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub _bufferline
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }
    
    
    sub _parseSam
    {
        my $self=shift;
        my $line=shift;
        my @a=split /\t/,$line;
        
        
        # readid, flag, chr, pos, mq, cigar, chrmate, posmate, distance, seq, qual, appendix
        my $entry=
        {
            readid=>$a[0],
            flag=>$a[1],
            chr=>$a[2],
            pos=>$a[3],
            mq=>$a[4],
            cigar=>$a[5],
            chrmate=>$a[6],
            posmate=>$a[7],
            distance=>$a[8],
            seq=>$a[9],
            qual=>$a[10],
            appendix=>$a[11]
        };
        
        return $entry;
    }
}

#GetOptions(
#    "fq1=s"             =>\$fastq1,
#    "fq2=s"             =>\$fastq2,
#    "sam1=s"	        =>\$sam1,
#    "sam2=s"            =>\$sam2,
#    "output=s"          =>\$output,
#    "debug"             =>\$debug,
#    "help"	        =>\$help

=head1 NAME

perl samro.pl - Merges two sam-files and crosslinks the paired end information (e.g.: when using bwa sw)

=head1 SYNOPSIS

 perl samro.pl --fq1 reads1.fastq --fq2 reads2.fastq --sam1 maped-reads1.sam --sam2 maped-reads2.sam --output merged.sam

=head1 OPTIONS

=over 4

=item B<--fq1>

The first read of a pair in fastq. Mandatory

=item B<--fq2>

The second read of a pair in fastq. Mandatory

=item B<--sam1>

the mapping results of the first read in sam. Mandatory

=item B<--sam2>

the mapping results of the second read in sam. Mandatory

=item B<--output>

The output file. [sam]  Mandatory.

=item B<--help>

Display the help

=back

=head1 DETAILS

=cut

