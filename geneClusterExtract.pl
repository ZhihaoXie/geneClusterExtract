#!/usr/bin/perl -w
use File::Basename;
use Bio::SeqIO;

my $usage = "USAGE:\n\tperl $0 query_seq db_seq out_prefix extend_length\n" .
"\n\tdb_seq format is:" .
"\n\t>start_gene direction\n\tfasta_seq" .
"\n\t>end_gene direction\n\tfasta_seq" .
"\n\n\tthe start_gene seq must in the beginning of db_seq, and direction must be 1 or -1" .
"\n\textend_length default 0\n";

if (@ARGV < 3) {
    print STDERR "$usage";
    exit;
}

my $query = $ARGV[0];
my $db_seq = $ARGV[1];
my $out_prefix = $ARGV[2];
my $extend_len = $ARGV[3];
unless ($extend_len) {
    $extend_len = 0;
}

# set begin and end flag(db_seq's id)
my $begin_flag;
my $end_flag;
my $begin_strand;
my $end_strand;

# get the length of query or db_seq
my %length_query;
my %length_db;
my %query_seqs;
my $q_obj = Bio::SeqIO->new(-file => "$query", -format=>'fasta');
while (my $seq = $q_obj->next_seq()) {
    my $id = $seq->display_id;
    $length_query{$id} = $seq->length;
    $query_seqs{$id} = $seq;
}

my $db_obj = Bio::SeqIO->new(-file => "$db_seq", -format=>'fasta');
while (my $seq = $db_obj->next_seq()) {
    my $id = $seq->display_id;
    my $desc = $seq->desc();
    if (!defined $desc || $desc eq "") {
        $desc = 1;
    } elsif ($desc ne "1" && $desc ne "-1") {
        $desc = 1;
    }
    if (!defined $begin_flag) {    # set begin and end flag(db_seq's id)
        $begin_flag = $id;
        $begin_strand = $desc;
    } else {
        $end_flag = $id;
        $end_strand = $desc;
    }
    $length_db{$id} = $seq->length;
}

# run blastn
my $bsn_out = $out_prefix . '.db.bsn';
unless (-e "$query.nsq") {
    system("formatdb -p F -i $query")==0 or die "Failed running formatdb!";
}
system("blastall -p blastn -d $query -i $db_seq -o $bsn_out -e 1e-5 -m 8 -a 10 -F F")==0 or die "Failed running blastn!";

# make the best hit for blast
my %hits;
my %tmp_score;
my %in_scaff;
open BSN, $bsn_out;
while (<BSN>) {
    chomp;
    my @line = split /\t/, $_;
    if ($line[3] >= 200) {    # set filt length is 200
        if (!defined $hits{$line[0]}) {
            $hits{$line[0]} = $_;
            $tmp_score{$line[0]} = $line[11];
            $in_scaff{$line[0]} = $line[1];
        } elsif ($line[11] > $tmp_score{$line[0]}) {
            $hits{$line[0]} = $_;
            $tmp_score{$line[0]} = $line[11];
            $in_scaff{$line[0]} = $line[1];
        }
    }
}
close BSN;

if ($in_scaff{$begin_flag} ne $in_scaff{$end_flag}) {

    open POSI, "> $out_prefix.partial.position";
    open OUT_FASTA, "> $out_prefix.partial.fasta";
    my @m = split /\t/, $hits{$begin_flag};
    my @n = split /\t/, $hits{$end_flag};
    my $str_begin;
    my $str_end;
    if ($begin_strand == 1) {
        my ($start, $end);
        if ($m[8] <= $m[9]) {
            $start = $m[8]-$m[6]-$extend_len;
            $end = $length_query{$m[1]};
            if ($start <= 1) {
                $start = 1;
            }
            $str_begin = $query_seqs{$m[1]}->subseq($start, $end);
            print POSI "$m[1]\t$start\t$end\tbegin\n";
        } else {
            $start = $m[8]+$m[6]+$extend_len;
            if ($start >= $length_query{$m[1]}) {
                $start = $length_query{$m[1]};
            }
            $end = 1;
            $str_begin = $query_seqs{$m[1]}->subseq($end, $start);
            my $x = Bio::Seq->new(-seq=>$str_begin, -id=>'1');
            $x = $x->revcom;
            $str_begin = $x->seq;
            print POSI "$m[1]\t$start\t$end\tbegin\n";
        }
    } elsif ($begin_strand == -1) {
        my ($start, $end);
        if ($m[8] <= $m[9]) {
            $start = $m[9]+($length_db{$m[0]}-$m[7])+$extend_len;
            $end = 1;
            if ($start >= $length_query{$m[1]}) {
                $start = $length_query{$m[1]};
            }
            $str_begin = $query_seqs{$m[1]}->subseq($end, $start);
            my $x = Bio::Seq->new(-seq=>$str_begin, -id=>'1');
            $x = $x->revcom;
            $str_begin = $x->seq;
            print POSI "$m[1]\t$start\t$end\tbegin\n";
        } else {
            $start = $m[9]-($length_db{$m[0]}-$m[7])-$extend_len;
            $end = $length_query{$m[1]};
            if ($start <= 1) {
                $start = 1;
            }
            $str_begin = $query_seqs{$m[1]}->subseq($start, $end);
            print POSI "$m[1]\t$start\t$end\tbegin\n";
        }
    }
    if ($end_strand == 1) {
        my ($start, $end);
        if ($n[8] <= $n[9]) {
            $start = 1;
            $end = $n[9]+($length_db{$n[0]}-$n[7])+$extend_len;
            if ($end >= $length_query{$n[1]}) {
                $end = $length_query{$n[1]};
            }
            $str_end = $query_seqs{$n[1]}->subseq($start, $end);
            print POSI "$n[1]\t$start\t$end\tend\n";
        } else {
            $start = $length_query{$n[1]};
            $end = $n[9]-($length_db{$n[0]}-$n[7])-$extend_len;
            if ($end <= 1) {
                $end = 1;
            }
            $str_end = $query_seqs{$n[1]}->subseq($end, $start);
            my $x = Bio::Seq->new(-seq=>$str_end, -id=>'1');
            $x = $x->revcom;
            $str_end = $x->seq;
            print POSI "$n[1]\t$start\t$end\tend\n";
        }
    } elsif ($end_strand == -1) {
        my ($start, $end);
        if ($n[8] <= $n[9]) {
            $start = $length_query{$n[1]};
            $end = $n[8]-$n[6]-$extend_len;
            if ($end <= 1) {
                $end = 1;
            }
            $str_end = $query_seqs{$n[1]}->subseq($end, $start);
            my $x = Bio::Seq->new(-seq=>$str_end, -id=>'1');
            $x = $x->revcom;
            $str_end = $x->seq;
            print POSI "$n[1]\t$start\t$end\tend\n";
        } else {
            $start = 1;
            $end = $n[8]+$n[6]+$extend_len;
            if ($end >= $length_query{$n[1]}) {
                $end = $length_query{$n[1]};
            }
            $str_end = $query_seqs{$n[1]}->subseq($start, $end);
            print POSI "$n[1]\t$start\t$end\tend\n";
        }
    }
    print OUT_FASTA ">$m[0] $m[1]\n$str_begin\n>$n[0] $n[1]\n$str_end\n";
    close POSI;
    close OUT_FASTA;

} else {
    my ($start, $end);
    my @m = split /\t/, $hits{$begin_flag};
    my @n = split /\t/, $hits{$end_flag};
    if ($m[8] <= $m[9] && $n[8] <= $n[9] && $m[9] <= $n[8]) {
        $start = $m[8]-$m[6]-$extend_len;
        $end = $n[9]+($length_db{$n[0]}-$n[7])+$extend_len;
    } elsif ($m[8] <= $m[9] && $n[8] > $n[9] && $m[9] <= $n[9]) {
        $start = $m[8]-$m[6]-$extend_len;
        $end = $n[8]+$n[6]+$extend_len;
    } elsif ($m[8] > $m[9] && $n[8] <= $n[9] && $m[8] <= $n[8]) {
        $start = $m[9]-($length_db{$m[0]}-$m[7])-$extend_len;
        $end = $n[9]+($length_db{$n[0]}-$n[7])+$extend_len;
    } elsif ($m[8] > $m[9] && $n[8] > $n[9] && $m[8] <= $n[9]) {
        $start = $m[9]-($length_db{$m[0]}-$m[7])-$extend_len;
        $end = $n[8]+$n[6]+$extend_len;
    } elsif ($m[8] <= $m[9] && $n[8] <= $n[9] && $m[8] >= $n[9]) {
        $start = $m[9]+($length_db{$m[0]}-$m[7])+$extend_len;
        $end = $n[8]-$n[6]-$extend_len;
    } elsif ($m[8] <= $m[9] && $n[8] > $n[9] && $m[8] >= $n[8]) {
        $start = $m[9]+($length_db{$m[0]}-$m[7])+$extend_len;
        $end = $n[9]-($length_db{$n[0]}-$n[7])-$extend_len;
    } elsif ($m[8] > $m[9] && $n[8] <= $n[9] && $m[9] >= $n[9]) {
        $start = $m[8]+$m[6]+$extend_len;
        $end = $n[8]-$n[6]-$extend_len;
    } elsif ($m[8] > $m[9] && $n[8] > $n[9] && $m[9] >= $n[8]) {
        $start = $m[8]+$m[6]+$extend_len;
        $end = $n[9]-($length_db{$n[0]}-$n[7])-$extend_len;
    }

    my $position_out = "$out_prefix.cluster.position";
    my $out_fasta = "$out_prefix.cluster.fasta";
    open POSI, "> $position_out";
    open OUT, "> $out_fasta";
    if ($start <= $end) {
        if ($start <=1) {
            $start = 1;
        }
        if ($end >= $length_query{$m[1]}) {
            $end = $length_query{$m[1]};
        }
        my $str = $query_seqs{$m[1]}->subseq($start, $end);
        print POSI "$m[1]\t$start\t$end\tcluster\n";
        print OUT ">cluster $out_prefix\n$str\n";
    } else {
        if ($start >= $length_query{$m[1]}) {
            $start = $length_query{$m[1]};
        }
        if ($end <= 1) {
            $end = 1;
        }
        my $str = $query_seqs{$m[1]}->subseq($end, $start);
        my $x = Bio::Seq->new(-seq=>$str,-id=>'1');
        $x = $x->revcom;
        $str = $x->seq;
        print POSI "$m[1]\t$start\t$end\tcluster\n";
        print OUT ">cluster $out_prefix\n$str\n";
    }
    close POSI;
    close OUT;
}
