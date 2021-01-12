#! /usr/bin/perl -w


if ($ARGV[0]=~/\-h/){

   print "OPTION -l gene_List file -ltype (draft)/genome -d NucSeqFile -outf m(multifasta)/s (-fix 0 1 2) -outn outputfile name
  (norm: Tempseq CDS_ID Dir Start End (codon_start: 1 2 3))
  (draft: multi-gbrown format in draft genome)
  (genome: single-gbrown format in complete genome)
";

   print "OPTION -o OUTPUT FORMAT Single(s) or Multi(m) FASTA\n";
   print "OPTION -h help\n";
   exit;
}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

###  GENE LIST

$optfix= exists $option->{-fix} ? $option->{-fix} : [1,2,3];

$list_type= exists $option->{-ltype} ? $option->{-ltype}->[0] : "draft";
($geneinfo)= &create_geneinfo($option->{-l}->[0],$list_type);

$out_mode= exists $option->{'-outf'} ? $option->{'-outf'}->[0] : 'm';

###  DTABASE
unless (exists $option->{'-d'}) {
     print "Where is TempSeq File ?";
     $tempfile=<STDIN>;  chomp $tempfile;
     $option->{'-d'}->[0]= $tempfile;
}

unless (-f $option->{'-d'}->[0]){
     print "Do not TempSeq File\n";
     exit;
}

($seqinfo)= &readseq_fasta($option->{-d});

### OUTPUT FILE

if ($out_mode eq "s") {
        system ('mkdir NUCs') unless (-e "NUCs");
        print "OUTPUT: ./NUCs/\n";
    
} elsif ($out_mode eq 'm'){
        
    if (exists $option->{'-outn'}) {
        $outn="$option->{'-outn'}->[0]";
        open (OUTN, ">${outn}_gene.fna");
        print "### OUTPUT ${outn}_gene.fna\n";

    } else {
        $outhead= (split/\//,$option->{-l}->[0])[-1];
        
 #       print "BB $outhead\n";
        if ($outhead=~ /\.\w+$/) { $outhead=~ s/\.\w+$//; } else {}
        
        
 #       print "AA $outhead\n";
        open (OUTN, ">${outhead}_gene.fna");
        print "### OUTPUT ${outhead}_gene.fna\n";
    }
}


###
$count= 0;

foreach $contig (@{$geneinfo->{contigs}}) {
    
    $ctgseq= $seqinfo->{seq}->{$contig};
    
    foreach $gene (@{$geneinfo->{genes}->{$contig}}) {
        ($type,$geneID,$dir,$st,$ed,$prod,$fix)= @{$geneinfo->{genes}->{$gene}}[2,3,5,6,7,8,11];
        
        if ($type =~ /CDS|RNA/) {}
        
        $g_nlen= $ed-$st+1;
        $g_seq= substr($ctgseq,$st-1,$g_nlen);
        $g_seq= &complement_seq($g_seq) if ($dir eq '-');

        $fixpattern= join '|',@$optfix;
        next if ($fix !~ /^($fixpattern)$/);

        print OUTN ">$geneID $prod\n".&multiline_seq($g_seq)."\n";
        
        $count++;
    }
    
}

    print "COUNT: $count\n";



##=== === === === === === === === === === === === === === === === ===

sub create_geneinfo {
    my($list,$list_type)=@_;
    my ($geneinfo);

    open(LI,$list);
    @lists=<LI>;
    foreach (@lists){
        chomp;
        next if(/^(\s*$|\#)/);
        @item_all=split/\t/;
        @item = @item_all > 15 ? @item_all[0..15] : @item_all[0..13];
        
        if ($list_type eq "draft"){
            push @{$geneinfo->{contigs}},$item[0] unless (exists $geneinfo->{genes}->{$item[0]});
            push @{$geneinfo->{genes}->{$item[0]}}, $item[3];
            $geneinfo->{genes}->{$item[3]}= [@item];
        }
        
    }
    return ($geneinfo);

}



##===============================================================================
##===============================================================================
sub count_stop {
  my ($seq)=@_;
  my ($stops);
  $_= $seq;
  $stops=s/\*//g;
  $stops= $stops > 0 ? $stops : 0;
  return $stops;
}

sub readseq_fasta {
    my ($seqfile)= @_;
    my ($seqinfo,$seqname,$fileone,@lines,$index,$dupinfo,$onefinfo);
  
    $seqinfo={};
    $dupinfo= {};
    foreach $fileone (@$seqfile) {
        open (SEQ,"$fileone");
        @lines= <SEQ>;
        $onefinfo= {};
        foreach(@lines){
            chomp;
            next if (/^(\#|\s*$)/);
            if (/^>(\S+)/){
                $seqname=$1;
                $index= $_;
                $dupinfo->{$seqname}++;
                $dupinfo->{total}++;
                if ($dupinfo->{$seqname} == 1) {
                    push @{$onefinfo->{tax}},$seqname;
                    $index=~ s/^>(\S+)\s*//;
                    $seqinfo->{index}->{$seqname}= $index;
                }
            }else{
                push @{$seqinfo->{seqline}->{$seqname}},$_ if ($dupinfo->{$seqname} == 1);
            }
        }
        
        foreach $seqname (@{$onefinfo->{tax}}) {
            if (exists $seqinfo->{seqline}->{$seqname}) {
                $seqinfo->{taxnum}++;
                push @{$seqinfo->{tax}},$seqname;
            } else {
                delete $seqinfo->{index}->{$seqname};
                print "   NOSEQ $seqname\n";
            }
        }
        
        print "FILE $fileone SEQ_NUMBER  $seqinfo->{taxnum}\/$dupinfo->{total}\n";
    }
  
    foreach $seqname (sort keys %{$dupinfo}) {
        next if ($seqname=~ /total/);
        print "  DUPLI $seqname $dupinfo->{$seqname}\n" if ($dupinfo->{$seqname} > 1);
    }
    
    ($seqinfo)= &join_seqline($seqinfo);
  

    return ($seqinfo);
  
}
sub join_seqline {
    my ($seqinfo)= @_;
    my ($sequence);
  
    foreach $seqname (@{$seqinfo->{tax}}) {
        $sequence= join "",@{$seqinfo->{seqline}->{$seqname}};
        $sequence=~ s/\s+//g;
        $sequence=~ tr/a-z/A-Z/;
        $sequence=~ tr/U/T/;
        $seqinfo->{seq}->{$seqname}= $sequence;
        $seqinfo->{len}->{$seqname}= length $sequence;
        delete $seqinfo->{seqline}->{$seqname};
        # print "$seqname $seqinfo->{len}->{$seqname}\n";
    }
    return ($seqinfo);
}

sub multiline_seq {
    my $seq= shift @_;
    my $unit = defined $_[0] ? $_[0] : 60;
    my ($len,@subs,$subseq,$x,$line_seq);
    $len=length $seq;
    @subs=();
    for ($x= 0; $x < $len; $x += $unit){
        $subseq=substr $seq,$x,$unit;
        push @subs,$subseq;
    }
    $line_seq=join "\n", @subs;
    return $line_seq;
}

sub complement_seq {
  my ($seq)= @_;
  my ($revseq);

  $seq=~ tr/ATGCatgc\-NBDHVRYKMSW/TACGtacg\-NVHDBYRMKSW/;
  $revseq = reverse($seq);
  return $revseq;
}




