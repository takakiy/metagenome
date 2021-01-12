#! /usr/bin/perl -w


if ($ARGV[0]=~/\-h/){
   print "OPTION -s NucSeqFile -gc (norm)/mit (-outn output name) --->  NT-->PEP\n";
   print "OPTION -l orf_List file -ltype (norm)/draft/genome -d NucSeqFile -gc (norm)/mit -outf m(multifasta)/s (-fix 0 1 2) -outn outputfile name
  (norm: Tempseq CDS_ID Dir Start End (codon_start: 1 2 3))
  (draft: multi-gbrown format in draft genome)
  (genome: single-gbrown format in complete genome)
";


   print "OPTION -h help\n";
   exit;
}

%{$geneticCode->{norm}} = qw(
    TTT F TTC F TTA L TTG L CTT L CTC L CTA L CTG L 
    ATT I ATC I ATA I ATG M GTT V GTC V GTA V GTG V 
    TCT S TCC S TCA S TCG S CCT P CCC P CCA P CCG P 
    ACT T ACC T ACA T ACG T GCT A GCC A GCA A GCG A 
    TAT Y TAC Y TAA * TAG * CAT H CAC H CAA Q CAG Q 
    AAT N AAC N AAA K AAG K GAT D GAC D GAA E GAG E 
    TGT C TGC C TGA * TGG W CGT R CGC R CGA R CGG R 
    AGT S AGC S AGA R AGG R GGT G GGC G GGA G GGG G 
);

%{$geneticCode->{mit}} = qw(
	TTT F TCT S TAT Y TGT C TTC F TCC S TAC Y TGC C
	TTA L TCA S TAA * TGA W TTG L TCG S TAG * TGG W
	CTT L CCT P CAT H CGT R CTC L CCC P CAC H CGC R
	CTA L CCA P CAA Q CGA R CTG L CCG P CAG Q CGG R
	ATT I ACT T AAT N AGT S ATC I ACC T AAC N AGC S
	ATA M ACA T AAA K AGA S ATG M ACG T AAG K AGG S
	GTT V GCT A GAT D GGT G GTC V GCC A GAC D GGC G
	GTA V GCA A GAA E GGA G GTG V GCG A GAG E GGG G
);


##===============================================================================

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};


#unless (exists $option->{'-outf'}) {
#   print "outfile format Single (s) or Multifasta (m) ?";
#   $form=<STDIN>;
#   chomp $form;
#   $option->{'-outf'}->[0]= $form;
#}
unlink "erro_log.txt" if ( -f "erro_log.txt");
open (ERR,">>erro_log.txt");

$genecode= exists $option->{-gc} ? $option->{-gc}->[0] : "norm";
$outform= exists $option->{'-outf'} ? $option->{'-outf'}->[0] : "m";

if (exists $option->{-l}){
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

    $optfix= exists $option->{-fix} ? $option->{-fix} : [1,2,3];
    $list_type= exists $option->{-ltype} ? $option->{-ltype}->[0] : "norm";
    
    $genes= &make_genelist_nuc($list_type,$option->{-l}->[0],$seqinfo,$optfix);

    ### OUTPUT FILE
    
    if ($outform eq "s") {
        system ('mkdir PEPs') unless (-e "PEPs");
    } elsif ($outform eq 'm'){
        
        if (exists $option->{'-outn'}) {
            $outn="$option->{'-outn'}->[0]";
  
            open (OUT, ">${outn}_cds.pep");
            open (OUTN, ">${outn}_cds.nuc");
            
            print "### OUTPUT ${outn}_cds.pep\n";

        } else {
            $outhead= $option->{-l}->[0];
        #    $outhead=~ s/^.*\///;
        #    $outhead=~ s/.\w+$//;

            open (OUT, ">${outhead}_cds.pep");
            open (OUTN, ">${outhead}_cds.nuc");
            
            print "### OUTPUT ${outhead}_cds.pep\n";

        }
        
    }

    
} elsif (exists $option->{-s}){
    
    if (-f $option->{-s}->[0]){
        ($seqinfo)= &readseq_fasta($option->{-s});
        @geness=();
        push @geness, [$_,"X",$seqinfo->{seq}->{$_}] foreach (@{$seqinfo->{tax}});
        $genes=[@geness];
    }else{
        print "Do not exists SeqFile\n";  exit;
    }
    
    $outhead= exists $option->{'-outn'} ? $option->{'-outn'}->[0] : (split/\//,$option->{-s}->[0])[-1];
    
    open (OUT, ">${outhead}.pep");
    open (OUTN, ">${outhead}.code.nuc");
    
    ###  FRAME TEST
    ($genes)= &control_frame_test($genes,$genecode);
    
    print "OUTPUT: ${outhead}.pep  ${outhead}.code.nuc\n";

}

##===============================================================================

$q_num=@$genes;
print "$q_num\n";

##===============================================================================
$|=1;  $c=0;

## OUTPUT
$outform= exists $option->{'-outf'} ? $option->{'-outf'}->[0] : 'm';

foreach (@$genes){
   $c++;
   ($name,$des,$orfseq)=@$_;
   ($cdsid)=(split/\s+/,$name)[0];
    
   if ($orfseq eq "NO_TEMPSEQ"){
      print ERR "$cdsid NO_TEMP_SEQ\n";
      next;
   }

    $dif=(length $orfseq)%3;
    
    if ($dif != 0) {
       print ERR "$cdsid\t\tERRO_LENG\n";
       next;
    }
    
    $start_codon=substr $orfseq,0,3;
    print ERR "$cdsid\t$start_codon\tSTART_ERRO\n" unless($start_codon=~/ATG|TTG|GTG/);
    $end_codon=substr $orfseq,-3,3;
    print ERR "$cdsid\t$end_codon\tEND_ERRO\n" unless ($end_codon=~/TAA|TAG|TGA/);
    
    
    ##  TO Amino acid
    $orf_aa= &nuc_to_amino($orfseq,$genecode);
    $orf_aaSeq= &multiline_seq($orf_aa);
    $stops= &count_stop($orf_aaSeq);
    if ($stops > 1) {
       print ERR "$cdsid $stops STOP_IN_SEQ_ERRO\n";
     ## next;
    }
   
    if ($outform eq 's'){
        $outfile='./PEPs/'.$cdsid.'.pep';
        open (OUT,">$outfile");
        
    } elsif ($outform eq 'm') {

        substr($orf_aaSeq,-1,1)= "";
        print OUT "\>$name  $des\n$orf_aaSeq\n";
        $orf_nuc=&multiline_seq($orfseq);
        print OUTN "\>$name  $des\n$orf_nuc\n";
        
   }

}

print "ERR: erro_log.txt\n";



##===============================================================================
##===============================================================================
sub control_frame_test {
   my ($genes,$genecode)= @_;
   my (@new_genes,$name,$des,$orfseq,$info);
   my ($bestinfo);
   
   unlink "conv_best.logs.txt" if (-e "conv_best.logs.txt");
   unlink "conv_best.fas" if (-e "conv_best.fas");
   open (COV,">out_conv_best.logs.txt");
   open (COS,">out_conv_best.fas");
   print COV "GENE\tLEN\tSTOP\tB_FRAM\tB_STOP\tB_LEN\n";
   
   @new_genes= ();
   foreach (@$genes) {
     ($name,$des,$orfseq)=@$_;
     if ($orfseq eq "NO_TEMPSEQ"){
        push @new_genes,$_;
       next;
     }
     $info= {};
     ($info)= &frame_test($orfseq,$genecode);

     print COV "$name\t$info->{ori}->{len}\t$info->{info}->{1}->{stop}\t";
     print COV "$info->{best}->{fram}\t$info->{best}->{stop}\t";
     print COV "$info->{best}->{seqlen}\n";
     print COS ">$name  $info->{best}->{fram}  $des\n$info->{best}->{seq}\n";
   
     push @new_genes,[$name," Frame $info->{best}->{fram}  $des",$info->{best}->{seq}];
   }
    
    print "OUT_BEST_LOG: out_conv_best.logs.txt\n";
    print "OUT_BEST_FASTA: out_conv_best.fas\n";

    
   return ([@new_genes]);

}
sub frame_test {
   my ($orfseq,$genecode)= @_;
   my ($info,$len,$cc,$tagseq,$seq,$mod,$subseq,$endcodon,$orf_aa,$stops,@seqord);
   
   $info= {};
   $len= length $orfseq;
   $info->{ori}->{len}= $len;
   $info->{ori}->{seq}= $orfseq;
   foreach $cc (1,2,3,-1,-2,-3) {
      $tagseq= $cc > 0 ? $orfseq : &complement_seq($orfseq);
      $seq= substr ($tagseq,abs($cc)-1,$len-abs($cc)-1);
      $mod= ($len-abs($cc)-1)%3;
      $subseq= substr($seq,0,($len-abs($cc)-1)-$mod);
      
      $endcodon= substr($subseq,-3,3);
      if ($endcodon =~ /TAA|TAG|TGA/) {
         substr($subseq,-3,3)= "";
      }
       $orf_aa= &nuc_to_amino($subseq,$genecode);
       $stops= &count_stop($orf_aa);
# print "$cc \n$orf_aa \n $stops\n";
      $info->{info}->{$cc}->{seq}= $subseq;
      $info->{info}->{$cc}->{stop}= $stops;
   }
   @seqord= sort {$info->{info}->{$a}->{stop}<=>
         $info->{info}->{$b}->{stop}}(keys %{$info->{info}});
   $info->{best}->{fram}= $seqord[0];
   $info->{best}->{seq}= $info->{info}->{$seqord[0]}->{seq};
   $info->{best}->{seqlen}= length $info->{info}->{$seqord[0]}->{seq};
   $info->{best}->{stop}= $info->{info}->{$seqord[0]}->{stop};
   
    
   return ($info);
}
##===============================================================================

sub make_genelist_nuc {
	my($list_type,$list,$seqinfo,$optfix)=@_;
	my(@lists,$cdsid,$dir,$st,$ed,$des,$note,$type,$fix);
	my(@genes);
	@genes=();

	open(LI,$list); 
	@lists=<LI>;
	foreach (@lists){
        chomp;
        next if(/^(\s*$|\#)/);
        @item=split/\t/;
		if ($list_type eq "norm"){
            ($cdsid,$tempseq,$dir,$st,$ed)=@item[1,0,2,3,4];
            $des= @item > 5 ? $item[6] : "";
            $type="CDS";
            $fix= 1;
            $des= "" unless(defined $des);
            if (@item == 6) {
                $frame=$item[5];
                $st= $dir eq '+' ? $st+$frame-1 : $st+($ed-$st+1)%3;
                $ed= $dir eq '+' ? $ed-($ed-$st+1)%3 : $ed-$frame+1;
            }
		} elsif ($list_type eq "draft"){
			($tempseq,$type,$cdsid,$dir,$st,$ed,$des,$fix)=@item[0,2,3,5,6,7,8,11];
		} elsif ($list_type eq "genome"){
			$tempseq= $seqinfo->{tax}->[0];
			($type,$cdsid,$dir,$st,$ed,$des,$fix)=@item[0,1,3..6,9];
		}

        next if ($type ne "CDS");

        $fixpattern= join '|',@$optfix;
		next if ($fix !~ /^($fixpattern)$/);

        #      print "$cdsid ";
        
		if (exists $seqinfo->{seq}->{$tempseq}){
            $stt= $st < $ed ? $st : $ed;
            $edd= $st < $ed ? $ed : $st;
			$subseq= substr($seqinfo->{seq}->{$tempseq},$stt-1,$edd-$stt+1);
			$subseq= &complement_seq($subseq) if ($dir eq '-');
		} else {
			$subseq= "NO_TEMPSEQ";
		}
		push @genes,[$cdsid,$des,$subseq,$frame];
	}

	return [@genes];
}


##===============================================================================
sub count_stop {
  my ($seq)=@_;
  my ($stops);
  $_= $seq;
  $stops=s/\*//g;
  $stops= $stops > 0 ? $stops : 0;
  return $stops;
}

sub nuc_to_amino {
    my ($orfseq,$genecode)=@_;
    my ($cdsid)= @_ > 1 ? $_[1] : "no";
    my ($len_orf,@aa_code,$cut,$code,$amin);
	my ($orf_aa);
	$len_orf = length $orfseq;
	@aa_code=();
	my $icuno= 0;

	my %genetic_code= %{$geneticCode->{$genecode}};

    for ($cut=0; $cut<$len_orf; $cut+= 3){
        $code = substr $orfseq,$cut,3;
        if (exists $genetic_code{$code}){
            $amin=$genetic_code{$code};
        } else {
            $amin="X";
            $icuno++;
        }
        #      $amin="M" if ($cut == 0);

        push @aa_code,$amin;
    }
  
    if ($icuno > 0) {
        print ERR "ICU $cdsid $icuno\n";
    }
    $orf_aa=join "",@aa_code;
    return $orf_aa;
}
##===============================================================================

sub request_file {
    print "Where is tempfile?\n";
    $i_file=<STDIN>;
    if (-f $i_file) {
        return ($i_file,1);
    }else{
        return ("No Tempfile\n",0);
    }
}

##===============================================================================

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
