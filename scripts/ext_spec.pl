#! /usr/bin/perl -w

if ($ARGV[0] =~ /\-h/){
  print "-l List(Id<\t>tempSeq<\t>Dir<\t>Start<\t>End(<\t>INFO))
   -d TempSeqFile -o (m/s)\n";
  print "OPTION\n";
  print "-e (numeric OR upstream_downstream), expand from both term\n"; 
  print "-ts (numeric OR upstream_downstream), extract both term\n";
  print "-tm (numeric OR upstream_downstream), merged both term\n";
  print "-et (numeric), output the flanking and term regions\n";
  print "-f (numeric OR upstream_downstream), output both flanking region\n";
  print "-info T/(F) adding infomation.....Only normal extraction\n";
  exit;
}
  
foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};
$opts= join " ",keys %{$option};

$option->{-info}->[0]= exists $option->{-info} ? $option->{-info}->[0] : "F";
$list = $option->{-l}->[0];

 ########################  OPTIONS  #####################
unless (exists $option->{'-o'}) {
   print "outfile format Single (s) or Multifasta (m) ?";
   $form=<STDIN>;  chomp $form;  $option->{'-o'}->[0]=$form;
}
$option->{'-o'}->[0]= "m" if ($opts=~ /(^|\s)-et(\s|$)/);

if ($option->{'-o'}->[0] eq "s") {
  system ('mkdir TEST') unless (-e "TEST") ;
} elsif ($option->{'-o'}->[0] eq "m") {
  if ($opts=~ /\be\b/){
    $output_m= "extra_exp.nuc";
  } elsif ($opts=~ /\bt(s|m)\b/){
    $output_m= "extra_ter.nuc";
  } elsif ($opts=~ /\bf\b/){
    $output_m= "extra_flank.nuc";
  } elsif ($opts=~ /\bet\b/){
    $output_m= "extra_Term.nuc";
  } else {
    $output_m= "extra.nuc";
  }
  open (OUT, ">$output_m");
}

##================================  MAIN LOOP ================================

open (LI, "$list");   @list=<LI>; 
open (ERR,">ext_spec.logs");

#%hasnt=();


($seqinfo)= &readseq_fasta($option->{-d});

$query_no= 0;
$extra_no= 0;

foreach (@list){
  chomp;
  next if (/^(#|\s*$)/);
  $query_no++;

  @items=split/\t/;
  $des= "";
  if ($items[2]=~/\+|\-/){
     ($id,$nt,$dir,$st,$ed)=@items[0..4];
     $des= $items[5] if (exists $option->{-info});
  }elsif ($items[3]=~/\+|\-/) {
     ($id,$nt,$dir,$st,$ed)=@items[0,1,3,4,5];
     $des= $items[6] if (exists $option->{-info});
  }
  $dir = $dir =~ /^(\-|Minus|C)$/ ? '-' : '+';

  if (exists $seqinfo->{seq}->{$nt}){
     $seq=$seqinfo->{seq}->{$nt};
		$extra_no++;
  }else{

		print "ERROR: $id Do not sequence: $nt;  STOP\n";
		print ERR "ERROR: $id Do not sequence: $nt;  STOP\n";
#		exit;
 #    ($seq)=&make_tempseq($nt);
  }
  next if ($seq eq "");

   ##================================  EXTRACT ================================
    print "$nt $ed\n";
    $stt=$st > $ed ? $ed : $st;
    $edd=$st > $ed ? $st : $ed;
    if ($opts =~ /\be\b/){
      $val= $option->{'-e'}->[0];
      ($valup,$valdw)= &check_vals($val); 
      $seqnts= &extract_core_flank($seq,$dir,$stt,$edd,$valup,$valdw);
      $seqnt= $$seqnts[1];
      $seqnt= "$$seqnts[0]\n$seqnt" if ($$seqnts[0] ne "");
      $seqnt= "$seqnt\n$$seqnts[2]\n" if ($$seqnts[2] ne "");
      $result= ">${id}_exp\n$seqnt\n";
      $outfile= "${id}_exp.nuc";
   } elsif ($opts =~ /\bf\b/){
      ($up_seq,$dw_seq)=("","");
      $val= $option->{'-f'}->[0];
      ($valup,$valdw)= &check_vals($val); 
      $seqnts= &extract_core_flank($seq,$dir,$stt,$edd,$valup,$valdw);
      $up_seq= ">${id}_flkUP\n$$seqnts[0]\n" if ($$seqnts[0] ne "");
      $dw_seq= ">${id}_flkDW\n$$seqnts[2]\n" if ($$seqnts[2] ne ""); 
      $result= $up_seq.$dw_seq;
      $outfile= "${id}_flk.nuc";
   } elsif ($opts =~ /\bt(m|s)\b/){
      ($up_seq,$dw_seq)=("","");  
      $val= exists $option->{'-tm'} ? $option->{'-tm'}->[0] : $option->{'-ts'}->[0] ;
      ($valup,$valdw)= &check_vals($val); 
      $seqnts= &extract_term($seq,$dir,$stt,$edd,$valup,$valdw);
      $seqnt= join "\n",@$seqnts;    $seqnt=~ s/^\n//;
      if ($opts =~ /\btm\b/){ 
          $result= ">${id}_ter\n$seqnt\n";
      }elsif ($opts =~ /\bts\b/){
          $up_seq= ">${id}_terUP\n$$seqnts[0]\n" if ($$seqnts[0] ne "");
          $dw_seq= ">${id}_terDW\n$$seqnts[1]\n" if ($$seqnts[1] ne ""); 
          $result= $up_seq.$dw_seq;
      }
      $outfile= "${id}_ter.nuc";
   } elsif ($opts =~ /\bet\b/){
       $val= $option->{'-et'}->[0];
      ($valup,$valdw)= &check_vals($val); 
       $seqnt_exp= &extract_core_flank($seq,$dir,$stt,$edd,$valup,$valdw);
       $seqnt_ter= &extract_term($seq,$dir,$stt,$edd,$valup,$valdw);
       $dw_ter_lu= &make_lineSeq($$seqnt_ter[0]);
       $dw_ter_ld= &make_lineSeq($$seqnt_ter[1]);
       $dw_ter_ldr= &reverse_seq($dw_ter_ld);
       $dw_exp_lu= &make_lineSeq($$seqnt_exp[0]);
       $dw_exp_ld= &make_lineSeq($$seqnt_exp[2]);
       $dw_exp_ldr= &reverse_seq($dw_exp_ld);
       $result= "$id\t$dw_exp_lu\t$dw_ter_lu\t$dw_ter_ld\t$dw_exp_ld\t";
       $result= $result."$dw_exp_ldr\t$dw_ter_ldr\n";
   }  else {

      $seqnts= &extract_core_flank($seq,$dir,$stt,$edd,0,0);
    #  $$seqnts[1]=~ tr/\-//d;
    #  $$seqnts[1]=~ tr/U/T/;
        
      $result= $option->{-info}->[0] ne "F" ? ">$id  $des\n".$$seqnts[1]."\n" :
               ">$id\n".$$seqnts[1]."\n";
      $outfile= "$id.nuc";
   }
   
   if ($option->{-o}->[0] eq "s") {
      open (OUT,">./TEST/$outfile");
      print OUT "$result";
   } else {
      print OUT "$result";
   }

}


print "QUERY: $extra_no\n";

##================================  SUB ROUTING ================================
sub make_tempseq {
    my ($nt)=@_;
    my ($ntf,$refseq,$pardir,$seq);
    my $ch=0;
    if (exists $option->{'-d'}) {
       $ntf= $option->{'-d'}->[0];
       if (-f $ntf){ 
         $refseq=&create_sequence($ntf);    %hasnt=@$refseq;
         $ch=1 if (exists $hasnt{$nt});
       }   
    }else {
       if (exists $ENV{'GENOME_DATABASE'}){
          $pardir=$ENV{'GENOME_DATABASE'};  $ntf=$pardir.'/'.$nt;
          if (-f $ntf){
             $refseq=&create_sequence($ntf);    %hasnt=@$refseq;  $ch=1;
          }
       }
    }
    if ($ch != 1){
       print "Sequence file of $nt?\n";  $ntf=<STDIN>;  chomp $ntf;
       if (-f $ntf){
          $refseq=&create_sequence($ntf);    %hasnt=@$refseq;
          $ch=1 if (exists $hasnt{$nt});
       }
    }
    if ($ch != 1){
       print "$id SEQUENCE FILE  ##  $nt ; DONT EXIST \n";
       $seq="";
    }else {
       $seq=$hasnt{$nt};
    }    
    return ($seq,$ch);
}

sub check_vals {
   my ($val)=@_;
   if ($val =~ /_/){
      $valup= (split/_/,$val)[0]; 
      $valdw= (split/_/,$val)[1];
   }else{
      $valup=$val;   $valdw=$val
   }
   return ($valup,$valdw);
}

sub extract_core_flank {
  my ($seq,$dir,$stt,$edd,$valup,$valdw)=@_;
  my ($seq_up,$seq_dw,$seq_core);

  $seq_up="";  $seq_dw="";  $seq_core="";
  if ($valup != 0){
		$seq_up= $dir eq '+' ? substr $seq,$stt-$valup-1,$valup : substr $seq, $edd,$valup;
		$seq_up= &reverse_seq($seq_up) if ($dir eq '-');
		$seq_up= &limit_seqline($seq_up);
  } 
  if ($valdw != 0){
    $seq_dw= $dir eq '+' ? substr $seq,$edd,$valdw : substr $seq, $stt-$valdw-1,$valdw;
    $seq_dw= &reverse_seq($seq_dw) if ($dir eq '-');
    $seq_dw= &limit_seqline($seq_dw);
  }
  $seq_core= substr $seq, $stt-1,$edd-$stt+1;  
  $seq_core= &reverse_seq($seq_core) if ($dir eq '-');
  $seq_core=&limit_seqline($seq_core);

  return [$seq_up,$seq_core,$seq_dw];
}

sub extract_term{
	my ($seq,$dir,$stt,$edd,$valup,$valdw)=@_;
	my ($up_term,$dw_term) = ("","");
	if ($valup != 0){
		if ($edd-$stt+1 >= $valup) {
			$up_term= $dir eq '+' ? substr $seq,$stt-1,$valup : substr $seq,$edd-$valup,$valup;
			$up_term= &reverse_seq($up_term) if ($dir eq '-');
			$up_term= &limit_seqline($up_term);
		} else {
			$up_term= substr $seq,$stt-1,$edd-$stt+1;
			$up_term= &reverse_seq($up_term) if ($dir eq '-');
			$up_term= &limit_seqline($up_term);
		}
  }
  if ($valdw != 0){  
		if ($edd-$stt+1 >= $valup) {
			$dw_term= $dir eq '+' ? substr $seq,$edd-$valup,$valup : substr $seq,$stt-1,$valup;  
			$dw_term= &reverse_seq($dw_term) if ($dir eq '-');
			$dw_term= &limit_seqline($dw_term);
		} else {
			$dw_term= substr $seq,$stt-1,$edd-$stt+1;  
			$dw_term= &reverse_seq($dw_term) if ($dir eq '-');
			$dw_term= &limit_seqline($dw_term);
		}
  }
  return [$up_term,$dw_term];
}


 ########################  Create sequence  ##################### 
sub create_sequence {
  my ($nt) =@_;
  my (@nts,$co,%has,@seqs,$seq,@seqlist,$ref);
  @seqlist=();
  open (INNT, "$nt");
  @nts=<INNT>;
  foreach(@nts){
    chomp;
    if (/^>(\S+)/){
      $co=$1;
    }else{
      push @{$has{$co}}, $_;
    }
  }
  foreach (sort keys %has){
    @seqs=@{$has{$_}};
    $seq=join "",@seqs;
    push @seqlist, $_;
    push @seqlist, $seq;
  }
  $ref=[@seqlist];
  return $ref;
  close (INNT);
}
 ########################  Line to one #####################
sub make_lineSeq {
  my ($seq)=@_;
  my @seqs=split/\n/,$seq;
  my $seq_l=join "",@seqs;
  return $seq_l;
}
 ########################  reverse #####################
sub reverse_seq {
  my ($seq)=@_;
  my ($sqa);
  local $_;
  $_=$seq;
  tr/ATGCUatgc/TACGAtacg/;
  $sqa = reverse ($_);
}
 ########################  limit_seqline #####################
sub limit_seqline {
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



##==================================================================
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
