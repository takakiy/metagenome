#!/usr/bin/perl -w
# GC_skew_content

if ($ARGV[0] =~/\-h/){
  print "ARGV;  -if seqfiles -o outfile [OPTION]\n";
  print "OPTION; 
           -p program GCcont\(gc\) GCprofile\(gcp\) skew\(skw\) tm\(tm\)
           -w window width  (only gcp skw t)
           -sk skip length\n";
           
  exit;
}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

##open (OLOG,">out_logs");

unless (exists $option->{-if}){
  print "Sequence ?\n";
  chomp($in=<STDIN>);
  $option->{-if}=[$in];
}

if (exists $option->{-p}) {
	if ($option->{-p}->[0] !~ /gc|gcp|skw|tm/) {
		print "program Select gc gcp skw t?\n";
		chomp($pro=<STDIN>);
		$option->{-p}->[0]= $pro;
	}	
} else {
  print "program GCcont\(gc\) GCprofile\(gcp\) skew\(skw\) tm\(tm\)?\n";
  chomp($pro=<STDIN>);
  $option->{-p}->[0]= $pro;
}

if ($option->{-p}->[0] =~ /(gcp|skw|tm)/){
  unless (exists $option->{-w}) {
		print "Window?\n";
		chomp($win=<STDIN>);
  		$option->{-w}->[0]= $win;
  }
  unless (exists $option->{-sk}) {
		print "Skip?\n";
		chomp ($skip=<STDIN>);
		$option->{-sk}->[0]= $skip;
  }  
}

$option->{-o}->[0]= exists $option->{-o} ? $option->{-o}->[0] : "out".$option->{-p}->[0];


&control_calculation($option);



##===============================================================================
##  MAIN
##===============================================================================
sub control_calculation {
	my ($option)= @_;
	my ($seqinfo,$one,$seqone,$seqlen,@total_val);
	my ($baseinfo,$gc_cont,$ag_cont,$win,$skip,$i,$start,$cent,$end);
	my ($gc_skew,$at_skew,$tm);
	my (($data_num,$av,$min,$max,$sd));
	
	($seqinfo)= &readseq_fasta($option->{-if});
	open (OUT, ">$option->{-o}->[0]");
	print OUT "## ID\tLen\tAMB\tGC\tAG\tINDEX\n" if ($option->{-p}->[0] eq 'gc');
	print OUT "## ID\tPROG\tST\tED\tVAL\n" if ($option->{-p}->[0] =~ /gcp|skw|tm/);
	
	foreach $one (@{$seqinfo->{tax}}) {
		$seqone= $seqinfo->{seq}->{$one};
		$seqlen= $seqinfo->{len}->{$one};
		$seqindex= $seqinfo->{index}->{$one};
		@total_val= ();
		if ($option->{-p}->[0] eq "gc") {
##print OLOG ">$one\n$seqone\n";
			($baseinfo)= &cal_base_count($seqone);
			$gc_cont= $baseinfo->{gc_rate};
			$ag_cont= $baseinfo->{ag_rate};
			print OUT "$one\t$seqlen\t$baseinfo->{n}\t$baseinfo->{gc_rate}\t$baseinfo->{ag_rate}\t$seqindex\n";
#			print "SEQ= $one $seqlen AMB= $baseinfo->{n} GC= $baseinfo->{gc_rate} AG= $baseinfo->{ag_rate}\n";
		} elsif ($option->{-p}->[0] =~ /gcp|skw|t/) {
			($win,$skip)= ($option->{-w}->[0],$option->{-sk}->[0]);
			$win= $win==0 ? $seqlen : $win;
			$skip= $skip==0 ? $seqlen : $skip;
	#		print "SEQ= $one SLEN= $seqlen PROG= $option->{-p}->[0] WIN= $win SKIP= $skip\n";
			for ($i=0; $i <= $seqlen-$skip; $i += $skip){
				$subseq= substr($seqone,$i,$win);
				($start,$cent,$end)= ($i+1,$i+int($win/2),$i+$win);
				($baseinfo)= &cal_base_count($subseq);
				if ($option->{-p}->[0] eq "gcp") {
					$gc_cont= $baseinfo->{gc_rate};
					$ag_cont= $baseinfo->{ag_rate};
					print OUT "gcp\t$start\t$end\t$gc_cont\t$one\t$skip\t$ag_cont\n";
					push @total_val,$gc_cont;
				} elsif ($option->{-p}->[0] eq "skw") {
					($gc_skew,$at_skew)= &cal_gc_at_skew($baseinfo);
					print OUT "skw\t$start\t$end\t$gc_skew\t$one\t$skip\n";
					push @total_val,$gc_skew;
				} elsif ($option->{-p}->[0] eq "tm") { 
#print "$one\n";
					($tm)= &cal_tm_NN($subseq);
					if ($win > 0) {
						print OUT "tm\t$start\t$end\t$tm\t$one\n";
						push @total_val,$tm;
					} else {
						print OUT "$one\t$seqlen\ttm\n";
					}						
				}
			}
		}
		if (@total_val > 1){
   		($data_num,$av,$min,$max,$sd)= &cal_av_max_min_sd([@total_val]);
  		 	print "## $one\t$data_num\t$av\t$min\t$max\t$sd\n";
		}
	}
	print "=============  FINISH  ============\n";
}



###  CALCULATION  ========================================================
sub cal_base_count {
	my ($seq)= @_;
	my ($baseinfo);
	my ($gc_num,$gc_cont,$ag_num,$ag_cont);
	
	$baseinfo= {};
	$_=$seq;
	$baseinfo->{a}= tr/Aa//;
	$baseinfo->{t}= tr/Tt//;
	$baseinfo->{g}= tr/Gg//;
	$baseinfo->{c}= tr/Cc//;
	$baseinfo->{n}= tr/[ATGCatgc]//c;
	$baseinfo->{base}= $baseinfo->{a}+$baseinfo->{t}+$baseinfo->{g}+$baseinfo->{c};
	
   $gc_num= $baseinfo->{g}+$baseinfo->{c};
   $gc_cont= $gc_num != 0 ? $gc_num/$baseinfo->{base} : 0;
   $baseinfo->{gc_rate}= sprintf ("%.4f",$gc_cont);
   $ag_num= $baseinfo->{g}+$baseinfo->{a};
   $ag_cont= $ag_num != 0 ? $ag_num/$baseinfo->{base} : 0;
   $baseinfo->{ag_rate}= sprintf ("%.4f",$ag_cont);	

	return ($baseinfo);
}

sub cal_gc_at_skew {
	my ($baseinfo)=@_;
	my ($gc_num,$gc_diff,$gc_skew);
	my ($at_num,$at_diff,$at_skew);
 
   $gc_num= $baseinfo->{g}+$baseinfo->{c};
   $gc_diff= $baseinfo->{g}-$baseinfo->{c};
   $gc_skew= $gc_num != 0 ? $gc_diff/$gc_num : 0;
   $gc_skew= sprintf ("%.4f",$gc_skew);
   $at_num= $baseinfo->{a}+$baseinfo->{t};
   $at_diff= $baseinfo->{a}-$baseinfo->{t};
   $at_skew= $at_num != 0 ? $at_diff/$at_num : 0;
   $at_skew= sprintf ("%.4f",$at_skew);
   return ($gc_skew,$at_skew);
}
sub cal_tm_wal {
	## Wallace formula 17-25 base
	my ($seq)= @_;
	my ($baseinfo,$gc_num,$at_num,$tm);
	
	($baseinfo)= &cal_base_count($seq);
	$gc_num= $baseinfo->{g}+$baseinfo->{c};
	$at_num= $baseinfo->{a}+$baseinfo->{t};
	$tm= sprintf ("%.2f",($gc_num*4+$at_num*2));

	return ($tm);
}
sub cal_tm_standard {
	## Current Protocols in Molecular Biology
	## 50mM KCL
	## 10mM Tris
	## 1.5mM MgCl2
	## Tm= 81.5+16.6*log[S]+0.41*(%GC)-(500/n)
	## K+ cont + Tris cont*0.65
	## [S]= 0.05+0.67*0.01= 0.0567 (M)
	## Tm= 60.8+0.41*(%GC)-(500/n)
	
	my ($seq)= @_;
	my ($baseinfo,$ag_cont,$gc_contP,$seqlen,$tm);
	($baseinfo)= &cal_base_count($seq);
	$gc_contP= $baseinfo->{gc_rate}*100;
	$seqlen= length $seq;
	$tm= sprintf ("%.2f",(60.8+0.41*$gc_contP-(500/$seqlen)));

	return ($tm);
}
sub cal_tm_NN {
	## Nearest Neignbor Method
	## Tm= (1000*sum_DeltaH)/(-10.8+sumDeltaS+R*ln(Oligo conc M/4))-273.15+16.6log[Na+]
	## Na conc 0.05 M / Oligo conc 0.5x10-6 M
	my ($seq)= @_;
	my ($leng,%twoH,%twoS);
	my ($sumH,$sumS,$y,$two,$na,$ct);
	
	$seq=~ tr/atgc/ATGC/;
	$leng= length $seq;
	%twoH= (
	'AA',-9.1,
	'TT',-9.1,
	'AT',-8.6,
	'TA',-6,
	'CA',-5.8,
	'TG',-5.8,
	'GT',-6.5,
	'AC',-6.5,
	'CT',-7.8,
	'AG',-7.8,
	'GA',-5.6,
	'TC',-5.6,
	'CG',-11.9,
	'GC',-11.1,
	'GG',-11,
	'CC',-11);
	
	%twoS= (
	'AA',-24,
	'TT',-24,
	'AT',-23.9,
	'TA',-16.9,
	'CA',-12.9,
	'TG',-12.9,
	'GT',-17.3,
	'AC',-17.3,
	'CT',-20.8,
	'AG',-20.8,
	'GA',-13.5,
	'TC',-13.5,
	'CG',-27.8,
	'GC',-26.7,
	'GG',-26.6,
	'CC',-26.6);

	($sumH,$sumS)= (0,0);
	for ($y=0; $y < $leng-1; $y++) {
		$two=substr $seq,$y,2;
		$sumH += $twoH{$two};
		$sumS += $twoS{$two};
	}

	$na= -1.30103;
	$ct= -15.8949521;

	$up=1000*$sumH/($sumS+(1.987*$ct)-10.8);

	$tm= sprintf ("%.2f",($up-273.15+(16.6*$na)));
	return ($tm);

}    

###  STATICS  ========================================================
sub cal_av_max_min_sd {
  my ($vals)=@_;
  my ($sum,@val_sort);
  my ($data_num,$av,$min,$max,$sd);
  
  @val_sort= sort {$a <=> $b}@$vals;
  $sum= 0;
  $sum += $_ foreach (@val_sort);
  $data_num= @val_sort;
##		print "$data_num  @val_sort\n";
  $av= sprintf ("%.4f",$sum/$data_num);
  $min= $val_sort[0];
  $max= $val_sort[-1];
  
  my ($diff,$diff_2);
  my $sumdiff=0;
  foreach $val (@val_sort){
     $diff= $val-$av;  $diff_2= $diff**2;
     $sumdiff += $diff_2;
  }
  $sd=sprintf ("%.4f",sqrt($sumdiff/$data_num));
  return ($data_num,$av,$min,$max,$sd);
}

##==============================Sequence Format ==============================
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








