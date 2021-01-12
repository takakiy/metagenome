#! /usr/bin/perl -w


$EM_MODE= 0;   ## 0: ECO mode 1: FASTA MODE
$setmode= 0;


 if ((@ARGV == 0) || ($ARGV[0] =~ /\-h/)) {

  print "

	##  EXTRACT GENOME SEQ / EXTRACT INDEX SEQ /CONVERT SINGLE FASTA FROM MULTIFASTA

	## SELECT EXEC MODE
	-e [o(org_seq)/s(single fasta)/e(extract)/c(check)/
       d(split_fasta)/r(removeGap)/x(cut)/n(one_line for glimmer)/
       v(rename)]
	-d database

	## EXTRACT OTION
		0) EXTRACT SEQ IN ORG
		ex.  manu_fasta.pl -e o -d dbseq -s org1 org2 .....

		## AND_SEARCH OF TWO KEYWORD
		1) DIRECT SPECIFIED
			SPECIFIED WITH OPTION -w1 KEYWORD KEYWORD -w2 KEYWORD KEYWORD

			IN '-e es'
				-w1 organism
				-w2 String

		2) LIST
			SPECIFIED WITH OPTION -l listfile
			### listfile
			>KWD1
			KEYWORD KEYWORD KEYWORD  ***  group
			>KWD2
			KEYWORD
			KEYWORD
			....
		## output
		-o s(single)/g(multifasta for genome)/[m(multifasta for all)]
			** s|g ---> EACH FOLDER

		ex.  manu_fasta.pl -e e -d dbseq -w1 ww ww -w2 ww ww  -o g
		ex.  manu_fasta.pl -e e -d dbseq -l list  -o m
	
	## DIVIDED MULTIFASTA ( by -unit number of seq & by -lsub number of getting file)
		ex. manu_fasta.pl -e d -d dbseq [-u 400] (-lsub 0) -outdir outdir

	## CONVERT SINGLE FASTA 
		ex.  manu_fasta.pl -e s -d dbseq ---> EACH FOLDER 

	## CHECK SEQ (REMOVE LENGTH 0 and Duplicate SEQ
		ex.  manu_fasta.pl -e c -d dbseq ---> rev_manuseq.fas
		ex.  manu_fasta.pl -e c -d dbseq -lmlen 100 ## Filter length---> rev_manuseq.fas

	## REMOVE GAP[-.]
		ex.  manu_fasta.pl -e r -d dbseq ---> rev_infile.fas

	## CUTTING SEQ x
		ex.  manu_fasta.pl -e x -d dbseq -a 1 200 (-b 20 15 both term)---> rev_infile.fas

	## RENAME SEQ
		ex.  manu_fasta.pl -e v -d dbseq -l listfile[ORINAME(ID) NEW_NAME] ---> rev_infile.fas

	## SPLIT one SEQ in a unit　(basepair)
		ex.  manu_fasta.pl -e m -d dbseq -u number ---> rev_stpseq.fas

	## CONCATENATE SEQ WITH GAP
		ex.  manu_fasta.pl -e sc -d dbseq -gap number ---> rename_scseq

	## COMPLEMENT
		ex.  manu_fasta.pl -e a -d dbseq ---> rev_infile.fas
		
	## STAT
		ex.  manu_fasta.pl -e stat -d dbseq  ---> out_manu_fasta_stat.logs
		
	## SAMPLING
		ex.  manu_fasta.pl -e sample -d dbseq -numb 100 ---> sampling_seq.fas

	## FILTER-OUT LENGTH
		ex.  manu_fasta.pl -e len -d dbseq -lmin 100 (-lmax xx)---> out_filter.fas
            
    ## NON-INTERLEAVE SEQ
            ex.  manu_fasta.pl -e noint -d dbseq ---> out_nointer.fas

		\n";

print "SETTING EM_MODE (Index only = 0, SEQ= 1) $EM_MODE\n";
print "SETTING setmode $setmode\n";


      exit;
}

#==========================================================================

## $in= $ENV{HOME}.'/Desktop/kegg/ftp.genome.jp/pub/kegg/genes/fasta/genes.pep';

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

open (LOG, ">manu_fasta.logs");

$outmode= exists $option->{-o} ? $option->{-o}->[0] : "m";
$exemode= exists $option->{-e} ? $option->{-e}->[0] : "s";

#==========================================================================
#  CONTROL EXEC MODE
#==========================================================================
$EM_MODE = 1 if ($exemode=~ /^([sdcrxnovma]|len)/);
$EM_MODE = 0 if ($exemode=~ /stat/);

($seqinfo)= &readseq_fasta($option->{-d}) if ($EM_MODE == 1);


if ($exemode eq 'e') {

	($selinfo)= &read_list($option->{-l}->[0]) if (exists $option->{-l});
	foreach (sort keys %{$option}) {
		if (/^\-w/) {
			$selinfo->{keynum}++;
			push @{$selinfo->{info}->{$selinfo->{keynum}}},@{$option->{$_}};
		}
	}
	if ($selinfo->{keynum} > 0) {
		print "SEL $_ @{$selinfo->{info}->{$_}}\n" foreach (sort keys %{$selinfo->{info}});
	}

	if ($EM_MODE == 1) {
		($extrainfo)= &extract_sel_seqname($seqinfo,$selinfo,$exemode);
		&output_extract_seq($seqinfo,$extrainfo,$outmode);
	} elsif ($EM_MODE == 0) {
		&output_extract_seq_eco2($option->{-d},$selinfo,$exemode,$outmode);
	}

} elsif ($exemode eq 'o') {

	$selinfo->{reg}= exists $option->{-s} ? (join "|",@{$option->{-s}}) : '\w+';

	&output_extract_orgseq($option->{-d},$selinfo);

} elsif ($exemode eq 's') {
	&output_single_fasta($seqinfo);
} elsif ($exemode eq "c") {
	open (OUT,">rev_manuseq.fas");
	open (OUT2,">rev_rev_manuseq.fas");

	$lmlen= exists $option->{-lmlen} ? $option->{-lmlen}->[0] : 1;
	$seqlimno= 0;
    $Ambase_no=0;
    $Nbase_no=0;
	foreach $tax (@{$seqinfo->{tax}}) {
        $seq= $seqinfo->{seq}->{$tax};
        
        $amb= 0;
        $amb++ if($seq=~ /[^ATGCatgcNn]/);
        $Ambase_no++ if ($amb > 0);

        $Ns= 0;
        $Ns++ if($seq=~ /[Nn]/);
        $Nbase_no++ if ($Ns > 0);

        $seqlen= $seqinfo->{len}->{$tax};
        next if ($seqlen < $lmlen);
        $seqlimno++;
        
        $seqline= &limit_seqline($seq);
		$seqline=~ tr/U/T/;
		$seqline=~ tr/./-/;

		###  INDEX ADDITIO
		print OUT2 ">$tax  $seqinfo->{index}->{$tax}\n$seqline\n";
		print OUT ">$tax\n$seqline\n";
#		print "$tax $seqinfo->{len}->{$tax}\n";
	}

	print "FILTERED SEQ: $seqlimno\n";
    print "Ns CONTAIN SEQ: $Nbase_no\n";
    print "AMB CONTAIN SEQ: $Ambase_no\n";
	print "OUTPUT: rev_manuseq.fas\n";
	print "OUTPUT2: rev_rev_manuseq.fas  ## contain index\n";

} elsif ($exemode eq "d") {
	$unit= exists $option->{-u}  ? $option->{-u}->[0] : 400;
	$basename= $option->{-d}->[0];
	$basename= $1 if ($basename=~ /(\w+)\.\w+$/);
	($count,$subno)= (0,0);

	$lsub= exists $option->{-lsub} ? $option->{-lsub}->[0] : 0;
	$outdir= exists $option->{-outdir} ? $option->{-outdir}->[0] : "DIV";

	mkdir $outdir unless (-d $outdir);

	foreach $tax (@{$seqinfo->{tax}}) {
		$count++;
		if ($count%$unit == 1) {
			$subno++;
			$outfile= "div_${basename}_${subno}.fasta";

			open (OUT,">${outdir}/$outfile");
		}
		$seqline= &limit_seqline($seqinfo->{seq}->{$tax});
		print OUT ">$tax  $seqinfo->{index}->{$tax}\n$seqline\n";

		last if ($subno > $lsub && $lsub > 0);
	}

	print "SEQ_NUM $count\n";
	print "LNUM: $lsub  OUTD: $outdir\n";

} elsif ($exemode eq "r") {
	$infile= $option->{-d}->[0];
	$inname= (split/\//,$infile)[-1];
	open (OUT,">rev_${inname}");

	foreach $tax (@{$seqinfo->{tax}}) {
		$index= exists $seqinfo->{index} ? $seqinfo->{index}->{$tax} : "";
		$seq= $seqinfo->{seq}->{$tax};
		$seq=~ tr/\-.//d;
		$seq=~ tr/U/T/;
		$seq= &limit_seqline($seq);
		print OUT ">$tax $index\n$seq\n";
	}
	print "OUTFILE rev_${inname}\n";

} elsif ($exemode eq "x") {
### CUTTING
	$infile= $option->{-d}->[0];
	$inname= (split/\//,$infile)[-1];
	
	($cinfo)=&cutting_sequence($seqinfo,$inname,$option);
	print "CUT_READ: $cinfo->{'y'} \/ $cinfo->{'c'}\n";

	print "OUTPUT: rev_cut_seq.fas\n";


} elsif ($exemode eq "n") {
	$infile= $option->{-d}->[0];
	$inname= (split/\//,$infile)[-1];
	open (OUT,">rev_${inname}");
	foreach $tax (@{$seqinfo->{tax}}) {
		$seq= $seqinfo->{seq}->{$tax};
		$len= $seqinfo->{len}->{$tax};

#		if (substr($seq,-3,3)=~ /TGA|TAA|TAG/) {
#			$subseq= substr($seq,0,$len-3);
#		} else {
			$subseq= $seq;
#		}
#		$subseq= &limit_seqline($subseq);
		print OUT ">$tax\n$subseq\n";
	}
} elsif ($exemode eq "v") {

	open (OUT,">rename_seq.fas");
	&rename_sequence($seqinfo,$option->{-l}->[0]);


} elsif ($exemode eq "m") {
	open (OUT,">rev_stpseq.fas");

	($total)= &mince_sequence($seqinfo,$option->{-u}->[0]);

	print "SPLIT NUM ==> $total\n";

} elsif ($exemode eq "sc") {
	open (OUT,">rename_scseq.fas");
	open (LOG,">rename_scseq.logs");

	($total)= &create_super_contigs($seqinfo,$option->{-gap}->[0]);

    print "SC CONTIGS ==> Total: $total\nOUTPUT: rename_scseq.fas rename_scseq.logs\n";

} elsif ($exemode eq "a") {
	open (OUT,">re_seq.fas");
	open (LOG,">re_seq.logs");

	($total)= &output_complement_seq($seqinfo);

	print "Complement seq ==> OUTPUT: re_seq.fas Total: $total\n";

} elsif ($exemode eq "stat") {
	$width= exists $option->{-width} ? $option->{-width}->[0] : 20;

	($info)= &statics_readlength($option->{-d},$width);

	print "#READS: $info->{count}\nTOTAL_BASE: $info->{sum}\n";
	print "MIN: $info->{min} MAX: $info->{max} AVE: $info->{ave} \n";
	print "LEN_DIST 1\n";

	print "#### READ 1\nLEN\tRATE\tCUMRate\tCOUNT\n";
		$sum_distrate= 0;
		foreach $inter (sort{$a<=>$b} keys %{$info->{len}}) {
			$inter_count= int($info->{len}->{$inter}/$info->{count}*100);
			$inter_rate= sprintf("%.3f",$info->{len}->{$inter}/$info->{count}*100);
			$sum_distrate+= $info->{len}->{$inter}/$info->{count};
			print "$inter\t".('*' x$inter_count)."  $inter_rate\t".sprintf("%4f",$sum_distrate)."\t$info->{len}->{$inter}\n";
		}

} elsif ($exemode eq "sample") {
	$sample_number= exists $option->{-numb} ? $option->{-numb}->[0] : 100;

  open (OUT,">sampling_seq.fas");

  $sample_no= 0;
  foreach $taxone (keys %{$seqinfo->{seq}}) {
		$sample_no++;
		$sampseq= $seqinfo->{seq}->{$taxone};
		$subline= &limit_seqline($sampseq);
		print OUT ">$taxone\n$subline\n";

		last if ($sample_no >= $sample_number);
	}

} elsif ($exemode eq "len") {
	$lmin= exists $option->{-lmin} ? $option->{-lmin}->[0] : 100;
    $lmax= exists $option->{-lmax} ? $option->{-lmax}->[0] : 0;

    $getno=0;
  	open (OUT,">out_filter.fas");
  	foreach $taxone (keys %{$seqinfo->{seq}}) {
		if (length($seqinfo->{seq}->{$taxone}) > $lmin) {
            if ($lmax > 0) {
                if (length($seqinfo->{seq}->{$taxone}) <= $lmax) {
                    print OUT ">$taxone\n$seqinfo->{seq}->{$taxone}\n";
                    $getno++;
                }
            } else {
                print OUT ">$taxone\n$seqinfo->{seq}->{$taxone}\n";
                $getno++;
            }
		}
	}
	print "## GET: $getno\n";
	print "OUTPUT: out_filter.fas\n";
} elsif ($exemode eq "noint") {
    open (OUT,">out_noint.fas");
    foreach $taxone (keys %{$seqinfo->{seq}}) {
        print OUT ">$taxone\n$seqinfo->{seq}->{$taxone}\n";
    }
    
    print "OUTPUT: out_noint.fas\n";
    
}



##==========================================================================

sub statics_readlength {
	my ($infiles,$width)= @_;
	my ($info);

	$info= {};

	open (LOG,">out_manu_fasta_stat.logs");
	open (IN,"$$infiles[0]");
	$info->{count}= 0;
	$info->{min}= 1000;
	$info->{max}= 0;
	$info->{sum}= 0;
	@lines=();

	while(<IN>) {
		chomp;
		if (/^>(\S+)/) {
			$info->{count}++;
			if ($info->{count} > 1) {
				$seq= join "",@lines;
				$seqlen= length $seq;
				$info->{min}= $seqlen if ($seqlen < $info->{min});
				$info->{max}= $seqlen if ($seqlen > $info->{max});
				$info->{sum}+= $seqlen;
				$interval= int($seqlen/$width)*$width;
				$info->{len}->{$interval}++;
				
				@lines=();
			}
			&display_couter($info->{count},[100000,10,3],$info->{count});
		} else {
			push @lines,$_

		}

	}

	$info->{ave}= sprintf("%.2f",$info->{sum}/$info->{count});
	return ($info);

}




##==========================================================================
##  RENAME
##==========================================================================
sub output_complement_seq {
	my ($seqinfo)= @_;


	$total= 0;
	foreach $taxone (@{$seqinfo->{tax}}) {
		$rev_seq= &complement_seq($seqinfo->{seq}->{$taxone});
		$rev_seq_inter=&multiline_seq($rev_seq);
		$index= $seqinfo->{index}->{$taxone};
		print OUT ">$taxone $index\n$rev_seq_inter\n";
		$total++;


	}
	return ($total);
}


sub create_super_contigs {
	my ($seqinfo,$gaplen)= @_;

    $gapseq= $gaplen > 0 ? 'N' x$gaplen : '';
	$sc_contig_seq="";
	$taxnum= $seqinfo->{taxnum};
	$sc_contig_len= 0;

	foreach $taxone (@{$seqinfo->{tax}}) {
		$seqlen= $seqinfo->{len}->{$taxone};
		if ($taxnum > 1) {
			$sc_contig_seq= $sc_contig_seq.$seqinfo->{seq}->{$taxone}.$gapseq;
			$start= $sc_contig_len+1;
			$end= $sc_contig_len+$seqlen;
			print LOG "$taxone\t+\t$start\t$end\t$seqlen\n";
		} else {
			$sc_contig_seq= $sc_contig_seq.$seqinfo->{seq}->{$taxone};
            $start= $sc_contig_len+1;
            $end= $sc_contig_len+$seqlen;
			print LOG "$taxone\t+\t$start\t$end\t$seqlen\n";
		}
		$taxnum--;
		$sc_contig_len= length($sc_contig_seq);
	}

	$sequence= &limit_seqline($sc_contig_seq);
	print OUT ">sc_contigs\n$sequence\n";

	return ($sc_contig_len);

}



sub rename_sequence {
	my ($seqinfo,$listfile)= @_;
	my (@line,$info,$regname,$seqname,$seq,$newname);

	open (LI,"$listfile");
	@line=<LI>;
	foreach (@line) {
		chomp;
		next if (/^(\#|\s*$)/);
		@item= split/\t/;
		$oriname= shift @item;
		$newname= shift @item;
		$index= @item > 0 ? join " ",@item : "";
		$info->{each}->{$oriname}= [$newname,$index];
		push @{$info->{all}},$oriname;
	}

#	open (OUTS,">rename_seq.fas");

	$regname= join "|",@{$info->{all}};

	$count= 0;

	foreach $seqname (@{$seqinfo->{tax}}) {
		$index= $seqinfo->{index}->{$seqname};		
		$sequence= &limit_seqline($seqinfo->{seq}->{$seqname});

		if (exists $info->{each}->{$seqname}) {
			($newname,$newindex)= @{$info->{each}->{$seqname}};
			$index= $newindex if ($newindex eq "");
			$seqname= $newname;
			$count++;
		} elsif ($seqname=~ /($regname)/) {
			($newname,$newindex)= @{$info->{each}->{$1}};
			$index= $newindex if ($newindex eq "");
			$seqname= $newname;
			$count++;
		}

		print OUT ">$newname  $index\n$sequence\n";

	}

	print "RENAME $count\n";

}

##==========================================================================
##  SELSEQ FASTA MODE 1
##==========================================================================
sub extract_sel_seqname {
	my ($seqinfo,$selinfo,$exemode)= @_;
	my ($one,@keyones,$seqname,$index,$hit,$cat);

	if ($selinfo->{keynum} > 0) {
		foreach $one (sort keys %{$selinfo->{info}}) {
			@keyones= ();
			foreach (@{$selinfo->{info}->{$one}}) {
				s/\|/\\|/g;
				s/\./\\./g;
				push @keyones,$_;
			}
			
			if ($exemode eq "es") {
				if ($one == 1) {
					$selinfo->{key}->{$one}= '^('.(join '|',@keyones).')\:';
				} else {
					$selinfo->{key}->{$one}= (join '|',@keyones);
				}
			} else {
				$selinfo->{key}->{$one}= (join '|',@keyones);
			}
 		}
	}

	foreach $seqname (@{$seqinfo->{tax}}) {
	#	$index= $seqname;
		$index= "$seqname $seqinfo->{index}->{$seqname}";

		($hit,$cat)= (0,'TMP');
		if ($selinfo->{keynum} == 0) {
			push @{$extrainfo->{$cat}},$seqname;
			next;
		}
			
		foreach $reg_key_no (sort keys %{$selinfo->{key}}) {
			$reg_key= $selinfo->{key}->{$reg_key_no};
			$hit++ if ($index=~ /$reg_key/);
			$cat= $& if ($hit == 1 && $reg_key_no == 1);
		}
		
		push @{$extrainfo->{$cat}},$seqname if ($hit == $selinfo->{keynum});
	}

	return ($extrainfo);
}
sub output_extract_seq {
	my ($seqinfo,$extrainfo,$outmode)= @_;
	my ($org,$one,$seq,$index,$seqline);

	open (OUT, ">extra_seq.fas") if ($outmode eq "m");

	if ($outmode=~ /g|s/) {
		mkdir ("EACH") unless (-d "EACH");
	}
	my $tohit= 0;
	foreach $org (sort keys %{$extrainfo}) {
		open (OUT, ">./EACH/extra_$org.fas") if ($outmode eq "g");
		$orghit= 0;
		foreach $one (@{$extrainfo->{$org}}) {
			open (OUT, ">./EACH/${one}.fas") if ($outmode eq "s");
			if (exists $seqinfo->{seq}->{$one}) {
				$seq= $seqinfo->{seq}->{$one};
				$index= $seqinfo->{index}->{$one};
				$seqline= &limit_seqline($seq);
				print OUT ">$one  $index\n$seqline\n";
				$orghit++;
				$tohit++;
			}
		}
		print LOG "ORG $org HIT $orghit\n";
		print "HITNUMBER 1ST_GP $org HIT $orghit\n";
	}
		print "TOTAL_HITNUMBER $tohit\n";
}


#==========================================================================
# SELSEQ ECO MODE
#==========================================================================
sub output_extract_orgseq {
	my ($infiles,$selinfo)= @_;
	my ($info,$onefile,$line);

	$info= {};
	open (OUT,">ext_seq.fas");

	$info->{total}= 0;
	$info->{to_hit}= 0;

	foreach $onefile (@$infiles) {
		open (IN,"$onefile");

		while (<IN>) {
			chomp;
			next if (/^(\#|\s*$)/);
			$line= $_;
			if (/^>(\S+)/) {
				$seqname= $1;
				$org= $1 if ($seqname=~ /^(\w+)\:/);
				$info->{total}++;
				if ($org=~ /^($selinfo->{reg})$/) {
					$info->{org_hit}->{$org}++;
					$info->{to_hit}++;
					print OUT "$line\n";
					$info->{hit} = 1;
				} else {
					$info->{hit} = 0;
				}
				&display_counter($info->{total},[10000,10,3],"$info->{to_hit}/$info->{total}");
			} else {
				print OUT "$line\n" if ($info->{hit} == 1);
			}
		}
	}

	print "SEL $selinfo->{reg}\n HIT $info->{to_hit} / TO $info->{total}\n";

}
sub output_extract_seq_eco2 {
	my ($infiles,$selinfo,$exemode,$outmode)= @_;
	my ($onefile,$seqname,$string,$tax);

	($selinfo)= &create_rev_selinfo($selinfo,$exemode);

	my $extrainfo= {};
	my $seqinfo= {};
	$count= 0;
	$hcount= 0;
	foreach $onefile (@$infiles) {
		open (IN,"$onefile");

		while (<IN>) {
			chomp;
			next if (/^(\#|\s*$)/);

			if (/^>(\S+)/) {
				$seqname= $1;
				$string= $_;
				$string=~ s/^\>//;
				($extrainfo)= &check_keyword($string,$selinfo,$extrainfo);
				if ($extrainfo->{hit} == 1) {
					push @{$seqinfo->{tax}},$seqname;
					$string=~ s/^$seqname\s*//;
					$seqinfo->{taxnum}++;
					$seqinfo->{index}->{$seqname}= $string;
					$hcount++;
				}
				$count++;
				&display_counter($count,[10000,10,3],"$hcount/$count");

			} else {
				push @{$seqinfo->{seqline}->{$seqname}},$_ if ($extrainfo->{hit} == 1);
			}
		}
	}
	($seqinfo)= &join_seqline($seqinfo);

	open (OUT,">rev_seqseq.fas");

	foreach $tax (@{$seqinfo->{tax}}) {
		$seqlin= &limit_seqline($seqinfo->{seq}->{$tax});
		$index= $seqinfo->{index}->{$tax};
		print OUT ">$tax  $index\n$seqlin\n";
	}

}

sub check_keyword {
	my ($string,$selinfo,$extrainfo)= @_;
	my ($hit,$cat,$reg_key_no,$reg_key);

	$extrainfo->{hit}= 0;

	$extrainfo->{hit}++ if ($selinfo->{keynum} == 0);
	$cat= 'TMP' if ($selinfo->{keynum} == 0);
	
	foreach $reg_key_no (sort keys %{$selinfo->{key}}) {
		$reg_key= $selinfo->{key}->{$reg_key_no};
		$extrainfo->{hit}++ if ($string=~ /$reg_key/);
		$cat= $& if ($extrainfo->{hit} == 1 && $reg_key_no == 1);
	}
	if ($extrainfo->{hit} >= $selinfo->{keynum}) {
		push @{$extrainfo->{grp}->{$cat}},$seqname;
		$extrainfo->{hit}= 1;
	}

	return ($extrainfo);
}
sub create_rev_selinfo {
	my ($selinfo,$exemode)= @_;
	my ($one,@keyones);

	if ($selinfo->{keynum} > 0) {
		foreach $one (sort keys %{$selinfo->{info}}) {
			@keyones= ();
			foreach (@{$selinfo->{info}->{$one}}) {
				s/\|/\\|/g;
				s/\./\\./g;
				push @keyones,$_;
			}
			
			if ($exemode eq "es") {
				if ($one == 1) {
					$selinfo->{key}->{$one}= '^('.(join '|',@keyones).')\:';
				} else {
					$selinfo->{key}->{$one}= (join '|',@keyones);
				}
			} else {
				$selinfo->{key}->{$one}= (join '|',@keyones);
			}
 		}
	}

	return ($selinfo);
}

##==========================================================================
##  MINCE(SPLIT) FASTA
##==========================================================================

sub mince_sequence {
	my ($seqinfo,$unit)= @_;
	my ($total,$tax,$seq,$len,$c,$x,$getlen,$subseq,$new_tax,$subline);

#	$limlen= 500;
	$total= 0;
	foreach $tax (@{$seqinfo->{tax}}) {
		$seq= $seqinfo->{seq}->{$tax};
		$len= $seqinfo->{len}->{$tax};
		$c= 0;
      for($x= 0; $x < $len; $x += $unit) {
			if ($len-$x < $unit) {
				$getlen= $len-$x;
			} else {
				$getlen= $unit;
			}
#         next if ($limlen > $getlen);

			$subseq= substr($seq,$x,$getlen);
			$c++;
			$new_tax= "${tax}_$c";
			$subline= &limit_seqline($subseq);
			print OUT ">$new_tax\n$subline\n";
			$total++;
		}
	}
	return ($total);
	
}




##==========================================================================
##  SINGLE FASTA
##==========================================================================
sub output_single_fasta {
	my ($seqinfo)= @_;
	my ($to,$co,$seqlin,$org,$name);

	mkdir ("EACH") unless (-d "EACH");
	($to,$co)= (0,0);
	foreach $one (@{$seqinfo->{tax}}) {
		$to++;
	##	next if ($seqinfo->{len}->{$one} > 100);

		if ($one=~/\:/) {
			($org,$name)= split/\:/,$one;
		} else {
			$name= $one;
		}


		$co++;
		$seqlin= &limit_seqline($seqinfo->{seq}->{$one});
		open (SIG,">./EACH/$name.fas");
		$index= exists $seqinfo->{index}->{$one} ? $seqinfo->{index}->{$one} : "";
		print SIG ">$name  $index\n$seqlin\n";

	}
	print "TOTAL $to GET $co\n";

}

##==========================================================================
##
##==========================================================================
sub read_list {
	my ($list)= @_;
	my ($selinfo);

	open (LI,"$list");
	$selinfo= {};
	$selinfo->{keynum}= 0;

	while(<LI>) {
		chomp;
   	next if (/^(\#|\s*$)/);
		if (/>(\S+)/) {
			$selinfo->{keynum}++;
		} else {
			push @{$selinfo->{info}->{$selinfo->{keynum}}},split/\s+/;
		}
  }

  return ($selinfo);
}

sub cutting_sequence {
	my($seqinfo,$inname,$option)= @_;
### CUTTING

	open (OUT,">rev_cut_seq.fas");

	$cinfo= {};

	foreach $tax (@{$seqinfo->{tax}}) {
		$index= exists $seqinfo->{index} ? $seqinfo->{index}->{$tax} : "";
		$seq= $seqinfo->{seq}->{$tax};
		$len= $seqinfo->{len}->{$tax};
		$cinfo->{c}++;
		if (exists $option->{-a}) {
			($st,$ed)= @{$option->{-a}};
			$sublen= $ed-$st+1;
			if ($len > $st) {
				$subseq= substr($seq,$st-1,$sublen);
				print OUT ">$tax\n$subseq\n";
 			  $cinfo->{y}++ if ($len >= $ed);
			} 
		} elsif (exists $option->{-b}) {
			($p5cut,$p3cut)= @{$option->{-b}};
			if ($len > $p5cut + $p3cut) {
				$subseq_1= substr($seq,0,$len-$p3cut);
				$subseq_2= substr($subseq_1,$p5cut);
				$sublen= length $subseq_2;
				print OUT ">$tax\n$subseq_2\n";
				$cinfo->{y}++;
			}
		}
	}

	return($cinfo);

}



##==========================================================================
##
##==========================================================================

sub readseq_fasta {
	my ($seqfile)= @_;
	my ($count,$preseqinfo,$seqinfo,$dupinfo);
	my ($fileone,@lines,$subtotal,$seqname,@index,$noseqnum);
  
	##  CHECK DUPLICATION SEQ AND 0 LEN SEQ
	$count= 0;
	$preseqinfo={};
	$seqinfo={};
	$dupinfo= {};
	$dupinfo->{total}= 0;
	$noseqnum= 0;
	foreach $fileone (@$seqfile) {
		open (SEQ,"$fileone");
		@lines= <SEQ>;
		$subtotal= 0;
		foreach(@lines){
			chomp;
			next if (/^(\#|\s*$)/);        
			if (/^>(\S+)/){
				$count++;
				&display_counter($count,[10000,10,2],"$count");
				$seqname=$1;
				@index= split/\s+/;
				shift @index;
				if (exists $preseqinfo->{seqline}->{$seqname}) {
					$dupinfo->{$seqname}++;
					$dupinfo->{total}++ if ($dupinfo->{$seqname}==2);
				} else {
					push @{$preseqinfo->{tax}},$seqname;
					$preseqinfo->{taxnum}++;
					$preseqinfo->{index}->{$seqname}= join " ",@index;
					$dupinfo->{$seqname}= 1;
				}
				$subtotal++;
			}else{
				push @{$preseqinfo->{seqline}->{$seqname}},$_ if ($dupinfo->{$seqname} == 1);
			}
		}
		print "\nSUMMARY\n";
		foreach $seqname (@{$preseqinfo->{tax}}) {

			if (exists $preseqinfo->{seqline}->{$seqname}) {
				push @{$seqinfo->{tax}},$seqname;
				$seqinfo->{taxnum}++;
				$seqinfo->{index}->{$seqname}= $preseqinfo->{index}->{$seqname};
				$seqinfo->{seqline}->{$seqname}= $preseqinfo->{seqline}->{$seqname};
			} else {
				print "   NOSEQ $seqname\n";
				$noseqnum++;
			}							 
		}		
		print "FILE $fileone SEQ_NUMBER  S= $seqinfo->{taxnum} D= $dupinfo->{total} T= $subtotal\n";
	}
  
	my $dupl= 0;
	foreach $seqname (sort keys %{$dupinfo}) {
		next if ($seqname=~ /total/);
		$dupl++ if ($dupinfo->{$seqname} > 1);
		print "  $dupl DUPLI $seqname $dupinfo->{$seqname}\n" if ($dupinfo->{$seqname} > 1);
	}
	print "SEQ $seqinfo->{taxnum} DUPLI $dupinfo->{total} NOSEQ $noseqnum\n";

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
	#	$sequence=~ tr/U/T/;
		$seqinfo->{seq}->{$seqname}= $sequence;
		$seqinfo->{len}->{$seqname}= length $sequence;
		delete $seqinfo->{seqline}->{$seqname};
		# print "$seqname $seqinfo->{len}->{$seqname}\n";
	}
	return ($seqinfo);
}
###=====================================================================
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


#==========================================================================
#=============================DISPLAY COUNTER =============================
sub display_counter {
  my ($count,$limits,$sign)=@_;
  my ($lim);

  ($subunit,$munit,$rep)= (@$limits);

  print "\-" if ($count%$subunit == 0);
  print "$sign" if ($count%($subunit*$munit) == 0);
  print "\n" if ($count%($subunit*$munit*$rep) == 0);
}


sub complement_seq {
  my ($seq)= @_;
  my ($revseq);
  
  $seq=~ tr/ATGCatgc\-NBDHVRYKMSW/TACGtacg\-NVHDBYRMKSW/;
  $revseq = reverse($seq);
  return $revseq;
}

sub display_couter {
  my ($count,$limits,$sign)=@_;
  my ($lim);

  ($subunit,$munit,$rep)= (@$limits);

  print "\-" if ($count%$subunit == 0);
  print "$sign" if ($count%($subunit*$munit) == 0);
  print "\n" if ($count%($subunit*$munit*$rep) == 0);
}

