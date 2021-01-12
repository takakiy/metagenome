#!/usr/bin/perl -w

#====================================  OPTIONS  =================================
## REV. 2015.11.12

if ($ARGV[0]=~/\-h/){
     print "-a GENE_ARRAY -d Database -t c\n";
     print "-l list(gene) -d Database -t c\n";
     print "-lo ortholog table -d Database\n";
     print "OPTIONS
         -t c: column / (r: row) / a: all 
         -m T/(F): MultiPattern
         -o (s; singleFile) /m: multiFatsa
         -c (e: each) / c: concatemer
         -n (ori)/same/new: NameConvert  (same is Only selected -t r)
         -w T/(F)  (Only top genes/ all genes in one box) ## -n same is -w T
         -info (F)/T\n";
     print "LIST; \>ID\n    Gene1 Gene2 Gene3\n";
     print "ORTH; ORID XX XX bcl	bha bpf O_NUM\n";
     exit;
}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

if (exists $option->{-lo}) {
  $unitmode= exists $option->{-t} ? $option->{-t}->[0] : "r";
} else {
  $unitmode= exists $option->{-t} ? $option->{-t}->[0] : "c";
}

$option->{-w}->[0]= exists $option->{-w} ? $option->{-w}->[0] : "F";
##===============================================================================
##  CREATE GENELIST
##===============================================================================
if (exists $option->{-a}) {
   $getgene_info->{extra}->{zzz}=[@{$option->{-a}}];
   ($orgs)= &extract_orgs($getgene_info);
} elsif (exists $option->{-l}) {
   open ("LI","$option->{-l}->[0]");
   @lines=<LI>;
   $getid= "getseq";
   foreach (@lines){
      chomp;
      next if (/^(\#|\s*$)/);
      if (/^\>(\S+)/) {
         $getid= $1;
      } else {
         push @{$getgene_info->{$getid}->{zzz}},grep{/\w+/}(split/\s+/);
         
      }
   }
   ($orgs)= &extract_orgs($getgene_info);
} elsif (exists $option->{-lo}) {
    $selorth= exists $option->{-selorth} ? $option->{-selorth}->[0] : '';
    
   ($orgs,$orthinfo)= &create_orthlog_info($option->{-lo}->[0],$selorth);
   if (exists $option->{-m}) {
     ($orthinfo)= &convert_multiset($orgs,$orthinfo) if ($option->{-m}->[0] eq "T");
   }
   if ($unitmode eq "r") {
     ($getgene_info)= &create_getgeneinfo_row($orgs,$orthinfo);
   } elsif ($unitmode eq "c") {
     ($getgene_info)= &create_getgeneinfo_column($orgs,$orthinfo);   
   } elsif ($unitmode eq "a") {
     ($getgene_info)= &create_getgeneinfo_all($orgs,$orthinfo);
   }
}

##  GLOBAL $orthinfo
##===============================================================================
##  CREATE SEQUENCE DATABASE
##===============================================================================
$|=1;

print "SELECT ORG @$orgs\n START CREATE SEQ DATABASE $option->{-d}->[0]\n";

## ($seqinfo)= &create_sequence($option->{-d},$orgs,$getgene_info);
($seqinfo)= &readseq_fasta($option->{-d});
print "END SEQ DATABASE\n";

print "$seqinfo->{tax}->[0]xxx\n";


##===============================================================================
##  OUTPUT 
##===============================================================================
use File::Basename;
($dbbase,$dbdirect,$dbsuf)= fileparse($option->{-d}->[0],('.nuc','.pep','.fas','.fna','.fa','.faa','.fasta'));
print "DBASE $dbbase $dbdirect $dbsuf\n";

if ((exists $option->{-l}) || (exists $option->{-lo})) {
   $l_base= basename($option->{-l}->[0]) if (exists $option->{-l});
   $l_base= basename($option->{-lo}->[0]) if (exists $option->{-lo});
} else {
   $l_base= "extra";
}

$outfile_mode= exists $option->{-o} ? $option->{-o}->[0] : "s";
if ($outfile_mode eq "s") {
    system ('mkdir EXT') unless (-e "EXT") ;
}

$output_mode= exists $option->{-c} ? $option->{-c}->[0] : "e";
$name_mode= exists $option->{-n} ? $option->{-n}->[0] : "ori";
if ((exists $option->{-a}) || (exists $option->{-l})) {
  $des_mode= exists $option->{-info} ? $option->{-info}->[0] : "T";
} else {
  $des_mode= exists $option->{-info} ? $option->{-info}->[0] : "F";
}

open (ER,">getgene.logs");
$options=[$outfile_mode,$output_mode,$name_mode,$des_mode,$l_base,$dbsuf];
print "@$options\n";
&control_output($orgs,$getgene_info,$seqinfo,$options);

print "OPTIONS $_ $option->{$_}->[0]\n" foreach (sort keys %{$option});

##===============================================================================
##  SUB ROUTING
##===============================================================================
sub extract_orgs {
   my ($getgene_info)=@_;
   my ($getid,%hasorg,@orgs);
   %hasorg=();
   foreach $getid (sort keys %{$getgene_info}) {
      foreach (@{$getgene_info->{$getid}->{zzz}}) {
         $hasorg{$1}=1 if (/^([a-z0-9]{3,})\:/);
      }
   }
   @orgs= sort keys %hasorg;
   @orgs=("zzz") if (@orgs == 0);
   return ([@orgs]);
}
##===============================================================================
sub create_orthlog_info {
  my ($orthfile,$selorth)=@_;
  my ($orthinfo,@lines,@items,@orgs,$precol,$col,$catcol,$orgc,$orgnum);
  my ($hpat,$hnum);
    
    $sel= 0;
    $selorths= {};
    if ($selorth ne '') {
        $sel= 1;
        open(OO,"$selorth");
        while (<OO>) {
            chomp;
            $selorths->{$_}++ foreach (split);
        }
    }
    
    
  open (ORF,"$orthfile");
  @lines=<ORF>;
  $orthinfo={};
  $rownum= 0;
  foreach (@lines){
    chomp;

    @items= split/\t/;
    if (/\#\#\s*(a_ortholog|ORID)/) {
       @orgs=();
       $precol=0; $orgnum=0;

		for ($col=0; $col<=$#items; $col++) {
			if ($items[$col]=~ /^Comp|Num$/) {
				## neglect
             $precol++ if ($orgnum==0);

			} elsif ($items[$col] eq "O_NUM") {
				$catcol= $col;
				last;
			} elsif ($items[$col]=~ /^[a-zA-Z\-\d]{3,}$/) {
             $orgnum++;
             push @orgs,$items[$col] ;
			} else {
             $precol++ if ($orgnum==0);
			}
		}
      print "PRE: $precol $catcol $orgnum\n";
    } elsif (/^(\#|\s*$)/) {
      next; 
    } else {
        $orthid= $items[0];
        
        if ($sel == 1) {
            next unless (exists $selorths->{$orthid});
        }
        
       for ($orgc= $precol; $orgc < $catcol; $orgc++) {
          if ($items[$orgc] ne "") {
             $orthinfo->{$items[0]}->{each}->{$orgs[$orgc-$precol]}=
                     [split/\s+/,$items[$orgc]];
          }
       }
       ($hpat,$hnum)= &make_hpat_hnum([@orgs],$orthinfo->{$items[0]});
       $orthinfo->{$items[0]}->{hpat}= $hpat;
       $orthinfo->{$items[0]}->{hnum}= $hnum;
       $orthinfo->{$items[0]}->{cat}= $items[$catcol] if (defined $catcol);
       
       $rownum++;
    }
  }
  
  print "ORG_NUM $orgnum ORTH_NUM $rownum\n";
  return ([@orgs],$orthinfo);
}

sub make_hpat_hnum {
  my ($orgs,$membinfo)=@_;
  my (@hpats,@hnums,$pat,$num);
  
  @hpats=("C");  @hnums=("N");
  foreach (@$orgs){
     $pat= exists $membinfo->{each}->{$_} ? 1 : 0;
     $num= $pat==1 ? @{$membinfo->{each}->{$_}} : 0;
     push @hpats,$pat;
     push @hnums,$num;
  }
  return ((join "",@hpats),(join "_",@hnums));
}
##===============================================================================
sub convert_multiset {
   my ($orgs,$orthinfo)=@_;
   
   foreach $orthid (sort keys %{$orthinfo}) {
     ($multiset)= &create_multiset($orthinfo->{$orthid});
     next if (@$multiset < 2);
     $subnum=1;
     foreach $setone (@$multiset) {
       $orthsubid= $orthid."_".$subnum;
       $orthinfo->{$orthsubid}= $setone;
       ($hpat,$hnum)= &make_hpat_hnum($orgs,$setone);
       $orthinfo->{$orthsubid}->{hpat}= $hpat;
       $orthinfo->{$orthsubid}->{hnum}= $hnum;
       $subnum++;
##       print "$orthsubid $hpat $hnum\n";
     }
      delete $orthinfo->{$orthid}
  }
  return ($orthinfo);
}
sub create_multiset {
   my ($membinfo)=@_;
   my (@tempset,@orgset,$org,$setnum);
   my (@multiset,$setone,$addone);
   
   @tempset=();
   @orgset= sort keys %{$membinfo->{each}};
   $org= shift @orgset;
   $setnum=0;
   foreach (@{$membinfo->{each}->{$org}}) {
     $tempset[$setnum]->{each}->{$org}=[$_];
     $setnum++;
   }
   foreach $org (@orgset) {
      $setnum=0;  @multiset=();
      foreach $setone (@tempset) {
         foreach $addone (@{$membinfo->{each}->{$org}}) {
            $multiset[$setnum]->{each}->{$_}= $setone->{each}->{$_} foreach (keys %{$setone->{each}});
            $multiset[$setnum]->{each}->{$org}= [$addone];
            $setnum++;
         }
      }
      @tempset= @multiset;
   }
   return ([@multiset]);
}
##===============================================================================
sub create_getgeneinfo_row {
   my ($orgs,$orthinfo)=@_;
   my ($getgene_info,$org,$orthid);
   $getgene_info={};
   foreach $orthid (sort keys %{$orthinfo}) {
      foreach $org (@$orgs) {
         if (exists $orthinfo->{$orthid}->{each}->{$org}) {
            if ($option->{-w}->[0] eq "T") {
               $getgene_info->{$orthid}->{$org}= [$orthinfo->{$orthid}->{each}->{$org}->[0]];
            } else {
               $getgene_info->{$orthid}->{$org}= $orthinfo->{$orthid}->{each}->{$org};
            }
         }
      }
   }
   return ($getgene_info);
}

sub create_getgeneinfo_column {
   my ($orgs,$orthinfo)=@_;
   my ($getgene_info,$org,$orthid,@orgmemb);
   
   $getgene_info={};
   foreach $org (@$orgs) {
      foreach $orthid (sort keys %{$orthinfo}) {
         if (exists $orthinfo->{$orthid}->{each}->{$org}) {
            if ($option->{-w}->[0] eq "T") {
              $getgene_info->{$org}->{$orthid}= [$orthinfo->{$orthid}->{each}->{$org}->[0]];
            } else {
              $getgene_info->{$org}->{$orthid}= $orthinfo->{$orthid}->{each}->{$org};
            }
         }
      }
   }
   return ($getgene_info);
}
sub create_getgeneinfo_all {
   my ($orgs,$orthinfo)=@_;
   my ($getgene_info,$org,$orthid);
   
   $getgene_info={};
   foreach $orthid (sort keys %{$orthinfo}) {
      foreach $org (@$orgs) {
         if (exists $orthinfo->{$orthid}->{each}->{$org}) {
            foreach (sort @{$orthinfo->{$orthid}->{each}->{$org}}) {
               $getgene_info->{$_}->{$orthid}= [$_];
            }
         }
      }
   }
   return ($getgene_info);
}
##===============================================================================
#  OUTPUT 
##===============================================================================
sub control_output {
   my ($orgs,$getgene_info,$seqinfo,$options)=@_;
   my ($outfile_mode,$output_mode,$name_mode,$des_mode,$base,$suffix)= @$options;
   my ($outfile,$mainunit,@subunits,@err_genes,$que_num);
   my ($subunit,@membs,$seqid,$sequence,$seqname);
   $suffix= $suffix =~ /\./ ? $suffix : '.fas';
   if ($outfile_mode eq "m") {
      $outfile= $base.$suffix;
      unlink $outfile if (-f $outfile);
      open (OUT, ">>$outfile");
       
       print "OUTPUT: $outfile\n"
   }
   @err_genes= ();
   $que_num= 0;
   foreach $mainunit (sort keys %{$getgene_info}) {
      if ($outfile_mode eq "s") {
        $outfile= $mainunit.$suffix;
        open (OUT,">\./EXT/$outfile");
          
      }
      print OUT "\>$mainunit\n" if ($output_mode eq "c");
      
      @subunits= $unitmode eq "r" ? @$orgs : (sort keys %{$getgene_info->{$mainunit}});
      
      foreach $subunit (@subunits) {
##   print "EE $mainunit $mainunit $subunit\n";

         if (exists $getgene_info->{$mainunit}->{$subunit}) {
            @membs= @{$getgene_info->{$mainunit}->{$subunit}};
         } else {
            next;
         }
         foreach $seqid (@membs) {

            $que_num++;
 ###           unless (exists $seqinfo->{$seqid}) {
            unless (exists $seqinfo->{seq}->{$seqid}) {
               print ER "ERR NOSEQ $mainunit $subunit $seqid\n";
               push @err_genes,$seqid;
               next;
            }
 ###           $sequence= &limit_seqline($seqinfo->{$seqid}->{seq});
            $sequence= &multiline_seq($seqinfo->{seq}->{$seqid});
            if ($output_mode eq "e") {
               $seqname= &name_convert($name_mode,$seqid,$mainunit,$subunit);
 ###            $desc= $seqinfo->{$seqid}->{desc};
               $desc= $seqinfo->{index}->{$seqid};
						if ($desc eq "") {
               		print OUT "\>$seqname\n$sequence\n";
						} else {
               		print OUT "\>$seqname $desc\n$sequence\n";
						}
            } elsif ($output_mode eq "c") {
               print OUT "$sequence\n";
            }
         }
      }
   }
   $err_num= @err_genes;
   print "\#\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\=\n";
   print "ERR $err_num QUERY_NUM $que_num\n";
   
   
}

sub name_convert {
   my ($name_mode,$seqid,$mainunit,$subunit)=@_;
   my ($seqname,$cat);
   
   if ($name_mode eq "ori") {
      $seqname= $seqid;
   } elsif ($name_mode eq "same") {
      $seqname= $subunit;
   } elsif ($name_mode eq "new") {
      $cat= $orthinfo->{$mainunit}->{cat};
      $seqname= $cat.$subunit."_".$mainunit;
   }
   
   return ($seqname);
}


##===============================================================================
#  EXTRACT SEQUENCE
##===============================================================================

sub create_sequence {
  my ($seqfiles,$orgs,$getgene_info) = @_;
  my ($seqinfo,$seqc,$seqid,@descs,$seqfileone);
  my ($selectreg,$select);
 
  $seqinfo={};
  $selectreg= join "\|",@$orgs;

	my @getgeneall= ();
	foreach $one (sort keys %{$getgene_info}) {
		foreach $two (sort keys %{$getgene_info->{$one}}) {
			push @getgeneall,(join "|",@{$getgene_info->{$one}->{$two}});
		}
	}
	
	my $getgene_reg= join "|",@getgeneall;
	my $hitc= 0;

	foreach my $seqfileone (@$seqfiles) {
		open (SEQ,$seqfileone);
		$seqc=0;
		print "$seqfileone\n"; 
		while (<SEQ>) {
			chomp;
			next if (/^(\#|\s*$)/);
			if (/^>(\S+)/){
				$seqid= $1;
				$select= 0;
				$select= $seqid=~ /^($selectreg)\:/ ? 1 : 0;
				$select= 1 if ($selectreg eq "zzz");

				$select= $seqid=~ /^($getgene_reg)$/ ? 1 : 0;

				if ($select == 1) {
					$seqinfo->{$seqid}->{desc}="";
					@descs= split/\s+/;
					$seqinfo->{$seqid}->{desc}= "  ".(join " ",@descs[1..$#descs]) if (@descs>1);
					$hitc++;
				}
				&display_counter($seqc,[50000,10,3],$seqc);
				$seqc++;
       
			} else {
				## tr/[a-z]/[A-Z]/;
				push @{$seqinfo->{$seqid}->{seqline}},$_ if ($select == 1);
			}
		}
#		close (SEQ);
		print "\n";

  }
  print "$seqc HIT $hitc\n";


  foreach $seqid (sort keys %{$seqinfo}){
    if (exists $seqinfo->{$seqid}->{seqline}) {
        $seqinfo->{$seqid}->{seq}= join "",@{$seqinfo->{$seqid}->{seqline}};
    } else {
        print "ERROR: $seqid  No Seq\n";
    }
    delete $seqinfo->{$seqid}->{seqline};
  }

  return ($seqinfo);
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


sub readseq_fasta {
	my ($seqfile)= @_;
	my ($seqinfo,$seqname,$fileone,@lines,$index,$dupinfo,$onefinfo);
  
	$seqinfo={};
	$dupinfo= {};
	$seqc= 0;

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
				$seqc++;
				&display_counter($seqc,[50000,10,3],$seqc);
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
	#	$sequence=~ tr/U/T/;
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


#==========================================================================















