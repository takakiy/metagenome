#! /usr/bin/perl -w
# GeneBank
##  update 07/12/26


if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
  print "ARG:  -i gbk_file  -o outputID\n";
  exit;
}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

##========================================================
use File::Basename;



&main_control($option);

##============================================================================
##  CREATE ENTRYINFO
##============================================================================
sub main_control {
	my ($option)= @_;
	my ($infile,$tempname,$symb,$entinfo,$entry);

	$infile= $option->{-i}->[0];
print "$infile\n";
	$tempname= $1 if ($infile=~ /\/*([\w\d_\-]+)(\.\w+|)$/);

	$pdir= dirname $option->{-i}->[0];
	print "$pdir\n";


	$symb= exists $option->{-o} ? $option->{-o}->[0] : $tempname;
	print "INFILE\= $tempname SYMB\= $symb\n";

	($entinfo)= &create_entry_info($infile);

    #    $output="${pdir}\/${symb}";
    $output="${symb}";
    open (ONUC,">${output}\.nt");
    open (OGENE,">${output}");
    open (OGP,">${output}_cds.pep");
    
    open (LOG,">${output}_reform.logs");

    
	$fh_nuc= *ONUC;
	$fh_gene= *OGENE;
	$fh_pep= *OGP;

	foreach $entry (@{$entinfo->{entry}}) {
        #print "ENTRY $entry\n";
		&control_output($symb,$pdir,$entry,$entinfo->{info}->{$entry},$fh_nuc,$fh_gene,$fh_pep);

	}


}
##============================================================================
##  CREATE ENTRYINFO
##============================================================================

sub create_entry_info {
	my ($infile)= @_;
	my ($entinfo,@line,$entry,$info,$item,$oneval,$subinfo,$entno);
	
	open (IN,"$infile");
	@line=<IN>;
	$entno= 0;

	$entinfo= {};

	foreach (@line) {
 		chomp;
		next if (/^(\#|\s*$)/);

 		if (/^LOCUS/../^\/\//) {

			if (/^LOCUS\s+(\S+)/) {
				$entry= $1;
				$info= {};
				push @{$entinfo->{entry}},$entry;
			}
			$onelin= $_;

			if ($onelin=~ /^\/\//) {

 				foreach $item (@{$info->{item}}) {
														
					foreach $one_val (@{$info->{val}->{$item}}) {
						($subinfo)= &create_subinfo_2nd($one_val,$item);
						push @{$info->{info}->{$item}},$subinfo;
 					}

 				}

				delete $info->{val};
				$entinfo->{info}->{$entry}= $info;

			} else {
				if ($onelin=~ /^([A-Z_]+)\s{1,}/) {
					$item= $1;
					$onelin=~ s/^$item\s+/  TOP  /;
					push @{$info->{item}},$item unless (exists $info->{val}->{$item});
					push @{$info->{val}->{$item}},[];
				}
				push @{$info->{val}->{$item}->[-1]},$onelin;
			}
		}

	}

	return ($entinfo);
}

sub create_subinfo_2nd {
	my ($one_val,$item1st)= @_;
	my ($info,$onelin,$item,$one_val_2,$subinfo);

	$info= {};

	foreach $onelin (@$one_val) {
		if ($onelin=~ /^\s{2,}(\w+)\s{2,}/) {
			$item= $1;
			$onelin=~ s/^\s{2,}(\w+)\s{2,}/      \/TOP\=/;
			push @{$info->{item}},$item unless (exists $info->{val}->{$item});
			push @{$info->{val}->{$item}},[];
		}
		push @{$info->{val}->{$item}->[-1]},$onelin;
	}

	foreach $item (@{$info->{item}}) { 						
 										
		foreach $one_val_2 (@{$info->{val}->{$item}}) {
			($subinfo)= &create_subinfo_3rd($one_val_2,$item1st,$item);
			push @{$info->{info}->{$item}},$subinfo;

	#		@subitem= @{$info->{info}->{$item}->[0]->{item}};
	#		print  "$_ ==> $info->{info}->{$item}->[0]->{form}->{$_}->[0]\n"  foreach (@subitem);
		}
	}
	
	delete $info->{val};

	return ($info);	

}

sub create_subinfo_3rd {
	my ($one_val,$item1st,$item2nd)= @_;
	my ($info,$onelin,$item,$one_val_2,$form);

	$info= {};

	foreach $onelin (@$one_val) {
		if ($onelin=~ /^\s{6,}\/(\w+)\=?/) {
			$item= $1;
			$onelin=~ s/^\s{6,}\/$item\=?//;
			push @{$info->{item}},$item unless (exists $info->{val}->{$item});
			push @{$info->{val}->{$item}},[];

		} elsif ($onelin=~ /^\s{6,}/) {
			$onelin=~ s/^\s{6,}//;
		}
		push @{$info->{val}->{$item}->[-1]},$onelin;
	}
	
	foreach $item (@{$info->{item}}) { 						
 										
		foreach $one_val_2 (@{$info->{val}->{$item}}) {
			$form= (join " ",@$one_val_2);
			$form=~ s/^\"//; 
			$form=~ s/\"$//; 
			push @{$info->{form}->{$item}},$form;
		}
	}
	
	delete $info->{val};

	return ($info);
}
##============================================================================
##  CREATE ENTRYINFO
##============================================================================
sub reform_class_locus {
	my ($subinfo)= @_;
	my ($entry,$len,$mol,$mtype,$date);

	$form= $subinfo->{info}->{TOP}->[0]->{form}->{TOP}->[0];

	($entry,$len,$mol,$mtype,$date)= (split/\s+/,$form)[0,1,3,4,6];

	return [$entry,$len,$mol,$mtype,$date];
}

sub reform_class_organism {
	my ($subinfo)= @_;
	my ($source,$organism,$orgname,$orgline);

	$source= $subinfo->{info}->{TOP}->[0]->{form}->{TOP}->[0];

	$organism= "";
	if (exists $subinfo->{info}->{ORGANISM}) {
		$organism= $subinfo->{info}->{ORGANISM}->[0]->{form}->{TOP}->[0];
	}
#	print "$organism\n";
	($orgname,$orgline)= ("","");
	($orgname,$orgline)= ($1,$2) if ($organism=~ /^([^\;]+)\s(\w+\; .+)$/);

	return ($orgname,$orgline);
}

sub reform_class_multiline {
	my ($subinfo)= @_;
	my (@subitems,@allform,$subitem,$reform);

	@subitems= qw (AUTHORS TITLE JOURNAL PUBMED);
	
	@allform= ();
	foreach $subitem (@subitems) {
		if (exists $subinfo->{info}->{$subitem}->[0]->{form}->{TOP}->[0]) {
			push @allform,$subinfo->{info}->{$subitem}->[0]->{form}->{TOP}->[0];
		} else {

		}
	}

	$reform= join " ",@allform;
	
	return $reform;

}

sub reform_class_oneline {
	my ($subinfo)= @_;
	my ($reform);

	$reform= $subinfo->{info}->{TOP}->[0]->{form}->{TOP}->[0];
	
	return $reform;
}


sub reform_subclass_multiline {
	my ($subinfo)= @_;
	my (@subitems,@allform,$subitem,$reform);

	@subitems= qw (AUTHORS TITLE JOURNAL PUBMED);
	
	@allform= ();
	foreach $subitem (@subitems) {
		if (exists $subinfo->{info}->{$subitem}->[0]->{form}->{TOP}->[0]) {
			push @allform,$subinfo->{info}->{$subitem}->[0]->{form}->{TOP}->[0];
		} else {

		}
	}

	$reform= join " ",@allform;
	
	return $reform;

}

##============================================================================
##  OUTPUT DATA 1
##============================================================================

sub output_test {
	my ($symb,$entry,$info)= @_;

	foreach $item (@{$info->{item}}) {
		last if ($item eq "ORIGIN");
		foreach $subinfo (@{$info->{info}->{$item}}) {			
			foreach $subitem (@{$subinfo->{item}}) {
				foreach $subsubinfo(@{$subinfo->{info}->{$subitem}}) {
					foreach $subsubitem (@{$subsubinfo->{item}}) {						
						foreach $form (@{$subsubinfo->{form}->{$subsubitem}}) {							
					#		print  "$item $subitem $subsubitem\n  ==> $form\n";
						}
					}
				}
			}
		}
	}		

}

##============================================================================
##  OUTPUT DATA 2
##============================================================================

sub control_output {
	my ($symb,$pdir,$entry,$info,$fh_nuc,$fh_gene,$fh_pep)= @_;
	my (@baseform,$loc_inform,$definit,$acc,$orgname,$orgline,$topref,$baseform,$nuc,$gnum);


	@baseform= ();
	$loc_inform= &reform_class_locus($info->{info}->{LOCUS}->[0]);
	$definit= &reform_class_oneline($info->{info}->{DEFINITION}->[0]);
    
    ($acc,$orgname,$orgline,$topref)= ("","","","");

    $acc= &reform_class_oneline($info->{info}->{ACCESSION}->[0]);

    ($orgname,$orgline)= &reform_class_organism($info->{info}->{SOURCE}->[0]);

    $topref= &reform_class_multiline($info->{info}->{REFERENCE}->[0]);

    $loc_inform_txt= join " ",@{$loc_inform}[0..3];
    $baseform= "\#LOCUS=\t$loc_inform_txt\n\#DEF=\t$definit\n\#ACC=\t$acc\n";
    $baseform= $baseform."\#ORG=\t$orgname\n\#LINE\=\t$orgline\n\#REF=\t$topref\n";

    $nuc= $info->{info}->{ORIGIN}->[0]->{info}->{TOP}->[0]->{form}->{TOP}->[0];
	$nuc=~ tr/[0-9 ]//d;
	$nuc= &multiline_seq($nuc);
	print $fh_nuc ">$$loc_inform[0] $symb $$loc_inform[1]\n$nuc\n";

	($gnum)= &control_FeatureTable($symb,$pdir,$baseform,$info->{info}->{FEATURES}->[0],$fh_gene,$fh_pep,$$loc_inform[0]);

}

##============================================================================
##  OUTPUT FEATURE
##============================================================================
sub control_FeatureTable {
	my ($symb,$pdir,$baseform,$subinfo,$fh_gene,$fh_pep,$locus)= @_;
	my ($organism,$isolation,$clone);
	my ($subitem_txt,$all_informs_ord,@all_informs,$allnum,$scale,$gno_ord);
	my ($gnum,$subsubinfo,$gno,$type,$posit,$geneforms,$translation);


    $clone= "";
    if (exists $subinfo->{info}->{source}) {
        ($organism,$isolation,$clone)= &reform_source($subinfo->{info}->{source}->[0]);
    }
    $baseform= $baseform."\#CLONE\=\t$clone\n";
    ##    print OGENE "$baseform\n";
    print LOG "$baseform\n";

    my @xx= keys %{$subinfo->{info}};
    #    print "@xx\n";  ## source TOP CDS
    my @yy=@{$subinfo->{item}};
    #    print "Y @yy\n";  ## TOP source
    my @zz=keys %{$subinfo->{info}->{source}->[0]->{form}};
    #    print "Z @zz\n";  ## mol_type strain TOP organis
    my $cc= $subinfo->{info}->{source}->[0]->{form}->{TOP}->[0];
    #    print "CC $cc\n";  ## 1..250
    my $clen= $1 if ($cc =~ /1\.\.(\d+)/);

    $geneidinfo= {};
    
	$subitem_txt= join "\|",(grep{/(CDS|RNA)/}@{$subinfo->{item}});

	if ($subitem_txt=~ /CDS|RNA/) {
		@title= ("#CONTIGS","LEN","#TYPE","ID","GENE","DIR","ST","ED","PRODUCT","NOTE","CAT","FIX",
             "EC","NOTE2","AA","locus_tag","protein_id","gene_no");
		print $fh_gene (join "\t",@title)."\n";

	#	open (OGP,">${pdir}\/${symb}_cds.pep") if ($subitem_txt=~ /CDS/);

		($all_informs_ord)= &create_allgene_set($subinfo);
		@all_informs= @$all_informs_ord;

		$allnum= @all_informs;
		$scale= int(log($allnum)/log(10));
		$gno_ord= 0 x$scale;

		$gnum= 0;
        

        
		foreach $subsubinfo (@all_informs) {
			$gnum++;
			$gno= $symb."_".$locus."_".substr($gno_ord.$gnum,length($gno_ord.$gnum)-$scale-1,$scale+1);

			$type= $subsubinfo->{form}->{TYPE}->[0];
			$posit= $subsubinfo->{form}->{TOP}->[0];
            #          print "C $gno $type $posit\n";

			($geneforms,$translation,$geneidinfo)= &reform_geneform($gno,$subsubinfo,$geneidinfo);

            $glen= 0;
            if ($$geneforms[0] =~ /CDS/) {
                $glen= ($$geneforms[5]-$$geneforms[4]-2)/3;
            } elsif ($$geneforms[0] =~ /RNA/) {
                $glen= $$geneforms[5]-$$geneforms[4]+1;
            }
			print $fh_gene (join "\t",($locus,$clen,@$geneforms))."\n";

			if ($translation=~ /\w+/) {
				$translation=~ s/\s+//;
				$translation= &multiline_seq($translation);
				print $fh_pep ">$$geneforms[1] $$geneforms[6]\n$translation\n";
			}
		}
	}
	
	return ($gnum);
}
##============================================================================
sub reform_source {
	my ($subsubinfo)= @_;
	my (@subsubitem,$organism,$isolation,$clone);

	@subsubitem= qw(organism isolation_source clone);

	$organism= exists $subsubinfo->{form}->{organism} ? $subsubinfo->{form}->{organism}->[0] : "";
	$isolation= exists $subsubinfo->{form}->{isolation_source} ? 
		$subsubinfo->{form}->{isolation_source}->[0] : "";
	$clone= exists $subsubinfo->{form}->{clone} ? $subsubinfo->{form}->{clone}->[0] : "";

	return ($organism,$isolation,$clone);

}

##============================================================================
sub create_allgene_set {
	my ($subinfo)= @_;
	my (@subitems,$geneinfo,$subsubinfo,$posit,$positform);
	my (@all_informs,$subitem,@all_informs_ord);

	@subitems= grep{/(CDS|RNA)/}@{$subinfo->{item}};

    ## COLLECT ALL GENE
	$geneinfo= {};
	foreach $subsubinfo (@{$subinfo->{info}->{gene}}) {
		$posit= $subsubinfo->{form}->{TOP}->[0];
		$positform= &reforme_position($posit);
		$subsubinfo->{form}->{TYPE}->[0]= "gene";
		$subsubinfo->{form}->{START}->[0]= $$positform[0];
        
        if ($$positform[2] eq '+') {
            $n_posit=$$positform[0]."..".$$positform[1];
        } elsif ($$positform[2] eq '-') {
            $n_posit='complement('.$$positform[0]."..".$$positform[1].')';
        }
        
		$geneinfo->{info}->{$n_posit}= $subsubinfo;
        # print "BB $posit\n";
	}
    ## COLLECT CDS/RNA & delete gene info

	@all_informs= ();
	foreach $subitem (@subitems) {
		foreach $subsubinfo (@{$subinfo->{info}->{$subitem}}) {
			$posit= $subsubinfo->{form}->{TOP}->[0];
			$positform= &reforme_position($posit);
            
            if ($$positform[2] eq '+') {
                $n_posit=$$positform[0]."..".$$positform[1];
            } elsif ($$positform[2] eq '-') {
                $n_posit='complement('.$$positform[0]."..".$$positform[1].')';
            }
            #         print "CC $n_posit \n";
			$subsubinfo->{form}->{TYPE}->[0]= $subitem;
			$subsubinfo->{form}->{START}->[0]= $$positform[0];
			push @all_informs,$subsubinfo;
			delete $geneinfo->{info}->{$n_posit};
		}
	}

    #       print "$_\n" foreach (sort keys %{$geneinfo->{info}});
    
	push @all_informs,$geneinfo->{info}->{$_} foreach (keys %{$geneinfo->{info}});
	@all_informs_ord= sort{$a->{form}->{START}->[0] <=> $b->{form}->{START}->[0]}@all_informs;

	return [@all_informs_ord];

}
##============================================================================

sub reform_geneform {
	my ($gno,$subsubinfo,$geneidinfo)= @_;
	my ($type,$posit,$locus_tag,$gene,$dir,$start,$end,$product,$note,$fix,$cat);
	my ($ori_post,$len,$comp,@geneform);

	$type= $subsubinfo->{form}->{TYPE}->[0];
	$posit= $subsubinfo->{form}->{TOP}->[0];
	($start,$end,$dir,$len,$comp)= @{&reforme_position($posit)};

	$locus_tag= exists $subsubinfo->{form}->{locus_tag} ? $subsubinfo->{form}->{locus_tag}->[0] : $gno;
	$gene= exists $subsubinfo->{form}->{gene} ? $subsubinfo->{form}->{gene}->[0] : "";
	$product= exists $subsubinfo->{form}->{product} ? $subsubinfo->{form}->{product}->[0] : "";
	$note= exists $subsubinfo->{form}->{note} ? join " ",@{$subsubinfo->{form}->{note}} : "";
	$translation= exists $subsubinfo->{form}->{translation} ? $subsubinfo->{form}->{translation}->[0] : "";
    $translation=~ tr/ //d;
    $protein_id= exists $subsubinfo->{form}->{protein_id} ? $subsubinfo->{form}->{protein_id}->[0] : "";

	$pseudo= exists $subsubinfo->{form}->{pseudogene} ? 1 : 0;
    if ($pseudo == 1) {
        $fix= 2;
    } elsif ($comp == 0) {
        $fix= 3;
    } else {
        $fix= 1;
    }

    $cat= "A";
    
    $glen= 0;
    $glen= $type eq "CDS" ? ($end-$start-2)/3 : $end-$start+1;

	$note= $posit."\:\:".$note if ($comp == 0);
	$note= "PSEUDO\:\:".$note if ($pseudo == 1);
    
    $gene_id="";
    #   print "$gno $type $pseudo\n";
    
    if ($type eq "CDS") {
        if ($fix =~ /[23]/) {
            $gene_id= $gno;
        } elsif ($fix == 1) {
            #            $geneidinfo->{id}->{$protein_id}++;
            #            $gene_id= $geneidinfo->{id}->{$protein_id} == 1 ? $protein_id : $gno;
            #            print "Dupli_ID: $protein_id\n" if ($geneidinfo->{id}->{$protein_id} > 1);

            $geneidinfo->{id}->{$locus_tag}++;
            $gene_id= $geneidinfo->{id}->{$locus_tag} == 1 ? $locus_tag : $gno;
            print "Dupli_ID: $locus_tag\n" if ($geneidinfo->{id}->{$locus_tag} > 1);

        }
    } elsif ($type =~ /RNA/) {
        $gene_id= $locus_tag ne "" ? $locus_tag : $gno;
    } elsif ($type eq "gene" && $fix == 2) {
        $gene_id= $gno;
    }
    
    @geneform= ($type,$gene_id,$gene,$dir,$start,$end,$product,$note,$cat,$fix,"","",
         $glen,$locus_tag,$protein_id,$gno);
    
	return ([@geneform],$translation,$geneidinfo);
}

sub reforme_position {
   my ($posit)= @_;
   my ($comp,$dir,@subposits,$start,$end,$len,$sub,$subpos,@areas);

    $comp= $posit=~ /\>|\<|\&lt\;|\&gt\;|join/ ? 0 : 1;
	$dir= $posit=~ /complement\(/ ? '-' : '+';
	
	$_= $posit;
	@subposits= /(?:\<|\&lt\;)?\d+(?:\.\.(?:\>|\&gt\;)?\d+)?/g;
	($start,$end,$len)= (0,0,0);
	for ($sub=0; $sub <= $#subposits; $sub++) {
   	    $subpos= $subposits[$sub];
        
        #      print "E $subpos\n";
		$subpos=~ s/\>|\<|\&lt\;|\&gt\;//g;
        #      print "F $subpos\n";

		@areas= split/\.\./,$subpos;
		$start= $areas[0] if ($sub == 0);
		$end= $areas[-1] if ($sub == $#subposits);
		$len += $areas[-1]-$areas[0]+1;
	}
	
   return [$start,$end,$dir,$len,$comp];

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


