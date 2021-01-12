#! /usr/bin/perl -w
# revise 07/12/28

if ($ARGV[0] =~ /-h/){
  print "
    ## EXEC pair alignment
      -i MultiFasta
      -p search_prog:  [(wa(water)),bl(bl2seq),fa(fasta),ne(needle),st(stretcher))]
           -go gapopen:  [(15)]
           -ge gapextension [(6)]
           -dist  distmode [iden/p/(jc)/kim]
           ## Jukes-Cantor correction
           ## Kimura 2 parameters
      -l search_list:  All query if this option is not specified  [file name]
         # Combination  >GroupName
                          Seq1 Seq2 .....
         # Pairwise     >pair
                          seq1 seq2
      -t  [(nuc)/pep]

       \*\*\* Output format
       -om output_matrix:  
        [M:  Matrix  P:  Each group of result (H:  Only output all pair)]
         w_xx_grpmat     w_xx_grppair              w_fa_pairhom

       -o basename of output file [(w_prog_)]
       -d number of results into one file [Number (0: one file)]
         (divided result into multiple files in directory data) 
    ## EXEC output matrix
      ****  Convert identity into evolution distace only in prog needle/water
      -rf file (read results of pair alignment)
      -rd directory (read results of pair alignment)
              [-ro  objlist file]
      -rl listfile for revise matrix
              [-rm template matrix file]

      \n";
     exit;
}

$starttime= time;

use File::Basename;
$|= 1;

#==========================================================================
#  PARAMETER SETTING
#==========================================================================
foreach (@ARGV){if (/\-/){$key= $_} else {push @{$option->{$key}},$_}};
$option_str= join " ",(keys %{$option});

%hasprog= qw(wa water fa fasta bl bl2seq ne needle st stretcher al malin);

$seqtype= exists $option->{-t} ? $option->{-t}->[0] : "nuc";
$quelist= exists $option->{-l} ? $option->{-l}->[0] : "all";
$outmode= exists $option->{-om} ? $option->{-om}->[0] : "H";
$outdiv= exists $option->{-d} ? $option->{-d}->[0] : "0";

$GAPOP= exists $option->{-go} ? $option->{-go}->[0] : 15;
$GAPEXT= exists $option->{-ge} ? $option->{-ge}->[0] : 6;
$DISTMOD= exists $option->{-dist} ? $option->{-dist}->[0] : "jc";
$PAIRDEL= 1;
$LIMLEN= 50;


$HOME=$ENV{HOME};

#==========================================================================
#  MAIN 
#==========================================================================
## $seqtype $quelist $outmode $outdiv $outbase $option_str

if ($option_str=~ /\-i/) {
   $pro= exists $option->{-p} ? $option->{-p}->[0] : "bl";
   $prog= $hasprog{$pro};
   $outbase= exists $option->{-o} ? $option->{-o}->[0] : "w_${pro}";
   open (LOG,">w_${pro}_logs");

   ($seqinfo)= &readseq_fasta($option->{-i});
   print "\*\* FINISH READ SEQUENCE NUM\= $seqinfo->{taxnum}\n";
   
   ###  CREATE PAIR INFORMATION
   @seqtaxs= @{$seqinfo->{tax}};
   ($queinfo)= &create_queryinfo([@seqtaxs],$quelist);


   &output_exec_queinfo($queinfo,$prog);

   print "\n\*\* START ALIGNMENT\n";
   $pnum= @{$queinfo->{allpair}};
   print LOG "\#\# $queinfo->{type} PAIR $pnum\n";
   print LOG "  PROG\= $prog GAPOP\= $GAPOP GAPEXT\= $GAPEXT\n\#\#\n\n";
 
   ($hominfo)= &control_hom_matrix($seqinfo,$queinfo,$pro,$seqtype,$outdiv,$outbase,$outmode);

   print "\n\*\* OUTPUT RESULT\n";
   if (($queinfo->{type} eq "comb") && ($outmode eq "M")) {
      $outfile= "${outbase}_grpmat";
      open (OUT,">$outfile");
      open (OUT_M,">${outfile}.squar");
      &output_matrix($queinfo,$hominfo,$prog,$outbase);
   } elsif ($outmode eq "G") {
      $outfile= "${outbase}_grppair";
      open (OUT,">$outfile");
      &output_pair($queinfo,$hominfo,$prog,$outbase);
      &output_pair_dist($queinfo,$hominfo,$prog,$outbase,$outfile);
       
       print "OUTDIST: ${outbase}_pairhome.dist\n";
       
   } elsif ($outmode eq "H") {
      $outfile= "${outbase}_pairhom";
      &output_pair_dist($queinfo,$hominfo,$prog,$outbase,$outfile);
       
       print "OUTDIST: ${outbase}_pairhome.dist\n";
   }

   print "   EXEC PROG\= $prog SEQTYPE\= $seqtype P_TYPE\= $queinfo->{type} PAIR\= $allpairnum\n";
   print "OUTPUT\= $outfile\n";

} elsif ($option_str=~ /\-r/) {
     open (LOG,">w_mat_logs");
     ($prog,$outfile)= &main_control_r_matrix();
     print "   MATRIX PROG\= $prog OUTPUT\= $outfile\n";
}

#==========================================================================
#==========================================================================
$endtime= time;
$cuosumt= $endtime-$starttime;
($user, $system, $cuser, $csystem) = times;
print "TIMES $user $system $cuser $csystem T $cuosumt\n";
print LOG "TIMES $user $system $cuser $csystem T $cuosumt\n";
#==========================================================================
#==========================================================================
sub output_exec_queinfo {
   my ($queinfo,$prog)= @_;
   my ($grpc,$grpone,$memnum,$pairnum);

   open (PAI,">w_${pro}_pair_log");
   open (OBJL,">w_${pro}_objlist");
   print OBJL "\#\# OBJ_LIST\n";

   $grpc= 0;
   print "   NO\tGP\tMEMB\tPAIR\n";
   $grpnum= @{$queinfo->{group}};
   foreach $grpone (@{$queinfo->{group}}) {
      $grpc++;
      $memnum= $queinfo->{each}->{$grpone}->{membnum};
      $pairnum= $queinfo->{each}->{$grpone}->{cnum};
      print "   $grpc\/$grpnum\t$grpone\t$memnum\t$pairnum\n";
      print LOG ">$grpone MEMB\= $memnum PAIR\= $pairnum\n";
      print LOG "@{$queinfo->{each}->{$grpone}->{memb}}\n\n";      
      print PAI "$grpone\t$_\n"  foreach (@{$queinfo->{each}->{$grpone}->{comb}});
      print OBJL ">$grpone\n";
      print OBJL "$_\n" foreach (@{$queinfo->{each}->{$grpone}->{memb}});
   }
   $allpairnum= @{$queinfo->{allpair}};
   print "   PAIR_TOTAL $allpairnum PROG $prog\n";

}
#==========================================================================
#  CREATE QUERY INFO
#==========================================================================
sub create_queryinfo {
    my ($seqtaxs,$quelist,$listtype)= @_;
    my ($gpid,$queinfo);
    my (@taxs,@lines,$gpone,@membs,$comb,$combnum,@allpair,@uniqpair,@uniqmemb);

    $gpid= "GP";
    $queinfo= {};
#    $listtype= "comb";
    if ($quelist eq "all") {
        @taxs= @$seqtaxs;
        $queinfo->{each}->{$gpid}->{memb}= [@taxs];
        push @{$queinfo->{group}},$gpid;
        $queinfo->{type}= "comb";
    } else {
        open (LI,"$quelist");
        @lines= <LI>;
        foreach (@lines) {
            chomp;
            next if (/^(\#|\s*$)/);
            if (/^>(\S+)/){
               $gpid= $1;
               push @{$queinfo->{group}},$gpid;
               @items= split/\s+/;
               $queinfo->{each}->{$gpid}->{desc}= join " ",@items[1..$#items];
               $listtype= $gpid eq "pair" ? "pair" : "comb";
               $queinfo->{type}= $listtype;
            } else {
               push @{$queinfo->{each}->{$gpid}->{memb}},(split/\s+/);
               if ($listtype eq "pair") {
                  $pair= join " ", sort((split/\s+/)[0,1]);
                  push @{$queinfo->{allpair}},$pair;
                  push @{$queinfo->{each}->{$gpid}->{comb}},$pair;
                  $queinfo->{each}->{$gpid}->{cnum}++;
#print "$pair\n";
               }
            }
        }
    }
#    $queinfo->{type}= $listtype;
    if ($queinfo->{type} eq "pair") {
       foreach $gpone (@{$queinfo->{group}}) {
           @uniqmemb = keys %{ {map{$_,1}@{$queinfo->{each}->{$gpone}->{memb}}}};
          $queinfo->{each}->{$gpone}->{memb}= [@uniqmemb];
          $queinfo->{each}->{$gpone}->{membnum}= @uniqmemb;
       }
    } elsif ($queinfo->{type} eq "comb") {
       @allpair= ();
       foreach $gpone (@{$queinfo->{group}}) {
           @membs= @{$queinfo->{each}->{$gpone}->{memb}};
           $queinfo->{each}->{$gpone}->{membnum}= @membs;
           ($comb,$combnum)= &combination_pair_high([@membs]);
           $queinfo->{each}->{$gpone}->{comb}= $comb;
           $queinfo->{each}->{$gpone}->{cnum}= $combnum;
           push @allpair,$_ foreach (@$comb);
       }
       $queinfo->{allpair}= [@allpair];
       unless ($quelist eq "all") {
          @uniqpair = keys %{ {map{$_,1}@allpair}};
          $queinfo->{allpair}= [@uniqpair];
       }
   }

     print "TYPE: $queinfo->{type}\n";
    
     return ($queinfo);
}
#==========================================================================
#  CONTROL EXEC and OUTPUT
#==========================================================================
sub control_hom_matrix {
    my ($seqinfo,$queinfo,$pro,$seqtype,$outdiv,$outbase,$outmode)= @_;
    my ($prog,$allpairnum,@titles,$title,$hominfo,$pairnum);
    my ($pairone,$resone,$objct,$sbjct,$pairinfos,$alinfile,$pairname);
    my ($outfile,$result);

    $prog= $hasprog{$pro};
    $allpairnum= @{$queinfo->{allpair}};
    @titles= qw(OBJ O_Len SBJ S_Len Iden AlLen Mis Gap O_st O_ed S_st S_ed Eval Scor);
    push @titles,qw(DIR SEGNo SEGNUM);
    $title= join "\t",@titles;
    if ($outdiv == 0) {
       open (DOUT,">${outbase}_pairhom");
       print DOUT "\#\# PAIR $allpairnum PROG $prog\n";
       print DOUT "\#\# $title\n";
   }

   $hominfo= {};
   $pairnum= 0;
   foreach $pairone (@{$queinfo->{allpair}}) {
      $resone= {};
      $pairnum++;
      &display_counter($pairnum,[1000,10,3],"$pairnum\/$allpairnum");
       
      ($objct,$sbjct)= split/\s+/,$pairone;
      $hominfo->{each}->{$objct}->{len}= $seqinfo->{len}->{$objct};
      $hominfo->{each}->{$sbjct}->{len}= $seqinfo->{len}->{$sbjct};
      $pairinfos= [$objct,$seqinfo->{len}->{$objct},$sbjct,$seqinfo->{len}->{$sbjct}];
      ($alinfile)= &exe_alignment($objct,$sbjct,$seqinfo,$pro,$seqtype);

      $pairname= join " ",sort ($objct,$sbjct);
#print "X $pairname\n";
      if ($prog =~ /water|needle|stretcher/) {
          ($resone)= &extract_alininfo_emboss($alinfile);
      } elsif ($prog eq "bl2seq") {
          ($resone)= &extract_alininfo_bl2seq($alinfile);
      } elsif ($prog eq "fasta") {
          ($resone)= &extract_alininfo_fasta($alinfile);
      }

#      $hominfo->{homo}->{$pairname}= $resone if ($outmode=~ /(M|G)/);
       $hominfo->{homo}->{$pairname}= $resone;

      if (($outdiv > 0) && ($pairnum%$outdiv== 1)) {
          $outfile= sprintf ("outSEQ_%s_%04d",$outbase,$pairnum);
          unlink "./data/$outfile" if (-e "./data/$outfile");
          open (DOUT,">>./data/$outfile");
          print DOUT "\#\# $pairnum\/$allpairnum PROG $prog\n";
          print DOUT "\#\# $title\n";
      }
       
      ($result)= &create_onepair_result($pairinfos,$resone,$prog);
      print DOUT "$result\n";

}
   print "\n   FINISH ALIGNMNET AND EXTRACT\n";
   return ($hominfo);
}

sub create_onepair_result {
   my ($pairinfos,$resone,$prog)= @_;
   my (@resultall,$resall,$segment,$seg_tonum,$subresone,$result)= @_;

   @resultall= ();
   if (exists $resone->{0}) {
      push @resultall,join "\t",(@$pairinfos,"nohit",0);
   } else {
     ($resall)= &reform_subresult($resone,$prog);
     $segment= 0;
     $seg_tonum= @$resall;
       
#   print "DD $_\n" foreach (@$resall);
       
     foreach $subresone (@$resall) {
        $segment++;
        next if ($segment > 1);
        push @resultall,join "\t",(@$pairinfos,$subresone,$segment,$seg_tonum);
         
     }
   } 
   return (join "\n",@resultall);
}
sub reform_subresult {
   my ($resone,$prog)= @_;
   my (@item1,@item2,@resall,@res_oneline);

   @item1= qw(iden allen mis gap);
   @item2= qw(objpos sbjpos stat);

   @resall= ();
   foreach $segone (sort{$a<=>$b}(keys %{$resone})) {
     @res_oneline= ();
     push @res_oneline,$resone->{$segone}->{$_} foreach (@item1);
     push @res_oneline,@{$resone->{$segone}->{$_}} foreach (@item2);
     push @res_oneline,$resone->{$segone}->{dir};
     push @res_oneline,@{$resone->{$segone}->{subst}} if ($prog =~ /needle|water|stretcher/);
     push @resall,(join "\t",@res_oneline);
   }
   return ([@resall]);
}

#==========================================================================
#  EXECUTE ALIGNMENT PROGRAM
#==========================================================================
sub exe_alignment {
    my ($objct,$sbjct,$seqinfo,$pro,$seqtype)= @_;
    my ($objseqf,$sbjseqf,$alnfile,$cmd,$subprog,$seqtypeopt);
 
    $EMBOSS="${HOME}/biotools/local/homology/EMBOSS-6.6.0/bin";

     ($objseqf,$sbjseqf)= ("w_${pro}_objct","w_${pro}_sbjct");
     $alnfile= "w_${pro}_output";
     open (QU,">$objseqf");
     print QU "\>$objct\n$seqinfo->{seq}->{$objct}\n";
     open (SB,">$sbjseqf");
     print SB "\>$sbjct\n$seqinfo->{seq}->{$sbjct}\n";
     close(QU);
     close(SB);

print "${objct}XX\n" unless(exists $seqinfo->{seq}->{$objct});

     $prog= $hasprog{$pro};
     if ($prog eq "water") {
       $cmd= "${EMBOSS}/water $objseqf $sbjseqf -outfile $alnfile -aformat3 srspair";
       $cmd= $cmd." -gapopen $GAPOP -gapextend $GAPEXT -auto";
     } elsif ($prog eq "needle") {
       $cmd= "${EMBOSS}/needle $objseqf $sbjseqf -outfile $alnfile -aformat3 srspair";
       $cmd= $cmd." -gapopen $GAPOP -gapextend $GAPEXT -auto";
     } elsif ($prog eq "stretcher") {
       $cmd= "${EMBOSS}/stretcher $objseqf $sbjseqf -outfile $alnfile -aformat3 srspair";
       $cmd= $cmd." -gapopen $GAPOP -gapextend $GAPEXT -auto";
     } elsif ($prog eq "bl2seq") {
         $subprog= $seqtype eq "nuc" ? "blastn" : "blastp";
         $cmd= "$ENV{HOME}/biotools/local/homology/ncbi-blast/bin/$subprog -query $objseqf -subject $sbjseqf -out $alnfile";
    #     $filter = $seqtype eq "nuc" ? '-dust yes' : '-seg yes';
         $filter = $seqtype eq "nuc" ? '-dust no' : '-seg no';
         $cmd= $cmd." -evalue 10 $filter -outfmt 6";
    #    $cmd= $cmd." -e 10 -F F -D 1 -W 11";
    ## -W  Word size, default if zero (blastn 11, megablast 28, all others 3)
    ## -F  Filter query sequence (DUST with blastn, SEG with others) [String]
    ## -D  Output format: 0 - traditional, 1 - tabular [Integer];
     } elsif ($prog eq "fasta") {
       $seqtypeopt= $seqtype eq "nuc" ? "-n" : "-p";
       $cmd= "/usr/local/fasta/bin/fasta36 -H $seqtypeopt -E 10 -m 9";
     #  $cmd= "/usr/local/fasta/bin/ssearch36 -H $seqtypeopt -E 10 -m 9";
       $cmd= $cmd. " -a $objseqf $sbjseqf > $alnfile";
       ## -E	Limit the number based on the expected number of scores
       ## -H	Omit histogram
       ## -m	Specify alignment type
## print "$cmd\n";

     }
## print "EXEC($objct,$sbjct)\n    $cmd\n";
     system ($cmd);
 #    print LOG "$objct $sbjct CMD $cmd\n";

     return ($alnfile);
}
#==========================================================================
#  EXTRACT RESULT
#==========================================================================
sub extract_alininfo_bl2seq {
    my ($outfile)=@_;
    my (@line,$info,$segment);
    my ($iden,$allen,$mis,$gap,$ost,$oed,$sst,$sed,$eval,$scor);

    open (AL,"$outfile");
    @line= <AL>;
    $info= {};
    $segment= 0;
    foreach (@line) {
      chomp;
      next if (/^(\#|\s*$)/); 
      ($iden,$allen,$mis,$gap,$ost,$oed,$sst,$sed,$eval,$scor)= (split/\t/)[2..11];
      $segment++;
      $info->{$segment}->{iden}= $iden;
      $info->{$segment}->{mis}= $mis;
      $info->{$segment}->{gap}= $gap;
      $info->{$segment}->{allen}= $allen;
      $info->{$segment}->{objpos}= [$ost,$oed];
      $info->{$segment}->{sbjpos}= [$sst,$sed];
      $info->{$segment}->{stat}= [$eval,$scor];
      $info->{$segment}->{dir}= '+';
      $info->{$segment}->{dist}= ($mis+$gap)/$allen;
   }
   $info->{$segment}->{iden}= "nohit" if ($segment==0);
   $info->{$segment}->{allen}= 0 if ($segment==0);
   close (AL);
   return ($info);
}
sub extract_alininfo_fasta {
    my ($outfile)=@_;
    my (@line,$info,$segment);
    my ($head,$tail,$dir);
    my ($eval,$iden,$simi,$scor,$allen,$ost,$oed,$sst,$sed,$objgap,$sbjgap);

    open (AL,"$outfile");
    @line= <AL>;
    $info= {};
    $segment= 0;
    foreach (@line) {
      chomp;
      next if (/^(\#|\s*$)/); 
      if (/^The best scores are:/../^\>\>\>/) {
        next if (/^(The best scores|\>\>\>)/);
        $segment++;
        ($head,$tail)= split/\]\s+/;
        $dir= $head=~/\[(f|r)/ ? $1 : "f";
        $dir= $dir eq "f" ? '+' : '-';
        ($eval,$iden,$simi,$scor,$allen,$ost,$oed,$sst,$sed,$objgap,$sbjgap)
                    = (split/\s+/,$tail)[2..6,7,8,11,12,15,16];
        $info->{$segment}->{iden}= sprintf ("%.2f",$iden*100);
        $info->{$segment}->{simi}= sprintf ("%.2f",$simi*100);
        $info->{$segment}->{mis}= &ceil((1-$iden)*$allen,0);
        $info->{$segment}->{gap}= $objgap+$sbjgap;
        $info->{$segment}->{allen}= $allen;
        $info->{$segment}->{objpos}= [$ost,$oed];
        $info->{$segment}->{sbjpos}= [$sst,$sed];
        $info->{$segment}->{stat}= [$eval,$scor];
        $info->{$segment}->{dir}= $dir;
        $info->{$segment}->{dist}= sprintf ("%.4f",
              ($info->{$segment}->{mis}+$objgap+$sbjgap)/$allen);
      }
   }
   close (AL);
   return ($info);
}

sub extract_alininfo_emboss {
    my ($outfile)=@_;
    my (@line,$info,$segment,@objpos,@sbjpos,$alininfo);
    my ($query,$sbjct,$allen,$hitnum,$simnum,$gapnum,$score);
    my ($pairname,$cor_allen,$cor_gapnum,$eval);

    $EMBOSS="/home/impact/biotools/local/homology/EMBOSS-6.6.0/bin";
    
    open (AL,"$outfile");
    @line= <AL>;
    close (AL);
    $info= {};
    $segment= 0;
    @objpos= ();
    @sbjpos= ();
    $alininfo= {};
    foreach (@line) {
       chomp;
       next unless (/^\# Aligned_/../^\#\-+/);
       if (/^\# Aligned_/../^\# Score\:/) {
           $query= $1 if (/^\# 1\: (\S+)/);
           $sbjct= $1 if (/^\# 2\: (\S+)/);
           $allen= $1 if (/^\# Length\: (\d+)/);
           $hitnum= $1 if (/^\# Identity\:\s+(\d+)\/\d+\s+\(\s*([\d\.]+)%\)/);
           $simnum= $1 if (/^\# Similarity\:\s+(\d+)\/\d+\s+\(\s*([\d\.]+)%\)/);
           $gapnum= $1 if (/^\# Gaps\:\s+(\d+)\/\d+\s+\(\s*([\d\.]+)%\)/);
           $score= $1 if (/^\# Score\:\s+([\-\d\.]+)/);
           $segment++ if (/^\# Score\:\s+([\-\d\.]+)/);
        } elsif (/^\#\=+/.../^\#\-+/) {
      #     next if (/^\s*$/);
           next if (exists $alininfo->{allgap});
            if (/^\#\=+/ || /^\#\-+/) {
                $lineno= 0;
            } elsif (/^\s*$/) {
                $lineno++; $linenosub= 0;
            } elsif (/^\S+/ && $linenosub == 0) {
                $linenosub++;

#                if ($query=~ /^$linetop/) {
                unless (exists $alininfo->{alin}->{$query}) {
                    push @{$alininfo->{objs}},$query;
                    $info->{$segment}->{objpos}->[0]= (split/\s+/)[1];
                }
                $info->{$segment}->{objpos}->[1]= (split/\s+/)[3];
                $qsubseq= (split/\s+/)[2];
#   print "$qsubseq\n";
                push @{$alininfo->{alin}->{$query}},$qsubseq;
#                $qseqpos= index($_,$qsubseq,0);
                my $qseqlen= length $qsubseq;
            } elsif (/^\S+/ && $linenosub == 1) {
#               } elsif ($sbjct=~ /^$linetop/) {
                unless (exists $alininfo->{alin}->{$sbjct}) {
                    push @{$alininfo->{objs}},$sbjct;
                    $info->{$segment}->{sbjpos}->[0]= (split/\s+/)[1];
                }
                $info->{$segment}->{sbjpos}->[1]= (split/\s+/)[3];
                push @{$alininfo->{alin}->{$sbjct}},(split/\s+/)[2];
#           } else {
            } elsif (/^ +[ |.]+/) {
              #unless ((/^\#\=+/) || (/^\#\-+/)) {
              #   $mseq= substr ($_,$qseqpos,$qseqlen);
                $mseq= (split/\s+/)[1];
                push @{$alininfo->{alin}->{match}},$mseq;
#   print "$mseq\n";
#              }
           }
           if (/^#\-+/) {
              $pairname= join " ",sort($query,$sbjct);
              ($alininfo)= &correct_alinlen($alininfo,$info->{$segment});
              my ($Ts,$Tv,$Gap)= &count_substitution($alininfo);
              $cor_allen= $allen-$alininfo->{allgap};
              $cor_gapnum= $gapnum-$alininfo->{allgap};
              $info->{$segment}->{allen}= $cor_allen;
              $info->{$segment}->{iden}= sprintf ("%.2f",$hitnum/$cor_allen*100);
              $info->{$segment}->{simi}= sprintf ("%.2f",$simnum/$cor_allen*100);
              $info->{$segment}->{mis}= $cor_allen-$hitnum-$cor_gapnum;
              $info->{$segment}->{gap}= $cor_gapnum;
              $eval= 0;
              $info->{$segment}->{stat}= [$eval,$score];
              $info->{$segment}->{dist}= sprintf ("%.4f",($cor_allen-$hitnum)/$cor_allen);
              $oridist= sprintf ("%.4f",($allen-$hitnum)/$allen);
              $info->{$segment}->{dir}= '+';
              $info->{$segment}->{objpos}= $alininfo->{objarea};
              $info->{$segment}->{sbjpos}= $alininfo->{sbjarea};
              $info->{$segment}->{subst}= [$Ts,$Tv,$Gap];
              print LOG "$pairname\tALIN $allen \-\> $cor_allen  $alininfo->{tgappat}\t";
              print LOG "DIST $oridist \-\> $info->{$segment}->{dist}\n";
               
           }
        }
     }
   return ($info);
}

sub correct_alinlen {
   my ($alininfo,$info)= @_;
   my ($name,$termFgap,$termRgap,$seq);

  # $alininfo->{gaps}= [];
   $alininfo->{allgap}= 0;
   foreach $name (@{$alininfo->{objs}}) {
      ($termFgap,$termRgap)= ("","");
      $seq= join "",@{$alininfo->{alin}->{$name}};
      $termFgap= $1 if ($seq=~/^(\-+)/);
      $termRgap= $1 if ($seq=~/(\-+)$/);
      $alininfo->{gaps}->{$name}= [length($termFgap),length($termRgap)];
      $alininfo->{allgap} += length($termFgap)+length($termRgap);

   }
   ($obj,$sbj)= @{$alininfo->{objs}};
   $alininfo->{tgappat}= join "_",(@{$alininfo->{gaps}->{$obj}},@{$alininfo->{gaps}->{$sbj}});
   $alininfo->{objarea}= [$info->{objpos}->[0]+$alininfo->{gaps}->{$sbj}->[0],
                          $info->{objpos}->[1]-$alininfo->{gaps}->{$sbj}->[1]];
   $alininfo->{sbjarea}= [$info->{sbjpos}->[0]+$alininfo->{gaps}->{$obj}->[0],
                          $info->{sbjpos}->[1]-$alininfo->{gaps}->{$obj}->[1]];

   return ($alininfo);
}

sub count_substitution {
    my ($alininfo)= @_;
    my ($obj,$sbj,$objseq,$sbjseq,$matseq,$substinfo,$allenori);
    my ($estmate_st,$estmate_ed,$matstr_pos,$objstr,$sbjstr,$substpat);
    my ($Ts,$Tv,$Gap);
    ($obj,$sbj)= @{$alininfo->{objs}};
    $objseq= join "",@{$alininfo->{alin}->{$obj}};
    $sbjseq= join "",@{$alininfo->{alin}->{$sbj}};
    $matseq= join "",@{$alininfo->{alin}->{match}};
    print LOG "$obj $objseq\n$sbj $sbjseq\n";
    $allenori= length $objseq;
    $substinfo= {};
    $estmate_st= (sort{$b<=>$a}($alininfo->{gaps}->{$obj}->[0],
                                $alininfo->{gaps}->{$sbj}->[0]))[0];
    $estmate_ed= $allenori-(sort{$b<=>$a}($alininfo->{gaps}->{$obj}->[1],
                                $alininfo->{gaps}->{$sbj}->[1]))[0]-1;
    $matstr_pos= 0;
    while($matseq=~ /([ \.])/g) {
       $matstr_pos= pos $matseq;
       next if (($matstr_pos < $estmate_st) || ($matstr_pos > $estmate_ed));
       $objstr= substr($objseq,$matstr_pos-1,1);
       $sbjstr= substr($sbjseq,$matstr_pos-1,1);
       $substpat= join "",sort($objstr,$sbjstr);
       $substinfo->{$substpat}++;
    }
    $Tv= 0;
    $Ts= 0;
    $Gap= 0;
    foreach $substpat (keys %{$substinfo}) {
       $Ts += $substinfo->{$substpat} if ($substpat=~/^(AT|AC|CG|GT)$/);
       $Tv += $substinfo->{$substpat} if ($substpat=~/^(AG|CT)$/);
       $Gap += $substinfo->{$substpat} if ($substpat!~/^(AT|AC|CG|GT|AG|CT)$/);
       print LOG "$estmate_st $estmate_ed SUBST $substpat $substinfo->{$substpat}\n";
    }

    return ($Ts,$Tv,$Gap);

}
#==========================================================================
#  OUTPUT
#==========================================================================

sub output_matrix {
    my ($queinfo,$hominfo,$prog,$outbase)= @_;
    my ($gpone,$comb,$combnum,@taxs,$taxnum,$matinfo);
    my ($query,@res_one,$result,@sbjcts,$sbj,$pairname);
    my ($same,$mis,$gap,$allen,$jcdist,$pdist);
  
    foreach $gpone (@{$queinfo->{group}}) {
               ##$comb= $queinfo->{each}->{$gpone}->{comb};
               ##$combnum= $queinfo->{each}->{$gpone}->{cnum};
        @taxs= @{$queinfo->{each}->{$gpone}->{memb}};
        $taxnum= $queinfo->{each}->{$gpone}->{membnum};
        ($comb,$combnum)= &combination_pair_low([@taxs]);
        $valinfo= {};
        $same= $DISTMOD =~ /^(p|jc|kim)$/ ? (sprintf ("%.4f",0)) : 100;
        foreach (@$comb) {
            ($obj,$sbj)= split/\s+/;
            ($objlen,$sbjlen)=($hominfo->{each}->{$obj}->{len},$hominfo->{each}->{$sbj}->{len});
            $pairname= join " ",sort($obj,$sbj);
            if (exists $hominfo->{homo}->{$pairname}) {
               if (exists $hominfo->{homo}->{$pairname}->{0}) {
                $dist = $DISTMOD =~ /^(p|jc|kim)$/ ? '2.0000' : ' -- ';
               } else {
                 $mis= $hominfo->{homo}->{$pairname}->{1}->{mis};
                 $gap= $hominfo->{homo}->{$pairname}->{1}->{gap};
                 $allen= $hominfo->{homo}->{$pairname}->{1}->{allen};
                 $subst= $hominfo->{homo}->{$pairname}->{1}->{subst};
                 $objcov= $allen/$objlen;
                 $sbjcov= $allen/$sbjlen;
                 if ($DISTMOD eq "iden") {
                    $dist= $hominfo->{homo}->{$pairname}->{1}->{iden};
                    if ($objcov < 0.8 && $sbjcov < 0.8) {
                        $dist= "${dist}*";
                    }
                 } elsif ($DISTMOD =~ /^(p|jc)$/) {
                    ($jcdist,$pdist)= &estimate_jc_dist($allen,$mis,$gap);
                    $dist= $DISTMOD eq "p" ? $pdist : $jcdist;
                    $dist= '2.0000' if ($allen < $LIMLEN);
                 } elsif ($DISTMOD =~ /^(kim)$/) {
                    $dist= &estimate_kim_dist($allen,$subst);
                    $dist= '2.0000' if ($allen < $LIMLEN);
                 }
               }
            } else {
               $dist = $DISTMOD =~ /^(p|jc|kim)$/ ? '2.0000' : ' -- ';
            }
            $valinfo->{each}->{$obj}->{$sbj}= $dist;
            push @{$valinfo->{all}->{$obj}},$sbj;
        }

      print OUT "\# $gpone  PROG $prog  DIST $DISTMOD  PAIRDEL $PAIRDEL\n";
      print OUT "\# ${taxnum}_$combnum\t".(join "\t",@taxs)."\n";
      print OUT "  ${taxnum}\n";
      
      foreach $raw (@taxs) {
         @res_one= ();
         if (exists $valinfo->{all}->{$raw}) {
            foreach $col (@{$valinfo->{all}->{$raw}}) {
               if (exists $valinfo->{each}->{$raw}->{$col}) {
                 push @res_one,$valinfo->{each}->{$raw}->{$col};
               } else {
                 print "ERR  $raw $col\n";
               }
            }
         }
         push @res_one,$same;
#         $resone_txt= join " ",@res_one;
#         printf OUT ("%-10s  %s\n",$raw,$resone_txt);
         print OUT (join "\t",($raw,@res_one))."\n";
print "$raw\n";
      }

###  COMPLETE MATRIX

      foreach $raw (@taxs) {
         @res_one= ();
 #        if (exists $valinfo->{all}->{$raw}) {
            foreach $col (@taxs) {
					if ($raw ne $col) {
               	if (exists $valinfo->{each}->{$raw}->{$col}) {
                 		push @res_one,$valinfo->{each}->{$raw}->{$col};
               	} elsif (exists $valinfo->{each}->{$col}->{$raw}) {
                 		push @res_one,$valinfo->{each}->{$col}->{$raw};
						} else {
  	               	print "ERR  $raw $col\n";
               	}
					} else {
						push @res_one,100;

					}
            }
#         }

         print OUT_M (join "\t",($raw,@res_one))."\n";
## print "RAW $raw\n";
      }
 
   }

}

               
sub output_pair_dist {
    my ($queinfo,$hominfo,$prog,$outbase,$outfile2)= @_;
    my ($gpone,$comb,$combnum,@taxs,$taxnum,$matinfo);
    my ($query,@res_one,$result,@sbjcts,$sbj,$pairname);
    my ($same,$mis,$gap,$allen,$jcdist,$pdist);
        
    open (OUT2,">${outfile2}.dist");
               
    foreach $gpone (@{$queinfo->{group}}) {
        $comb= $queinfo->{each}->{$gpone}->{comb};
        $valinfo= {};
        $same= $DISTMOD =~ /^(p|jc|kim)$/ ? (sprintf ("%.4f",0)) : 100;
        foreach (@$comb) {
            ($obj,$sbj)= split/\s+/;
            ($objlen,$sbjlen)=($hominfo->{each}->{$obj}->{len},$hominfo->{each}->{$sbj}->{len});
            $pairname= join " ",sort($obj,$sbj);
 
            if (exists $hominfo->{homo}->{$pairname}) {
#print "X $pairname\n";
                if (exists $hominfo->{homo}->{$pairname}->{0}) {
                    $dist = $DISTMOD =~ /^(p|jc|kim)$/ ? '2.0000' : ' -- ';
                } else {
                    $mis= $hominfo->{homo}->{$pairname}->{1}->{mis};
                    $gap= $hominfo->{homo}->{$pairname}->{1}->{gap};
                    $allen= $hominfo->{homo}->{$pairname}->{1}->{allen};
                    $subst= $hominfo->{homo}->{$pairname}->{1}->{subst};
                    $objcov= $allen/$objlen;
                    $sbjcov= $allen/$sbjlen;
                    if ($DISTMOD eq "iden") {
                        $dist= $hominfo->{homo}->{$pairname}->{1}->{iden};
                        if ($objcov < 0.8 && $sbjcov < 0.8) {
                            $dist= "${dist}*";
                        }
                    } elsif ($DISTMOD =~ /^(p|jc)$/) {
                        ($jcdist,$pdist)= &estimate_jc_dist($allen,$mis,$gap);
                        $dist= $DISTMOD eq "p" ? $pdist : $jcdist;
                        $dist= '2.0000' if ($allen < $LIMLEN);
                    } elsif ($DISTMOD =~ /^(kim)$/) {
                        $dist= &estimate_kim_dist($allen,$subst);
                        $dist= '2.0000' if ($allen < $LIMLEN);
                    }
                }
            } else {
#print "Y $pairname\n";
               $dist = $DISTMOD =~ /^(p|jc|kim)$/ ? '2.0000' : ' -- ';
            }
            $valinfo->{each}->{$obj}->{$sbj}= $dist;
            push @{$valinfo->{all}->{$obj}},$sbj;
               
               print OUT2 "$obj\t$sbj\t$dist\n";
               
        }
                     
    }

}
               
sub output_pair {
   my ($queinfo,$hominfo,$prog)= @_;
   my (@titles,$gpone,$comb,$combone,$obj,$sbj,$pairinfos,$pairname,$result);
   
    @titles= qw(OBJ O_Len SBJ S_Len Iden AlLen Mis Gap O_st O_ed S_st S_ed Eval Scor);

   foreach $gpone (@{$queinfo->{group}}) {
      $comb= $queinfo->{each}->{$gpone}->{comb};
      print OUT "\#\# GP\= $gpone PROG\= $prog $queinfo->{each}->{$gpone}->{desc}\n";
      print OUT "\#\# ".(join "\t",(@titles,qw(DIR SEG))),"\n";
      foreach $combone (@$comb) {
         ($obj,$sbj)= split/\s+/,$combone;
         $pairinfos= [$obj,$hominfo->{each}->{$obj}->{len},
                      $sbj,$hominfo->{each}->{$sbj}->{len}];
         $pairname= join " ",sort($obj,$sbj);
         ($result)= &create_onepair_result($pairinfos,$hominfo->{homo}->{$pairname});
         print OUT "$result\n";

      }
   }

}


#==========================================================================
sub readseq_fasta {
  my ($seqfile)= @_;
  my ($seqinfo,$seqname,$fileone,@lines);
  
  $seqinfo={};
  foreach $fileone (@$seqfile) {
     open (SEQ,"$fileone");
     @lines= <SEQ>;

     foreach(@lines){
       chomp;
       if (/^>(\S+)/){
         $seqname=$1;
         push @{$seqinfo->{tax}},$seqname;
         $seqinfo->{taxnum}++;
       }else{
         push @{$seqinfo->{seqline}->{$seqname}},$_;
       }
     }
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
     $sequence=~ tr/U/T/;
     $sequence=~ tr/a-z/A-Z/;
     $seqinfo->{seq}->{$seqname}= $sequence;
     $seqinfo->{len}->{$seqname}= length $sequence;
     delete $seqinfo->{seqline}->{$seqname};
#    print "$seqname $seqinfo->{len}->{$seqname} $sequence\n";
   #           print "$seqname $sequence\n";

               
  }
  return ($seqinfo);
}

##=============================  SEQ LINE =============================
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
##===============================================================================
#  CREATE MATRIX
##===============================================================================
sub combination_pair_low {
   my ($grps)= @_;
   my (@comb,@agrp,@bgrp,@subjs,$obj,$sbj);
   
   @comb= ();
   @agrp= @$grps;
   @bgrp= @$grps;
   shift @agrp;
   @subjs= ();
   while (@agrp) {
     $obj= shift @agrp;
     push @subjs,shift @bgrp;
     foreach $sbj (@subjs) {
        ## push @comb,[$obj,$sbj];
        push @comb,(join " ",($obj,$sbj));
     }
   }
   $combnum= @comb;
   return ([@comb],$combnum);
}
sub combination_pair_high {
   my ($grps)= @_;
   my (@group,@comb,$combnum,$aone,$bone,$pairname);
   
   @group= @$grps;
   $combnum= 0;
   @comb= ();
   open (LI,">w_pairlist0");
   while (@group > 1) {
      $aone= shift @group;
      foreach $bone (@group) {
         $combnum++;
     #    push @comb,[$aone,$bone];
         $pairname= join " ",($aone,$bone); 
         push @comb,$pairname;
         &display_counter($combnum,[5000,10,3],"$combnum");
         if ($combnum%10000== 1) {
           ##open (LI,">w_pairlist$combnum");
           print LI ">pair\n";
         }
         print LI "$pairname\n";
      
      }
   }
   print " TOTAL $combnum\n";
   return ([@comb],$combnum);
}

##===============================================================================
#  CREATE MATRIX
##===============================================================================
##  UNIQ
##  @uniq = sort keys %{ {map{$_,1}@list}};

sub create_matrix_low {
  my ($array,$k,$mode)=@_;
  my (@revs,$comb,$combnum,$i,%hasord,@combin,@comninate);

#  @revs= reverse (@revs);
   @revs= @$array;
  ($comb,$combnum)= &permute_combinate([@revs],2,"comb");
  $i=0;
  foreach (@$array) {
     $hasord{$_}= $i;  $i++;
  }
  @combin= map{my @x = @$_;
           [$hasord{$x[0]} > $hasord{$x[1]} ? @x[0,1] : @x[1,0]];
           }@$comb;

  @comninate= sort {
       $hasord{$a->[0]} <=> $hasord{$b->[0]} or
       $hasord{$a->[1]} <=> $hasord{$b->[1]}} @combin;
  return ([@comninate],$combnum);
}

sub create_matrix_high {
  my ($array,$k,$mode)=@_;
  my (@revs,$comb,$combnum,$i,%hasord,@combin,@comninate);

#  @revs= reverse (@revs);
   @revs= @$array;
  ($comb,$combnum)= &permute_combinate([@revs],2,"comb");
  $i=0;
  foreach (@$array) {
     $hasord{$_}= $i;  $i++;
  }
  @combin= map{my @x = @$_;
           [$hasord{$x[0]} < $hasord{$x[1]} ? @x[0,1] : @x[1,0]];
           }@$comb;

  @comninate= sort {
       $hasord{$a->[0]} <=> $hasord{$b->[0]} or
       $hasord{$a->[1]} <=> $hasord{$b->[1]}} @combin;
  return ([@comninate],$combnum);
}


sub permute_combinate {
  my ($array,$k,$mode)=@_;
  my (@result,$numnum,$num,@combinates,$combnum);
  my (@allperms,$one,@permutes,$permnums);

  # use Math::Combinatorics;
  
  $num= @$array; 
  @result=();
  $numnum=""; 
  if ($mode eq "comb") {
     @combinates = combine($k,@$array);
     $combnum= @combinates;
  #  $combnum= int(factorial($num)/factorial($k)/factorial($num-$k));
     @result= @combinates;
     $numnum= @combinates;
  } elsif ($mode eq "perm") {
     @combinates = combine($k,@$array);
     @allperms=();
     foreach $one (@combinates) {
       @permutes = permute(@$one);
       push @allperms,@permutes;
     }
     @result= @allperms;
     $permnums= @result;
  #  $permnums= int(factorial($num)/factorial($num-$k));
     $numnum= $permnums;
  }
  return ([@result],$numnum);

#   my $c = Math::Combinatorics->new(count => 2, data => [@arrays]);
#   while (my @comb = $c->next_combination){
#     print join (' ',@comb)."\n";
#   }
  
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
#==========================================================================
#=============================  CEIL and TRUNC =============================
sub ceil {
	my $val = shift;	
	my $col = shift;	
	my $r = 10 ** $col;
	if ($val > 0) {
		my $tmp = $val * $r;
		if ($tmp == int($tmp)) {
			return $tmp / $r;
		} else {
			return int($tmp + 1) / $r;
		}
	} else {
		return int($val * $r) / $r;
	}
}

sub trunc {
	my $val = shift;
	my $col = shift;
	my $r = 10 ** $col;
	if ($val > 0) {
		return int($val * $r) / $r;
	} else {
		my $tmp = $val * $r;
		if ($tmp == int($tmp)) {
			return $tmp / $r;
		} else {
			return int($tmp - 1) / $r;
		}
	}
}
#==========================================================================
#==========================================================================
#  CREATE DIST MATRIX FROM HOMOPAIR OR TEMP MATRIX
#==========================================================================
#==========================================================================
sub main_control_r_matrix {
   my ($inopt,$filepaths,$distmat,$r_objinfo,$pro,$outfile,$outfiles,$outfile_r);
   
   $outfiles= "";
   if ((exists $option->{-rf}) || (exists $option->{-rd})) {
      ###=====   CREATE INFILE LIST  ============
      $inopt= exists $option->{-rf} ? '-rf' : '-rd';
      ($filepaths)= &collect_infiles($inopt,"outSEQ");
      ($distmat,$r_objinfo)= &control_create_r_matrix($filepaths);
      $pro= substr($r_objinfo->{prog},0,2);
      $outfile= "w_${pro}r_grpmat";
      open (OUT,">$outfile");
      print OUT "  $r_objinfo->{objnum}  $r_objinfo->{pnum}\n$distmat\n";
      $outfiles="$outfile ";
   } elsif (exists $option->{-rm}) {
      ($distmat,$r_objinfo)= &read_r_matrix($option->{-rm}->[0]);
      $pro= substr($r_objinfo->{prog},0,2);
   }

   ###=====   REVISE MATRIX  ============
   if (exists $option->{-rl}) {
        $outfile_r= "w_${pro}r_grpmat_rev";
        print "\n  \#\# REVISE HOMO PAIR BY LIST $option->{-rl}->[0]\n";
        &control_revise_matrix($option->{-rl}->[0],$r_objinfo,$distmat,$outfile_r);       
        $outfiles= $outfiles.$outfile_r;
   }
   
   return ($r_objinfo->{prog},$outfiles);
}

#==========================================================================
#   COLLECTION FILES
#==========================================================================
sub collect_infiles {
   my ($inopt,$selform)= @_;
   my ($inmode,$inpath,$base,$path,@filepaths,$indir,@datanames);

   $inmode= $inopt=~ /\-\wd/ ? "directory" : "file";
   $inpath= $option->{$inopt}->[0];

   @filepaths= ();
   if ($inmode eq "file"){
      ($base,$path)= fileparse($inpath);
      @filepaths=("${path}$base");
   } elsif ($inmode eq "directory"){
      $indir= $inpath;
      opendir (QD,"$indir");
      @datanames= readdir QD;
      ($base,$path)= fileparse($inpath);
      @filepaths= 
           map{$path.$base.'/'.$_}
           grep {/$selform/}@datanames;
   }

   return ([@filepaths]);
}
#==========================================================================
#   CREATE MATRIX FROM PAIR HOMO
#==========================================================================
sub control_create_r_matrix {
    my ($filepaths)= @_;
    my ($r_objinfo,@objorder,$tempmat,$combnum,$distmat,$prog,$diststac);
 
      ###=====   CREATE OBJ OINFO============
      print "  \n\#\#  CREATE OBJLIST\n";
    if (exists $option->{-ro}) {
        ($r_objinfo)= &create_r_objlist($option->{-ro}->[0]);
    } else {
        ($r_objinfo)= &extract_obj_pairfile($filepaths);
    }
      ###=====   CREATE TEMP MATRIX  ============
      print "  \n\#\#  CREATE MATRIX TABLE  OBJ\= $r_objinfo->{objnum}\n";
    @objorder= 0..$r_objinfo->{objnum}-1;
    ($tempmat,$combnum)= &create_matrix_lowtable($r_objinfo->{obj},[@objorder]);
    $r_objinfo->{pnum}= $combnum;

     ###=====   CREATE CONVERT DATA  ============
    print "  \n\#\# READ PAIRHOMO FILE AND CREATE MATRIX  PAIR\= $combnum\n";
    ($distmat,$prog,$dstatic)= 
            &convert_matrix_dist($filepaths,$tempmat,$r_objinfo->{no});
    $r_objinfo->{dnum}= $$dstatic[0];
    $r_objinfo->{ncnum}= $$dstatic[1];
    $r_objinfo->{prog}= $prog;
    print "  \n\#\# OUTPUT HOMO MATRIX\n";

    print LOG "\n\#\# MAIN MATRIX FROM DATA\n";
    print LOG "   OBJ\= $r_objinfo->{objnum} PAIR\= $r_objinfo->{pnum} ";
    print LOG "PRO\= $r_objinfo->{prog}\n";
    print LOG "   DIST\= $DISTMOD DISTNO\= $$dstatic[0] DISTNC\= $$dstatic[1]  ";
    print LOG "MIN-MAX-AVE @$dstatic[5..7]\n";
    print LOG "   ALLEN\n";
    print LOG "      \< L60\= $$dstatic[2]\n      L60\-80\= $$dstatic[3]\n";
    print LOG "      \> L80\= $$dstatic[4]\n\n";

    return ($distmat,$r_objinfo);
}
#==========================================================================
sub create_r_objlist {
   my ($fileone)= @_;
   my ($r_objinfo,@list);
   $r_objinfo= {};
   open (OBJLIST,"$fileone");
   @list= <OBJLIST>;

   foreach (@list) {
      chomp;
      next if (/^(\#|\s*$)/);
      next if (/^\>/);
      push @{$r_objinfo->{obj}},split/\s+/;
   }

   $r_objinfo->{objnum}= @{$r_objinfo->{obj}};

   for ($objno= 0; $objno < $r_objinfo->{objnum}; $objno++) {
     $r_objinfo->{no}->{$r_objinfo->{obj}->[$objno]}= $objno;
   }
   return ($r_objinfo);
}
sub extract_obj_pairfile {
   my ($files)= @_;
   my ($r_objinfo,$fileone,@list,$obj,$objlen,$sbj,$sbjlen);
   
   $r_objinfo= {};

   foreach $fileone (@$files) {
      print "$fileone\n";
      open (OBJLIST,"$fileone");
      @list= <OBJLIST>;
      foreach (@list) {
         chomp;
         next if (/^(\#|\s*$)/);
         next if (/^\>/);
         ($obj,$objlen,$sbj,$sbjlen)= /^(\S+)\t(\d+)\t(\S+)\t(\d+)\t/;
         foreach ([$obj,$objlen],[$sbj,$sbjlen]) {
            push @{$r_objinfo->{obj}},$$_[0] unless (exists $r_objinfo->{each}->{$$_[0]});
            $r_objinfo->{each}->{$$_[0]}= $$_[1];
         }
      }
   }
   $r_objinfo->{objnum}= @{$r_objinfo->{obj}};

   for ($objno= 0; $objno < $r_objinfo->{objnum}; $objno++) {
     $r_objinfo->{no}->{$r_objinfo->{obj}->[$objno]}= $objno;
   }
   return ($r_objinfo);
}
#==========================================================================
#==========================================================================
sub convert_matrix_dist {
   my ($files,$matrix,$objnos)= @_;
   my ($distno,$distnc,@matrixlin,@dist_stac,$fileone,$r_homoinfo,$index,$lino,$value);
   my ($obj,$noobj,$dist,$objname,$distmat);

   $distno= 0;
   $distnc= 0;
   @matrixlin= split/\n/,$matrix;
   @dist_stac= ();
   my $alleninfo= {};
   $alleninfo->{mindist}= 2;
   $alleninfo->{maxdist}= 0;
   foreach $fileone (@$files) {
      ($r_homoinfo)= &extract_homo_pairfile($fileone,$objnos);
      foreach $index (keys %{$r_homoinfo->{index}}) {
         $lino= $1 if ($index=~ /(\d+)\-(\d+)/);
         $value= $r_homoinfo->{index}->{$index};

         $alleninfo->{allen}->{$r_homoinfo->{allen}->{$index}}++;
         $alleninfo->{mindist}= $value if ($value <= $alleninfo->{mindist});
         $alleninfo->{maxdist}= $value if ($value >= $alleninfo->{maxdist});
         $alleninfo->{sumdist} += $value;
         $alleninfo->{numb}++;

         $distno++ if ($value < 2);
         $distnc++ if ($value == 2);
         $matrixlin[$lino]=~ s/ $index / $value /;

         &display_counter($distno,[10000,10,3],"E $distno/$index");
      }
      print "  $fileone $distno\n";
   }
 
   foreach $obj (keys %{$objnos}) {
      $noobj= $objnos->{$obj};
      $dist= sprintf ("%.4f",0);
      $objname= sprintf ("%-10s",$obj);
      $matrixlin[$noobj]=~ s/\b$noobj\-$noobj$/$dist/;
      $matrixlin[$noobj]=~ s/^$noobj\b/$objname/;
   }
   
   ##  Sum Short ALlen
   my ($low_al6,$low_al8,$low_al9)= (0,0,0);

   foreach (sort{$a<=>$b}(keys %{$alleninfo->{allen}})) {
      $low_al6 += $alleninfo->{allen}->{$_} if ($_ <= 60);
      $low_al8 += $alleninfo->{allen}->{$_} if (($_ > 60) && ($_ <= 80));
      $low_al9 += $alleninfo->{allen}->{$_} if ($_ > 80);
   }
   my @short_al= ($low_al6,$low_al8,$low_al9);
   my $avedist= sprintf ("%.4f",$alleninfo->{sumdist}/$alleninfo->{numb});
   my $dstatic= [$distno,$distnc,$low_al6,$low_al8,$low_al9,
                 $alleninfo->{mindist},$alleninfo->{maxdist},$avedist];

   $distmat= join "\n",@matrixlin;

   return ($distmat,$r_homoinfo->{prog},$dstatic);
}

sub extract_homo_pairfile {
   my ($fileone,$objnos)= @_;
   my ($pairno,@list,$r_homoinfo,@resone,$obj,$objlen,$sbj,$sbjlen);
   my ($iden,$allen,$mis,$gap,$objno,$sbjno,$diff,$dist,$pairone);

   open (PAIRHOM,"$fileone");
#   @list= <OBJLIST>;
   $r_homoinfo= {};
   $pairno= 0;
   while (<PAIRHOM>) {
       chomp;
       $r_homoinfo->{prog}= $1 if (/\#\#.+PROG (\w+)/);
       next if (/^(\#|\s*$)/);
       $pairno++;
       @resone= split/\t/;
       ($obj,$objlen,$sbj,$sbjlen,$iden,$allen)= @resone[0..5];
       next unless ((exists $objnos->{$obj}) && (exists $objnos->{$sbj}));
       $objno= $objnos->{$obj};
       $sbjno= $objnos->{$sbj};

       if ($iden ne "nohit") {
          ($mis,$gap,$ts,$tv,$tgap)= @resone[6..7,15..17];
       } else {
          ($mis,$gap,$ts,$tv,$tgap)= ("","","","","");
       }

      if ($allen < $LIMLEN) {
         printf LOG ("%-10s %-10s \(%3d %3d\)  ALLEN\= %3d GAP\= %3s\n",
                        $obj,$sbj,$objlen,$sbjlen,$allen,$gap);
      }
      if ($iden ne "nohit") {
           if ($DISTMOD eq "iden") {
              $dist= $iden;
           } elsif ($DISTMOD =~ /^(p|jc)$/) {
              ($jcdist,$pdist)= &estimate_jc_dist($allen,$mis,$gap);
              $dist= $DISTMOD eq "p" ? $pdist : $jcdist;
              $dist= '2.0000' if ($allen < $LIMLEN);
           } elsif ($DISTMOD =~ /^(kim)$/) {
              $dist= &estimate_kim_dist($allen,[$ts,$tv,$tgap]);
              $dist= '2.0000' if ($allen < $LIMLEN);
           }
      } else {
         $dist= '2.0000';
         $allen= $objlen > $sbjlen ? $objlen : $sbjlen;
      }
     
      $pairone= join '-',sort{$b<=>$a}($objno,$sbjno);
      $r_homoinfo->{index}->{$pairone}= $dist;
      $r_homoinfo->{allen}->{$pairone}= $allen;
    
      &display_counter($pairno,[10000,10,3],"$pairno\/$fileone");
   }
   print " SUMP $pairno\n";
   return ($r_homoinfo);

}

#==========================================================================
#   CREATE TEMP MATRIX
#==========================================================================
sub create_matrix_lowtable {
   my ($membs,$orders)= @_;
   my (@memb,@order,$combnum,@allpair,$apo,$bpo,$pairname);
   
   @memb= @$membs;
   @order= @$orders;

   $combnum= 0;
   @allpairs= ();
   $pos= 0;
   foreach $apo (@order) {
      $bpo= $apo;
      @subpairs= ($apo);
      foreach $bpo (@order[0..$pos]) {
         $combnum++;
         $pair= join "-",sort{$b<=>$a}($apo,$bpo);
         push @subpairs,$pair;
      }
      push @allpairs,(join " ",@subpairs);
      $pos++;      
   }
   $matrix= join "\n",@allpairs;

##   $matrix=~ s/(\d+)\-(\d+)/$memb[$1]\-$memb[$2]/g;
##   print "$matrix\n";   
   return ($matrix,$combnum);
}
#==========================================================================
#   REVISE MATRIX 
#==========================================================================
sub read_r_matrix {
   my ($rmatfile)= @_;
   my ($r_objinfo,@lines,$top,$objnum,$pnum,$rawno,$obj,$distmat);

   open (RMAT,"$rmatfile");
   @lines= <RMAT>;
   $r_objinfo= {};
   $top= shift @lines;
   ($objnum,$pnum)= ($1,$2) if ($top=~ /(\d+)\s+(\d+)/);
   $r_objinfo->{objnum}= $objnum;
   $r_objinfo->{pnum}= $pnum;
   $r_objinfo->{prog}= "xxx";
   $rawno= 0;
   foreach (@lines) {
      chomp;
      $obj= $1 if (/^(\S+)\b/);
      push @{$r_objinfo->{obj}},$obj;
      $r_objinfo->{each}->{$obj}= 100;
      $r_objinfo->{no}->{$obj}= $rawno;
      $rawno++;
   }
   $distmat= join "\n",@lines;

   print LOG "\#\# MAIN MATRIX FROM MAT\n";
   print LOG "   OBJ\= $objnum PAIR\= $pnum PRO\= $r_objinfo->{prog}\n\n";

   return ($distmat,$r_objinfo);
}

#==========================================================================
#==========================================================================
sub control_revise_matrix {
   my ($revfile,$r_objinfo,$orimat,$outfile)= @_;
   my ($rg_objinfo,@orimatrixs,$grp,$grpmemb,$grporder,$rev_matrix,$grp_combnum);
   my (@rev_matrixs,$rawone,@pairs,$rawone_index,$rawone_name);
   my ($pairone,$rawno,$colno,$orival,$rev_distmat);
   my (@dist_stac);
  
   open (OUTR,">$outfile");

   ($rg_objinfo)= &create_rg_objinfo($revfile,$r_objinfo);
   
   print LOG "\#\# GRP REV MATRIX\n";

   @orimatrixs= split/\n/,$orimat;

   foreach $grp (sort keys %{$rg_objinfo}) {
      $grpmemb= $rg_objinfo->{$grp}->{obj};
      $grporder= $rg_objinfo->{$grp}->{order};
      ($rev_matrix,$grp_combnum)= &create_matrix_lowtable($grpmemb,$grporder);
      @rev_matrixs= split/\n/,$rev_matrix;
      @dist_stac= ();
      $distno= 0;
      for($rawone=0; $rawone <= $#rev_matrixs; $rawone++) {
          $_= $rev_matrixs[$rawone];
          @pairs= /\b(\d+\-\d+)\b/g;
          $rawone_index= $1 if (/^(\d+)\b/);
          $rawone_name= sprintf ("%-10s",$rg_objinfo->{$grp}->{obj}->[$rawone]);
          foreach $pairone (@pairs) {
             ($rawno,$colno)= ($1,$2) if ($pairone=~ /(\d+)\-(\d+)/);
             $orival= (split/\s+/,$orimatrixs[$rawno])[$colno+1];
             ## $orival= &estimate_jc_dist($orival);

             $rev_matrixs[$rawone]=~ s/\b$pairone\b/$orival/;

             if ($rawno ne $colno) {
                @dist_stac= (sort{$a<=>$b}(@dist_stac,$orival))[0,-1];
                $distno++ if ($orival != 1);
             }
          }
          $rev_matrixs[$rawone]=~ s/^$rawone_index\b/$rawone_name/;
      }

      $rev_distmat= join "\n",@rev_matrixs;
      print OUTR "\#\# $grp\n";
      print OUTR "  $rg_objinfo->{$grp}->{objnum}  $grp_combnum\n";
      print OUTR "$rev_distmat\n\n";

      print LOG "   GRP $grp\n";
      print LOG "   OBJ\= $rg_objinfo->{$grp}->{objnum} PAIR\= $grp_combnum ";
      print LOG "PRO\= $r_objinfo->{prog}\n";
      print LOG "   DISTNO\= $distno MIN-MAX @dist_stac\n\n";
   }

}

sub create_rg_objinfo {
   my ($fileone,$r_objinfo)= @_;
   my ($rg_objinfo,@list,$grp,$obj,@sortord);

   $rg_objinfo= {};
   open (OBJLIST,"$fileone");
   @list= <OBJLIST>;
   foreach (@list) {
      chomp;
      next if (/^(\#|\s*$)/);
      if (/^\>(\S+)/) {
         $grp= $1;
      } else {
         push @{$rg_objinfo->{$grp}->{obj}},split/\s+/;
      }
   }
   
   foreach $grp (keys %{$rg_objinfo}) {
      my @grpmembs= @{$rg_objinfo->{$grp}->{obj}};
      $rg_objinfo->{$grp}->{obj}= [];
      foreach $obj (@grpmembs) {
         if (exists $r_objinfo->{no}->{$obj}) {
            push @{$rg_objinfo->{$grp}->{obj}},$obj;
            push @{$rg_objinfo->{$grp}->{order}},$r_objinfo->{no}->{$obj};
         } else {
            print LOG  "\#\# ERR_OBJ $obj\n";
         }
      }
      $rg_objinfo->{$grp}->{objnum}= @{$rg_objinfo->{$grp}->{obj}};

      ## @sortord= sort{$a<=>$b}(@{$r_objinfo->{$grp}->{order}});
      ## $r_objinfo->{$grp}->{order}= [@sortord];
   }
   return ($rg_objinfo);
}
#==========================================================================
#==========================================================================
sub estimate_jc_dist {
  my ($allen,$mis,$gap)= @_;
  my ($diff,$pvalnorm,$pvaldel,$pdist,$jcdist);

  $diff= $mis+$gap;

  $pdistnorm= ($mis+$gap) == 0 ? 0 : ($mis+$gap)/$allen;

  $pdistdel= $mis == 0 ? 0 : $mis/($allen-$gap);

  $pdist= $PAIRDEL == 1 ? $pdistdel : $pdistnorm;

  if ($pdist < 0.75) {
     $jcdist= 0-3/4*log(1-(4/3*$pdist));
     $jcdist= sprintf ("%.4f",$jcdist);
  } else {
     $jcdist= 2;
  }
  $jcdist= sprintf ("%.4f",$jcdist);
  $pdist= sprintf ("%.4f",$pdist);

  return ($jcdist,$pdist);

}
sub estimate_kim_dist {
  my ($allen,$subst)= @_;
  my ($al_len,$ts,$tv,$gap,$p_ts,$p_tv,$kimdist);
  
  ($ts,$tv,$gap)= @$subst;
  $al_len= $PAIRDEL == 1 ? $allen-$gap : $allen;

  $p_ts= $ts > 0 ? $ts/$al_len : 0;  ## P
  $p_tv= $tv > 0 ? $tv/$al_len : 0;  ## Q

  my $w1= 1-2*$p_ts-$p_tv;
  my $w2= 1-2*$p_tv;

  if (($w1 > 0) && ($w2 > 0)) {
     $kimdist= 0-(1/2*log($w1))-(1/4*log($w2));
  } else {
     $kimdist= 2;
  }
  $kimdist= sprintf ("%.4f",$kimdist);

  return $kimdist;

}

