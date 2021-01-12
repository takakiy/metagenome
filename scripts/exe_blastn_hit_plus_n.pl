#! /usr/bin/perl -w
#



if ($ARGV[0] =~ /-h/){
   print "EXEC OPTION;
            exec blastall; -e query -d database
            exec make_hitlist; -r blastout file_path
            exec make_hitlist_folder; -rd blastout_file_path
             \(exp. ~/folder/.*out\)\n";
   print "FILTER OPTION;           
            -n max_number(5)               ** 0 = do not limited
                   
            -s remove same sequence (T)/F
            
            -f filttering (T)/F
               -cq cover_rate for query % (0) ** 0 = do not limited
               -cs cover_rate for sbjct % (0) ** 0 = do not limited
               -lh homology% (80)              ** 0 = do not limited
               -lp Probability (1e-5)          ** 10 = do not limited
            
            -l length (100)                ** 0 = do not limited
           \n";
                
   print "\nOPTION;   -org T/F
                 -o output \n";
  
   exit;
}elsif ($ARGV[0] !~ /\-[er]/){
   print "ERROR;    Retry in adding option (-e -r -rd)\n";
  # exit;
}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

#================================================================================
#      CONTROL EXECUTIVE MODE 
#================================================================================
$blastmode= "blastn";
use File::Basename;

	$blast_home="/home/bioresource/biotools/local/homology/blast";
	$blast_home_plus="/home/bioresource/biotools/local/homology/ncbi-blast";
	$show_gis= "";

@filepaths=();
if (exists $option->{-e}){
	$infile= $option->{-e}->[0];
	$basename= (split/\//,$infile)[-1]; $head= (split/\./,$basename)[0];
	$blastout_file="temp_blastnout".$head;

	$algn_num= exists $option->{-algn_num} ? $option->{-algn_num}->[0] : 1000;

	if (exists $option->{-d}){
		$database=$option->{-d}->[0];
	}else{
		print "Database ?\n";
		$database=<STDIN>;
		chomp $database;
	}
	if ($database =~ /^(nr|kegg|nt)$/) {
		$database= "${blast_home}/db/$database";
	}
	if (exists $option->{-show_gis}) {
		if ($option->{-show_gis}->[0] == 1) {
			$show_gis= '-show_gis';
		}
	}
    
	# COMMAND
	$cmd= "${blast_home_plus}/bin/${blastmode} -query $infile -db $database -out $blastout_file";
	$cmd= $cmd." -num_descriptions $algn_num -num_alignments $algn_num -dust no -num_threads 6 -evalue 1e-5 $show_gis";

	print "$cmd\n";    

#   open (IN, $cmd."|");    
    system ("$cmd");
    @filepaths= ($blastout_file);
    
}elsif (exists $option->{-r}){
    $infile=$option->{-r}->[0];
    if (-f $infile){
      @filepaths=($infile);
      $head= substr $infile,-3,3;
    }else {
      print "ERROR   Not exists input file\n";
      exit;
    }
    
}elsif (exists $option->{-rd}){
    $sel= exists $option->{-sel} ? $option->{-sel}->[0] : '^\w+';
    if (-d $option->{-rd}->[0]){
       ($base,$path)= fileparse($option->{-rd}->[0]);
       $pdir= $path.$base.'/';
       $head= $base;  
       opendir (QD,"$option->{-rd}->[0]");
       @datafile= readdir QD;
       @filepaths= map{$pdir.$_}
      ##     grep{/pout$/}@datafile;
        grep{/${sel}$/}@datafile;
    }else {
      print "ERROR   Not exists input Directory\n";
      exit;
    }    
}


#================================================================================
($lnum,$filters)= &setting_parameter($option);
$llen= exists $option->{-l} ? $option->{-l}->[0] : 100;

$resultout= exists $option->{-o} ? $option->{-o}->[0] : "hitlist_blastn".$head;
open (OUT, ">$resultout");
print OUT "\#Query\tQ_Len\tQ_Des\tSbjct\tS_Len\tS_Des\tScor\tExpec";
print OUT "\tQ_Cov\tS_Cov\tCLen\tIden\tGAP\tQ_st\tQ_ed\tS_st\tS_ed\tFram\n";

#================================================================================
#     MAIN LOOP
#================================================================================
$loginfo= {};
foreach $infile_path (@filepaths){
   print "PATh $infile_path\n";

	$volume_type = 1;
    if ($volume_type == 1) {
        ### open (IN,"$infile_path");
        ### undef $/;
        ###        $blastout=<IN>;
	} else {
        ###		($blastout)= &doing_smaller_hittxt($infile_path);
	}

    #   $error_que= &main_control($blastmode,$blastout,$lnum,$filters,$llen);
       $error_que= &main_control($blastmode,$infile_path,$lnum,$filters,$llen);
   push @{$loginfo->{error}},@$error_que;
}

print "FILTER; 0;DUPLI($$filters[0])
        1;Query_COVERRAGE($$filters[1]) 2;Sbjct_COVERRAGE($$filters[2])
        3;SIMILARITY%($$filters[3]) 4;PROBABILITY($$filters[4])\n";
print "$$filters[5]\n\n";
print "NUMBER($lnum)\n";
print "LENGTH($llen)\n";
print "OUTPUTFILE; $resultout\n";

$errornum= @{$loginfo->{error}};
### print "ERROR QUERY $errornum\n\n@{$loginfo->{error}}\n";
 print "ERROR QUERY $errornum\n";

#================================================================================
#      SETTING PARAMETER
#================================================================================
sub setting_parameter {
    my ($option)=@_;
    my ($lnum,$dupl,$larea_q,$larea_s,$lhomo,$lexp,$exefil);
    my ($subreg,$subfilter,$filters);
    
    $lnum= exists $option->{-n} ? $option->{-n}->[0] : 5;
    $dupl= exists $option->{-s} ? $option->{-s}->[0] : "T";
    
    $larea_q= exists $option->{-cq} ? $option->{-cq}->[0] : 0;
    $larea_s= exists $option->{-cs} ? $option->{-cs}->[0] : 0;
    $lhomo= exists $option->{-lh} ? $option->{-lh}->[0] : 80;
    $lexp= exists $option->{-lp} ? $option->{-lp}->[0] : '1e-5';
    $exefil= exists $option->{-f} ? $option->{-f}->[0] : "T";
    
    if ($exefil eq "T") {
        $subfilter="";
        foreach ($larea_q,$larea_s,$lhomo,$lexp){
            $subreg= $_ == 0 ? '[TF]' : '[T]';
            $subfilter= $subfilter.$subreg;
        }
    } else {
        ##  $lnum= 0;
        $dupl= "F";
        $subfilter= '[TF][TF][TF][TF]';
    }
    
    $filters= [$dupl,$larea_q,$larea_s,$lhomo,$lexp,$subfilter];
    
    return ($lnum,$filters);
}

#================================================================================
#     MAIN EXTRACTION AND SUMMARIZE CONTROL
#================================================================================ 
sub main_control{
    my ($blastmode,$infile_path,$lnum,$filters,$llen)=@_;
    my ($subfasta,$count,$one_fastaout,$refhit,@error_query);
    my ($query_info,$hit,$hitresult);
    my ($hitsbj_sum,$sel_hitsbj);
  
    #    $subfasta= &ToEach_blastout($blastout);


    $|=1; ## release buffer
    
    open (IN,"$infile_path");
    @qline= ();
    $errors= [];
    $count= 1;
    while (<IN>) {
        chomp;
        if (/Effective search space used\: \d+/) {
            $one_blastout= join "\n",@qline;
            ($errors)= &sub_main_control($blastmode,$one_blastout,$lnum,$filters,$llen,$errors);
            @qline= ();
            $count++;
            &display_counter($count,[100,10,4],$count);
        } else {
            push @qline, $_;
        }
    }
    #    $one_blastout= join "\n",@qline;
    #    ($errors)= &sub_main_control($blastmode,$one_blastout,$lnum,$filters,$llen,$errors);
    
    return $errors;
}

sub sub_main_control {
    ($blastmode,$one_blastout,$lnum,$filters,$llen,$errors)= @_;
    #  foreach $one_blastout (@$subfasta){
    
    @error_query= @$errors;
    if ($one_blastout !~ /Searching\.+done/) {
        push @error_query, $1 if ($one_blastout=~ /Query= (\S+)/); 
    }

    $refhit= &addtag_blastout($one_blastout);

    ($query_info,$hit,$hitresult)= &collect_hit_blastn([@$refhit]);

    print OUT "\#\# $$query_info[0] $hit\n";

    if ($hit eq "HIT"){
        ($hitsbj_sum)= &create_sbjct_suminfo($query_info,$hitresult);
        ($sel_hitsbj)= &control_filter_each($query_info,$hitsbj_sum,$filters);
        ($sel_hitsbj)= &filter_complen_each($sel_hitsbj,$llen) if (@$sel_hitsbj>0);

        $hit= @$sel_hitsbj > 0 ?  "HIT" : "NO_HIT";
                    
    } elsif ($hit eq "NO_DATA"){
        $sel_hitsbj= [];
    } elsif ($hit eq "NO_HIT"){
        $sel_hitsbj= [];
    }
      
    $lnum= &make_table_blastn($query_info,$hit,$sel_hitsbj,$lnum);

    #    print "$blastmode Reform FINISH No \n";
    return ([@error_query]);

}

#================================================================================ 
sub display_counter {
  my ($count,$limits,$sign)=@_;
  my ($lim);

  ($subunit,$munit,$rep)= (@$limits);

  print "\-" if ($count%$subunit == 0);
  print "$sign" if ($count%($subunit*$munit) == 0);
  print "\n" if ($count%($subunit*$munit*$rep) == 0);
}

#================================================================================
#     COLLECT BLAST HIT
#================================================================================ 
sub collect_hit_blastn {
   my ($list)=@_;
   my (@qus,$quelen,$qutxt,@quss,$query,$qu_des,$z,@query_info,$hit);
   
   my (@db_des,$db,$db_des,$dblen,$dbdes);
   my (@scorlist,$scor,$expec,$i);
   my ($identH,$identC,$identP);
   my ($str,$s_start,$s_end);
   
   my ($qs,$qe,$qss,$ss,$se,$sss);
   my ($partq,$parts);
   my ($result);
   my (%hasresult,@result_db);

 ###  Main Loop   #####
	$hit="NO_DATA";
   foreach (@$list) {
     chomp;
      if (/^Query=/../^OUTT/){   ##=====  Summary  =======

           if(/^Query=/../^Length\=/){
#print "$_\n";
              @query_info=() if (/^Query=/); 
              unless (/^Length\=/){
                 push @qus,$_;
              }else{
                 /^Length\=(\d+)/;
						$quelen=$1;
						$quelen=~ s/\,//g;

						$qutxt=join " ",@qus; $qutxt=~s/ {2,}/ /g;
						@quss=split / +/, $qutxt; shift @quss; $query= shift @quss;
						$qu_des=@quss > 0 ? join " ",@quss : "x";
						@query_info=($query,$quelen,$qu_des);
						$z=0;                        
						@qus=();
						my $info= {};
              }
          } elsif (/No hits found/) {
              $hit="NO_DATA";
          }
      }elsif (/^>/.../^OUTT/) {   ##=====  Alignment  =======
           $hit="HIT";
          if (/^>/.../Length/) {   ##==== Sbjct Info  ======
              if (/^> *(\S+)/){
                 $db= $1;
                 @db_des= split;
                  if ($db_des[0]=~ /^>$db/) {
                      shift @db_des;
                  } else {
                      shift @db_des;
                      shift @db_des;
                  }
                 $dbdes="x";
            
              } elsif (/Length=(\d+)/) {
                 $dblen=$1;
                 $dbdes= join " ", @db_des if (@db_des > 0);
                 $dbdes=~s/ {2,}//g;

              } else {
                 push @db_des,$_;
              }
          }elsif (/Score/../^TEST/){    ##==== Score Info  ======
              if (/Score/) {
                 @scorlist = split;
                 $scor = $scorlist[2]; $expec = $scorlist[7];
                 $expec=$expec eq  '0.0'  ? 0 : $expec; $expec=~s /^e\-/1e\-/;
                 $i=1;

	       		} elsif (/Identities/) {
                 /Identities \= (\d+)\/(\d+) \((\d+)%\)/;
                 $identH=$1; $identC=$2; $identP = $3;
                 if (/Gaps \= (\d+)\/(\d+) \((\d+)%\)/){
                    $gapH=$1; $gapC=$2; $gapP = $3;  $gapval= "$1\/$2";
                 }else {
                    $gapH="x"; $gapC="x"; $gapP = "x";  $gapval= "0\/$identC";
                 }
              }elsif (/Strand/){
                    /Strand=Plus\/(\S+)/;  $str=$1;
                    $str= $str eq "Plus" ? '+' : '-';
              } elsif (/^Query /) {
                  ($qs,$qe) = (split)[1,3]; $qss=$qs if ($i==1);
              } elsif (/^Sbjct /) {
                ($ss,$se) = (split)[1,3]; $sss=$ss if ($i==1);  ++$i; 
              } elsif (/^TEST/){
                  $s_start= $str eq '+' ? $sss : $se;
                  $s_end= $str eq '+' ? $se : $sss;
                  $partq= ($qe-$qss+1)/($quelen)*100;  $partq= int $partq;
                  $parts= ($s_end-$s_start+1)/($dblen)*100;   $parts= int $parts;
                  $result=[$db,$dblen,$dbdes,$scor,$expec,
                                                      $partq,$parts,$qe-$qss+1,$identP,$gapval,
                                                      $qss,$qe,$s_start,$s_end,$str];
                  push @{$hasresult{$db}},$result;  ## Collect hit in a db
                  ($identH,$identC,$qs,$qe,$ss,$se,$qss,$sss)=(0,0,0,0,0,0,0,0);
             }
           }
       }
    }  ## End Foreach
    
    @result_db=();
    if ($hit eq "HIT") {
        @result_db=values %hasresult;
		#	print "$_\n" foreach (@{$info->{db}});
		#	push @result_db,$hasresult{$_} foreach (@{$info->{db}});
    }

    my($qnum,$res_num);
        $qnum=@query_info;
        $res_num= @result_db;
    #        print "QQQQQQ $qnum $res_num\n";
    
    
    return ([@query_info],$hit,[@result_db]);
}

#================================================================================
#     SUB CREATE SBJCT SUMMARY
#================================================================================
sub create_sbjct_suminfo {
    my ($query_info,$hitresult)=@_;
    my ($query,$quelen,@hitsbjs);
    my ($sbjone,$sbj,$sbjlen,@plus_sbjone,@minus_sbjone);
    my ($sbjsubone,$que_covsum,$que_sted,$sbj_covsum,$sbj_sted);
    my ($qcov_sump,$scov_sump,@sbjsubnew,$max_scor,$max_expec,$max_homo,$frame);
    my ($sbjsuminfo,$sbjsubone_num);
    
    ($query,$quelen)= @$query_info[0,1];

	my @sort_sbj= sort{$b->[0]->[3] <=> $a->[0]->[3]}@$hitresult;

    @hitsbjs=();
    foreach $sbjone (@sort_sbj){

       ($sbj,$sbjlen)= @{$sbjone->[0]}[0,1];
#print "  $sbj\n";
       @plus_sbjone=();
       @minus_sbjone=();
       @plus_sbjone= grep {$$_[14]=~ /\+/}@$sbjone;
       @minus_sbjone= grep {$$_[14]=~ /\-/}@$sbjone;

		$sumobjcov= 0;
       foreach $sbjsubone ([@plus_sbjone],[@minus_sbjone]) {
          if (@$sbjsubone > 0){
            ($que_covsum,$que_sted,$sbj_covsum,$sbj_sted)= &cal_coverarea($sbjsubone);
 					$sumobjcov += $que_covsum;
#print "     $que_covsum ";
             $qcov_sump= sprintf ("%.1f",$que_covsum/$quelen*100);
             $scov_sump= sprintf ("%.1f",$sbj_covsum/$sbjlen*100);
             @sbjsubnew= sort {$b->[3] <=> $a->[3] or $a->[4] <=> $b->[4]}@$sbjsubone;
             $max_scor= $sbjsubnew[0]->[3];
             $max_expec= $sbjsubnew[0]->[4];
             $max_homo= $sbjsubnew[0]->[8];
             $frame= $sbjsubnew[0]->[14];
             $mmgap= $sbjsubnew[0]->[9];
             $sbjsuminfo= [$sbj,"x","x",$max_scor,$max_expec,$qcov_sump,$scov_sump,
                          $que_covsum,$max_homo,$mmgap,$$que_sted[0],$$que_sted[1],
                          $$sbj_sted[0],$$sbj_sted[1],$frame];
             @sbjsubnew= sort {$a->[10] <=> $b->[10] or $a->[11] <=> $b->[11]}@sbjsubnew;
             unshift @sbjsubnew,$sbjsuminfo;
             $sbjsubone_num= @sbjsubnew;
      ##       print "$query @$sbjsuminfo\t\=> $sbjsubone_num\n";
             push @hitsbjs,[@sbjsubnew];       
          }
      }
#		print " \== $sumobjcov\n";
   }
   return [@hitsbjs];  
}
sub cal_coverarea {
    my ($onedb) =@_;
    my (@onedb_ques,@onedb_sbjs,@cover);
    my ($one,$top,$m_st,$m_ed,$start,$end);
    my (@mlens,$cov_sum,$subdb,$st,$ed);

    @onedb_ques= map {[@$_[10,11]]}
                 sort {$a->[10] <=> $b->[10]}@$onedb;
    @onedb_sbjs= map {[@$_[12,13]]}
                 sort {$a->[12] <=> $b->[12]}@$onedb;
    @cover=();
    foreach $one ([@onedb_ques],[@onedb_sbjs]) {
       $top= shift @$one;
       ($m_st,$m_ed)= @$top;
       $start= $m_st;
       @mlens=(); $cov_sum=0;
       while (@$one) {
         $subdb= shift @$one;
         ($st,$ed)= @$subdb;
        if (($m_ed > $st) && ($m_ed < $ed)) {
           $m_ed = $ed;
        } elsif ($m_ed < $st) {
           push @mlens,$m_ed-$m_st+1;
           $m_st= $st; $m_ed= $ed;
        }
      }
      $end= $m_ed;
      push @mlens,$m_ed-$m_st+1;
      $cov_sum += $_ foreach (@mlens);
      push @cover,$cov_sum,[$start,$end];
    }

    return (@cover);
}    
#================================================================================
#     SUB FIRST FILTER
#================================================================================ 
sub control_filter {
    my ($query_info,$hitsbj_sum,$filters)=@_;
    my ($query,$quelen,@sel_hitdb_sum,$sbjone,$sbjsuminfo,$res_filter);
    
    ($query,$quelen)= @$query_info[0,1];
    
    @sel_hitsbj_sum=();
    foreach $sbjone (@$hitsbj_sum){
       $sbjsuminfo= $sbjone->[0];
       $res_filter= &filterring_sbjsuminfo($query,$sbjsuminfo,$filters);
       push @sel_hitsbj_sum,$sbjone if ($res_filter==1);  

#print "SUBSUB     @$_\n" foreach (@$sbjone);
     }

     return ([@sel_hitsbj_sum]);
}
sub control_filter_each {
    my ($query_info,$hitsbj_sum,$filters)=@_;
    my ($query,$quelen,@sel_hitsbj,$sbjone,@sel_hitsbj_sub,$hitone,$res_filter);
    
    ($query,$quelen)= @$query_info[0,1];
 
    @sel_hitsbj=(); 
    foreach $sbjone (@$hitsbj_sum){
       @sel_hitsbj_sub=();
       foreach $hitone (@$sbjone) {
          $res_filter= &filterring_sbjsuminfo($query,$hitone,$filters);
          push @sel_hitsbj_sub,$hitone if ($res_filter==1);
# print "@$hitone\n" if ($res_filter==1);
#print "$res_filter @$hitone \n";
       }
       push @sel_hitsbj,[@sel_hitsbj_sub] if (@sel_hitsbj_sub>0);
     }
     return ([@sel_hitsbj]);
}

sub filterring_sbjsuminfo {
   my ($query,$sbjsuminfo,$filters)=@_;
   my ($res_filter,$sbj,@filter,$subfilter);
   my ($max_expec,$qcov_sump,$scov_sump,$max_homo);
   my ($pass);

   @filter= @$filters;
   $res_filter= 0;
   $sbj= $sbjsuminfo->[0];
   if ($filter[0] eq "T") {   ## same sample
       return $res_filter if ($query eq $sbj);
   }   
   $subfilter= pop @filter;
   
   ($max_expec,$qcov_sump,$scov_sump,$max_homo)= @{$sbjsuminfo}[4,5,6,8];
   $pass="";  
   $pass= $qcov_sump > $filter[1] ? $pass."T" : $pass."F";   ### 1 Query_COVERRAGE
   $pass= $scov_sump > $filter[2] ? $pass."T" : $pass."F";   ### 2 Sbjct_COVERRAGE
   $pass= $max_homo > $filter[3] ? $pass."T" : $pass."F";    ### 3 SIMILARITY
   $pass= $max_expec < $filter[4] ? $pass."T" : $pass."F";   ### 4 PROBABILITY
   $res_filter= 1 if ($pass =~ /^$subfilter/);    
# print "HOMO $max_homo $filter[3] E $max_expec $filter[4]  $pass $subfilter  $res_filter  ===> "; 


   return $res_filter;        
}

#================================================================================
#     SUB FILETER COMPARED LENGTH
#================================================================================
sub filter_complen {
   my ($sel_hitsbj,$llen)=@_;
   my (@sel_hitsbj_llen,$sbjone,@sbj_llen,$sbjsumone);
   
   @sel_hitsbj_llen=();
   foreach $sbjone (@$sel_hitsbj) {
      @sbj_llen=();
      $sbjsumone= shift @$sbjone;
		if ($$sbjsumone[7] >= $llen) {
      	push @sbj_llen,$_ foreach (@$sbjone);
		}
      push @sel_hitsbj_llen,[$sbjsumone,@sbj_llen] if (@sbj_llen > 0);
   }

   return ([@sel_hitsbj_llen]);
}
sub filter_complen_each {
   my ($sel_hitsbj,$llen)=@_;
   my (@sel_hitsbj_llen,$sbjone,@sbj_llen,$sbjsumone);
   
   @sel_hitsbj_llen=();
   foreach $sbjone (@$sel_hitsbj) {
      @sbj_llen=();
      $sbjsumone= shift @$sbjone;
      foreach (@$sbjone) {
         push @sbj_llen,$_ if ($_->[7] >= $llen);
      }
      push @sel_hitsbj_llen,[$sbjsumone,@sbj_llen] if (@sbj_llen > 0);
   }

   return ([@sel_hitsbj_llen]);
}
 
#================================================================================
#     SUB MAKE TABLE
#================================================================================ 
sub make_table_blastn { 
  my ($query_info,$hit,$sel_hitsbj,$lnum)=@_;
  my ($qgene,$qlen,$qdes);
  my (@sel_hitsbj_sort,$str_count,$one_sbj,$frame);
  my ($hitorgs,$hitorg);
  my ($res);
  
   ($qgene,$qlen,$qdes)= @$query_info;      
   @sel_hitsbj_sort=sort {$a->[0]->[14] cmp $b->[0]->[14] or
           $b->[0]->[3] <=> $a->[0]->[3] or
           $a->[0]->[4] <=> $b->[0]->[4]}@$sel_hitsbj;  #  Score or EXPECTATION

   $str_count={};
   $str_count->{"+"}= 0;
   $str_count->{"-"}= 0;
	$output= 0;

   if ($hit eq "HIT"){
       foreach $one_sbj (@sel_hitsbj_sort){
          $frame= $one_sbj->[0]->[14];
          shift @$one_sbj;
          if (exists $option->{-org}) {
              $hitorgs= &ext_hitorg($one_sbj->[0]->[2]);
              $hitorg=join "\t", @$hitorgs;
              push @$_, $hitorg foreach (@$one_sbj);
          }          
				$hsp_num= @$one_sbj;
			foreach (@$one_sbj) {
				$comp_len= $$_[7];
				last if (($lnum > 0) && ($lnum <= $str_count->{$frame}));
				$res=join "\t",@$_;

				next if ($comp_len < $llen);

				$str_count->{$frame}++;
				$output++;
				print OUT "$qgene\t$qlen\t$qdes\t$res\t$str_count->{$frame}\t$hsp_num\n";


###print "$res\n";
			}
		}

		if ($output == 0) {
			print OUT "$qgene\t$qlen\t$qdes\tno_hit\n";
		}
    } elsif ($hit eq "NO_DATA"){
       print OUT "$qgene\t$qlen\t$qdes\tno_data\n";
    } elsif ($hit eq "NO_HIT"){
       print OUT "$qgene\t$qlen\t$qdes\tno_hit\n";
    }else {
       print "$qgene\tNONONO\n";
    }
    print OUT "\n";
    return $lnum;
}

sub ext_hitorg {
  my ($des)=@_;
  my (@hitorgs);
  $_=$des;  @hitorgs=();
  @hitorgs=/\[(\w+\s.+?)\]/g;
  if (@hitorgs < 2){
     push @hitorgs,"","";
  }
  @hitorgs=@hitorgs[0,1];
  return [@hitorgs];
}
#================================================================================
#     SUB OTHERS
#================================================================================ 
#====================  SUB BLASTHIT  ===============================
sub ToEach_blastout {
  my ($in)=@_;
  my (@outputs);

$lenlen= length $in;
  $_=$in;
  @outputs= split/Effective search space used\: \d+/;

$lenlen1= length $outputs[0];

$numnum= @outputs;
print "QUERY NUM: $numnum   $lenlen $lenlen1\n";

##  shift @outputs;
  pop @outputs;
  
  return [@outputs];
}

sub doing_smaller_hittxt {
	my ($blastoutfile)= @_;
	my (@hitone,@hit_all,$hitonetxt,$allhittxt);

	open (BL,"$blastoutfile");
	@hitone= ();
	@hit_all= ();
	$cc= 0;
	while (<BL>) {
		chomp;
		if (/BLAST[PXN] 2\.2/) {
			if (@hitone > 0) {
				$hitonetxt= join "\n",@hitone;
				push @hit_all, $hitonetxt if ($hitonetxt !~ /\Q***** No hits found ******\E/);
				$cc++ if ($hitonetxt !~ /\Q***** No hits found ******\E/);
			}
			
			@hitone= ($_);
		} else {
			push @hitone,$_;
		}
	}
			if (@hitone > 0) {
				$hitonetxt= join "\n",@hitone;
				push @hit_all, $hitonetxt if ($hitonetxt !~ /\Q***** No hits found ******\E/);
			}

	$allhittxt= join "\n",@hit_all;

	return ($allhittxt);
}


#=============  SUB ADD TAG ====================================
sub addtag_blastout {
  my ($in)=@_;
  my ($out,@lists);
  $_=$in;
  s/\n\>/\nTEST\nOUTT\n\>/g;
  s/\n Score/\nTEST\n Score/g;

  $out=$_;
  $out= $out."\nTEST\nOUTT\n";
  @lists = split(/\n/,$out);
  return [@lists];
}



     
