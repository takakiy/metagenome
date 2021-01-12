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
               -lh Similarity% (40)               ** 0 = do not limited
               -lp Probability (1e-3)           ** 10 = do not limited
            
            -z select strand (both)/+/-      **  do not limited\n";
            
   print "\nOPTION;   -org T/F
                 -o output \n";
  
   exit;
}

foreach (@ARGV){if (/^\-\w/){$key= $_} else {push @{$option->{$key}},$_}};

if ((join " ",keys %{$option}) !~ /-[er]/) {
   print "ERROR;    Retry in adding option (-e -r -rd)\n";
   exit;
}

open (LOG,">out_exe_blastx_hit.logs");

#================================================================================
#      CONTROL EXECUTIVE MODE 
#================================================================================
$blastmode= "blastx";
use File::Basename;

	$blast_home="/home/bioresource/biotools/local/homology/blast";
	$blast_home_plus="/home/bioresource/biotools/local/homology/ncbi-blast";
	$show_gis= "";


@filepaths=();
if (exists $option->{-e}){
     $infile= $option->{-e}->[0];
    if (exists $option->{-d}){
        $database=$option->{-d}->[0];
##			$database= "${blast_home}/db/$database" if ($database =~ /^nr|kegg|uniprot$/);
			$database= "${blast_home_plus}/db/$database" if ($database =~ /^nr|kegg|uniprot$/);
    }else{
      print "Database ?\n";
      $database=<STDIN>;
      chomp $database;
			$database= "${blast_home}/db/$database" if ($database =~ /^nr|kegg|uniprot$/);
    }
    $head= (split/\//,$infile)[-1]; $head= (split/\./,$head)[0];
    $blastout_file="temp_blastxout".$head;
    # NORMAL

		$cmd= "${blast_home_plus}/bin/blastx -query $infile -db $database";
		$cmd= $cmd." -num_descriptions 200 -num_alignments 200 -seg no -num_threads 12 -evalue 1e-5 $show_gis -out $blastout_file";

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
    if (-d $option->{-rd}->[0]){
       ($base,$path)= fileparse($option->{-rd}->[0]);
       $pdir= $path.$base.'/';
       $head= $base;  
       opendir (QD,"$option->{-rd}->[0]");
       @datafile= readdir QD;
       @filepaths= map{$pdir.$_}
      #     grep{/blast/}@datafile;
     #      grep{/^KKK/}
           grep{/^\w+/}@datafile;
		@filepaths= sort @filepaths;
    }else {
      print "ERROR   Not exists input Directory\n";
      exit;
    }    
}

#================================================================================
($lnum,$filters)= &setting_parameter($option);
$lfram= exists $option->{-z} ? $option->{-z}->[0] : '\+|\-';

$resultout= exists $option->{-o} ? $option->{-o}->[0] : "hitlist_blastx".$head;
open (OUT, ">$resultout");
@title= ('#Query');
push @title,qw(Q_Len Q_Des Sbjct S_Len S_Des Scor Expec);
push @title,qw(Q_Cov S_Cov Iden Simi GAP Q_st Q_ed S_st S_ed Fram GP_N HIT_N SUB_N);
print OUT (join "\t",@title)."\n";
#================================================================================
#     MAIN LOOP
#================================================================================
$loginfo= {};
foreach $infile_path (@filepaths){
   print "PATh $infile_path\n";
   open (IN,"$infile_path");
   undef $/;
   $blastout= <IN>;
   
	print LOG "\#\# $infile_path\n";

   $error_que= &main_control($blastmode,$blastout,$lnum,$filters,$lfram);
   push @{$loginfo->{error}},@$error_que;
}

print "FILTER; 0;DUPLI($$filters[0])
        1;Query_COVERRAGE($$filters[1]) 2;Sbjct_COVERRAGE($$filters[2])
        3;SIMILARITY%($$filters[3]) 4;PROBABILITY($$filters[4])\n";
print "$$filters[5]\n\n";
print "NUMBER($lnum)\n";
print "STRAND($lfram)\n";
print "OUTPUTFILE; $resultout\n";

$errornum= @{$loginfo->{error}};
## print "ERROR QUERY $errornum\n\n@{$loginfo->{error}}\n";

print LOG "ERROR QUERY $errornum\n\n@{$loginfo->{error}}\n";

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
   $lhomo= exists $option->{-lh} ? $option->{-lh}->[0] : 40;
   $lexp= exists $option->{-lp} ? $option->{-lp}->[0] : '1e-3';
   $exefil= exists $option->{-f} ? $option->{-f}->[0] : "T";

   if ($exefil eq "T") {
      $subfilter="";
      foreach ($larea_q,$larea_s,$lhomo,$lexp){
         $subreg= $_ == 0 ? '[TF]' : '[T]';
         $subfilter= $subfilter.$subreg;
      }
   } else {
  ##    $lnum= 0;
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
  my ($blastmode,$blastout,$lnum,$filters,$lfram)=@_;
  my ($subfasta,$count,$one_fastaout,$refhit,@error_query);
  my ($query_info,$hit,$hitresult);
  my ($hitsbj_sum,$sel_hitsbj);
  
  $subfasta= &ToEach_blastout($blastout);

  $count=1;
  $|=1; ## release buffer
  @error_query=();

  foreach $one_blastout (@$subfasta){

     if ($one_blastout !~ /Searching\.+done/) {
        push @error_query, $1 if ($one_blastout=~ /Query= (\S+)/); 
     }
    
     $refhit= &addtag_blastout($one_blastout);
     ($query_info,$hit,$hitresult)= &collect_hit_blastx([@$refhit]);
     
     if ($hit eq "HIT"){
        ($hitsbj_sum)= &create_sbjct_suminfo($query_info,$hitresult);
       
        ($sel_hitsbj)= &control_filter($query_info,$hitsbj_sum,$filters);

        ($sel_hitsbj)= &filter_frame($sel_hitsbj,$lfram);
        
         $hit= @$sel_hitsbj > 0 ?  "HIT" : "NO_HIT";
                    
     } elsif ($hit eq "NO_DATA"){
         $sel_hitsbj= [];     
     }
      
     $lnum= &make_table_blastx($query_info,$hit,$sel_hitsbj,$lnum);

     $count++;
     &display_counter($count,[100,10,4],$count);
     
    # open (RES,">temp_addTag");
    # print RES join "\n",@$refhit;
  }
  print "$count\n";
  print "$blastmode Reform FINISH No \n";
  
  return [@error_query];
}

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
sub collect_hit_blastx {
   my ($list)=@_;
   my (@qus,$quelen,$qutxt,@quss,$query,$qu_des,$z,@query_info);
   
   my (@db_des,$db,$db_des,$dblen,$dbdes);
   my (@scorlist,$scor,$expec,$i);
   my ($identH,$identC,$identP,$simiH,$simiC,$simiP);
   my ($gapH,$gapC,$gapP);
   my ($fram);
   
   my (@qlist,$qs,$qe,$qss,@slist,$ss,$se,$sss);
   my ($pass,$partq,$parts,$qstart,$qend);
   my ($result);   
   my (%hasresult,@result_db);
   my ($sel_hitdb,$mark_hitdb,$hit);

	my $resultinfo= {};

   @query_info=(); 
 ###  Main Loop   #####
	$hit="NO_DATA";

   foreach (@$list) {
     chomp;
      if (/^Query=/../^OUTT/){   ##=====  Summary  =======
				if(/^Query=/../Length\=/){
              @query_info=() if (/^Query=/); 
              unless (/Length\=/){
                 push @qus,$_;
              }else{
						/Length\=(\d+)/;
						$quelen=$1;
						$quelen=~tr/,//d;
						$qutxt=join " ",@qus;
						$qutxt=~s/ {2,}/ /g;
						@quss=split / +/, $qutxt;
						shift @quss;
						$query=shift @quss;
						$qu_des=@quss > 0 ? join " ",@quss : "xx";
						@query_info=($query,$quelen,$qu_des);
						$z=0;                        
						@qus=();

						print LOG "  $query\n";
#						print "  $query\n";

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

                  push @{$resultinfo->{sbjs}},$db;

              } elsif (/Length=(\d+)/) {
                  $dblen=$1;
                  $dbdes=join " ", @db_des;
                  $dbdes=~s/ {2,}/ /g;

                  $resultinfo->{sbj}->{$db}->{len}= $dblen;
                  $resultinfo->{sbj}->{$db}->{desc}= $dbdes;

              } else {
                 push @db_des,$_;
              }
			}elsif (/Score/../^TEST/){    ##==== Score Info  ======

                if (/Score/) {
                    @scorlist = split;
                    $scor = $scorlist[2]; $expec = $scorlist[7];
                    $expec=~ s/^e\-/1e\-/; $expec=~ s/\,//;
                    $expec=$expec eq  '0.0'  ? 0 : $expec;
                    $i=1;
				} elsif (/Identities/) {
                    /Identities \= (\d+)\/(\d+) \((\d+)%\)/;
                    $identH=$1; $identC=$2; $identP = $3;
                    /Positives \= (\d+)\/(\d+) \((\d+)%\)/;
                    $simiH=$1; $simiC=$2; $simiP = $3;
                    if (/Gaps \= (\d+)\/(\d+) \((\d+)%\)/){
                        $gapH=$1; $gapC=$2; $gapP = $3;  $gapval= "$1\/$2";
                    }else {
                        $gapH="x"; $gapC="x"; $gapP = "x";  $gapval= "0\/$identC";
                    }
				} elsif (/Frame \= ([\+\-]\d)/){
                    $fram=$1;
				} elsif (/^Query /) {
                    ($qs,$qe) = (split)[1,3];
                    ($qs,$qe)= map{tr/,//d; $_}($qs,$qe);
                    $qss=$qs if ($i==1);
				} elsif (/^Sbjct /) {
                    ($ss,$se)= ($1,$2) if (/^Sbjct *(\d+).+ +(\d+)/);
                    ($ss,$se)= map{tr/,//d; $_}($ss,$se);
                    $sss= $ss if ($i==1);
                    ++$i;
				} elsif (/^TEST/){
                    $qstart = $fram=~/\+/ ? $qss : $qe;
                    $qend = $fram=~/\+/ ? $qe : $qss;
                    $partq= ($qend-$qstart+1)/($quelen)*100;  $partq= int $partq;
                    $parts = ($se-$sss+1)/($dblen)*100; $parts=int $parts;
                    ## 0 1 2 3 4 - 5 6 7 8 9 - 10 11 12 13 14
                    $result=[$db,$dblen,$dbdes,$scor,$expec,
                                                  $partq,$parts,$identP,$simiP,$gapval,
                                                  $qstart,$qend,$sss,$se,$fram];

                  push @{$resultinfo->{sbj}->{$db}->{hits}},$result;

							push @{$hasresult{$db}},$result;
                                  
                  ($identH,$identC, $simiH,$simiC,$gapH,$gapC,$gapP)=qw(0 0 0 0 0 0 0);
                  ($qs, $qe, $ss, $se, $qss, $sss)=(0,0,0,0,0,0);
              }    
          }
       }        
   }  ## End Foreach 
   
    @result_db=();
    if ($hit eq "HIT") {
##       @result_db=values %hasresult;

		push @result_db,$resultinfo->{sbj}->{$_}->{hits} foreach (@{$resultinfo->{sbjs}});

    }
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
    @hitsbjs=();
    foreach $sbjone (@$hitresult){
       ($sbj,$sbjlen)= @{$sbjone->[0]}[0,1];
       @plus_sbjone=();
       @minus_sbjone=();
       @plus_sbjone= grep {$$_[14]=~ /\+/}@$sbjone;
       @minus_sbjone= grep {$$_[14]=~ /\-/}@$sbjone;
##$snum= @$sbjone;
## print "$sbj $snum\n";
       foreach $sbjsubone ([@plus_sbjone],[@minus_sbjone]) {
          if (@$sbjsubone > 0){
             ($que_covsum,$que_sted,$sbj_covsum,$sbj_sted)= &cal_coverarea($sbjsubone);
             $qcov_sump= sprintf ("%.1f",$que_covsum/$quelen*100);
             $scov_sump= sprintf ("%.1f",$sbj_covsum/$sbjlen*100);
             @sbjsubnew= sort {$b->[3] <=> $a->[3] or $a->[4] <=> $b->[4]}@$sbjsubone;
             $max_scor= $sbjsubnew[0]->[3];
             $max_expec= $sbjsubnew[0]->[4];
             $max_homo= $sbjsubnew[0]->[8];
             @sbjsubnew= sort {$a->[10] <=> $b->[10] or $a->[11] <=> $b->[11]}@sbjsubnew;
             $frame= $sbjsubnew[0]->[14];
             $mmgap= $sbjsubnew[0]->[9];
             $sbjsuminfo= [$sbj,"x","x",$max_scor,$max_expec,$qcov_sump,$scov_sump,
                          $que_covsum,$max_homo,$mmgap,$$que_sted[0],$$que_sted[1],
                          $$sbj_sted[0],$$sbj_sted[1],$frame];
             unshift @sbjsubnew,$sbjsuminfo;
             $sbjsubone_num= @sbjsubnew;
# print LOG "$query @$sbjsuminfo\t\=> $sbjsubone_num\n";
# print LOG "   @$_\n" foreach (@sbjsubnew);
             push @hitsbjs,[@sbjsubnew];

#print "\n";   
          }
      }
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
   return $res_filter;        
}
#================================================================================
#     SUB FILETER COMPARED LENGTH
#================================================================================
sub filter_frame {
   my ($sel_hitsbj,$lfram)=@_;
   my (@sel_hitsbj_lfram,$sbjone,@sbj_lfram,$sbjsumone);
   
   @sel_hitsbj_lfram=();
   foreach $sbjone (@$sel_hitsbj) {
      @sbj_lfram=();
      $sbjsumone= shift @$sbjone;

 $sbj= $sbjone->[0]->[0];
## $snum= @$sbjone;
##print "$sbj $snum $lfram\n";

      foreach (@$sbjone) {
         push @sbj_lfram,$_ if ($_->[14] =~ /$lfram/);
      }

## $ssnum= @sbj_lfram;
## print "$sbj $ssnum $lfram\n";

      push @sel_hitsbj_lfram,[$sbjsumone,@sbj_lfram] if (@sbj_lfram > 0);
   }

   return ([@sel_hitsbj_lfram]);
}
 
 
#================================================================================
#     SUB MAKE TABLE
#================================================================================ 
sub make_table_blastx { 
  my ($query_info,$hit,$sel_hitsbj,$lnum)=@_;
  my ($qgene,$qlen,$qdes);
  my ($grpcount,$sbj_subsort,$sub_sbjs,$hitcount,$one_sbj,$frame);
  my ($hitorgs,$hitorg);
  my ($res);
  
   ($qgene,$qlen,$qdes)= @$query_info; 
    $grpcount=0;
    if ($hit eq "HIT"){
       ($sbj_subsort)= &sort_qarea_blastx($sel_hitsbj,$qgene);
       foreach $sub_sbjs (@$sbj_subsort){  ##  Each Area in a Strand
          $hitcount=0;
          $grpcount++;
          foreach $one_sbj (@$sub_sbjs) {   ## Each Sbj of Each Area in a Strand
my  $onenum= @$one_sbj;
my  $sbj= $one_sbj->[0]->[0];
# print "$sbj $onenum\n";
             $frame= $one_sbj->[0]->[14];
             shift @$one_sbj;
             $hitcount++;       
             if (exists $option->{-org}) {
                $hitorgs= &ext_hitorg($one_sbj->[0]->[2]);
                $hitorg=join "\t", @$hitorgs;
                push @$_, $hitorg foreach (@$one_sbj);
             }
             $subhitcount=0;
             foreach (@$one_sbj) {
                last if (($lnum > 0) && ($lnum < $hitcount)); 
                $subhitcount++;                     
                $res=join "\t",@$_;
                print OUT "$qgene\t$qlen\t$qdes\t$res\t$grpcount\t$hitcount\t$subhitcount\n";
             }
          }
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

sub sort_qarea_blastx {
   my ($sel_hitsbj,$qgene)= @_;
   my ($sbj_strgrp,@sbj_subsort);
   my ($strand,@strgrp_one,$subgrpid,$one_sbj,$qstart,$qend,$qterm);
   my (@sub_sbjs);
   
   $sbj_strgrp={};
   ## GRP for STRAND
	foreach (@$sel_hitsbj) {
		$mstrand= $_->[0]->[14] > 0 ? '+' : '-';
   	push @{$sbj_strgrp->{$mstrand}},$_;
	}
## print "EE $_->[0]->[14]\n" foreach (@$sel_hitsbj);

   @sbj_subsort=();
   foreach $strand (sort keys %{$sbj_strgrp}) {
      @strgrp_one= sort {$a->[0]->[10] <=> $b->[0]->[10] or
           $a->[0]->[11] <=> $b->[0]->[11]}@{$sbj_strgrp->{$strand}};
      delete $sbj_strgrp->{$strand};
    
      $subgrpid=0;
      $qterm=0;
      foreach $one_sbj (@strgrp_one){
#print "@{$one_sbj->[0]}\n";
	##		($one_sbj)= &merge_over_subgrp($one_sbj);
			($one_sbj)= &merge_over_subgrp_2($one_sbj);

         $qstart= $one_sbj->[0]->[10];
         $qend= $one_sbj->[0]->[11];

         if ($qterm < $qstart) {
            $subgrpid++;
            $sbj_strgrp->{$strand}->{$subgrpid}= [$one_sbj];
            $qterm= $qend;
## print "EEE $qstart $qend\n" if ($qgene eq "KASF001_A05");
         } elsif ($qterm >= $qstart) {
            push @{$sbj_strgrp->{$strand}->{$subgrpid}},$one_sbj;
            $qterm= $qend if ($qterm < $qend);
         }
##  print "$one_sbj->[0]->[0] $qstart $qend N $qterm $subgrpid $dd $strand\n" if ($qgene eq "KASF001_A05");
      }
      @subgrpids= sort {$a <=> $b}(keys %{$sbj_strgrp->{$strand}}); 
      foreach $subgrpid (@subgrpids) {
         @sub_sbjs= sort {$b->[0]->[3] <=> $a->[0]->[3] or
             $a->[0]->[4] <=> $a->[0]->[4]}@{$sbj_strgrp->{$strand}->{$subgrpid}};
         push @sbj_subsort,[@sub_sbjs];
      }
        
   }
   
   return ([@sbj_subsort]);
}

sub merge_over_subgrp {
	my ($one_sbjhit)= @_;
	my ($sbjone_sum,$sub_grpno_info,@sort_sbjone,$sub_grp_top,$m_st,$m_ed,$t_st,$t_ed,$sub_covrate);

			$sbjone_sum= shift @$one_sbjhit;
			
			$sub_grpno_info= {};

			@sort_sbjone= sort {$a->[10] <=> $b->[10] or
           $a->[11] <=> $b->[11]}@$one_sbjhit;

			$sub_grp_top= shift @sort_sbjone;
			($m_st,$m_ed)= @{$sub_grp_top}[10,11];
			$sub_grpno= 1;
			push @{$sub_grpno_info->{$sub_grpno}},$sub_grp_top;
			foreach (@sort_sbjone) {
				($t_st,$t_ed)= @{$_}[10,11];
				if ($m_ed > $t_st) {
					if ($m_ed <= $t_ed) {
						push @{$sub_grpno_info->{$sub_grpno}},$_;
					} else {
						$sub_covrate= ($m_ed-$t_st+1)/($m_ed-$m_st+1);
						if ($sub_covrate > 0.7) {
							push @{$sub_grpno_info->{$sub_grpno}},$_;
							($m_st,$m_ed)= ($m_st,$t_ed);
						} else {
							$sub_grpno++;
							push @{$sub_grpno_info->{$sub_grpno}},$_;
							($m_st,$m_ed)= ($m_st,$t_ed);
						}
					}
				} else {
					$sub_grpno++;
					push @{$sub_grpno_info->{$sub_grpno}},$_;
					($m_st,$m_ed)= ($t_st,$t_ed);
				}

## print "$$sbjone_sum[0] $m_st,$m_ed <=> $t_st,$t_ed $sub_grpno\n";

			}


	@new_sbjone= ($sbjone_sum);
	foreach $sub_grpno (sort{$a<=>$b}(keys %{$sub_grpno_info})) {
		@sub_grpno_hits= sort{$b->[3]<=>$a->[3]}@{$sub_grpno_info->{$sub_grpno}};
		push @new_sbjone,$sub_grpno_hits[0];
	}

	return ([@new_sbjone]);

}

sub merge_over_subgrp_2 {
	my ($one_sbjhit)= @_;
	my ($sbjone_sum,$sub_grpno_info,@sort_sbjone,$sub_grp_top,$m_st,$m_ed,$t_st,$t_ed,$sub_covrate);

	$sbjone_sum= shift @$one_sbjhit;
			
	$sub_grpno_info= {};

	# score sort
	@sort_sbjone= sort{$a->[4] <=> $b->[4]}@$one_sbjhit;

	$sub_grp_top= shift @sort_sbjone;
	push @{$sub_grpno_info->{selhit}->{1}},$sub_grp_top;
	push @{$sub_grpno_info->{areas}},[@{$sub_grp_top}[10,11]];

	$sub_grpno= 1;
	$lim_ov_rate= 0.3;

	foreach $hitseg (@sort_sbjone) {
		($sbj,$t_st,$t_ed)= @{$hitseg}[0,10,11];
		$uniq_frg= 1;
print LOG "$sbj\n";
		foreach $oneseg (@{$sub_grpno_info->{areas}}) {
			($m_st,$m_ed)= @$oneseg;
			$m_len= $m_ed-$m_st+1;
			$left_ov= $t_ed-$m_st;
			$right_ov= $m_ed-$t_st;
			$left_rate= $left_ov/$m_len;
			$right_rate= $right_ov/$m_len;
print LOG "  $m_st,$m_ed  $t_st,$t_ed $left_rate $right_rate\n";
			if ($left_rate < 0 && $right_rate > 1) {
				$uniq_frg= 1;   ## LEFT overlap
			} elsif ($left_rate > 1 && $right_rate < 0) {
				$uniq_frg= 1;   ## RIGHT overlap
			} elsif ($left_rate < $lim_ov_rate && $right_rate > 1) {
				$uniq_frg= 1;   ## Left over					
			} elsif ($left_rate > 1 && $right_rate < $lim_ov_rate) {
				$uniq_frg= 1;	 ## Right over
			} else {
				$uniq_frg= 0;  ## Overlap
				last;
			}
print LOG "$uniq_frg $sub_grpno\n";
		}
		if ($uniq_frg == 1) {
			$sub_grpno++;
			push @{$sub_grpno_info->{selhit}->{$sub_grpno}},$hitseg;
			push @{$sub_grpno_info->{areas}},[$t_st,$t_ed];
		}
	}

	@new_sbjone= ($sbjone_sum);
	foreach $sub_grp (1..$sub_grpno) {
		push @new_sbjone,@{$sub_grpno_info->{selhit}->{$sub_grp}};
	}

	return ([@new_sbjone]);

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

