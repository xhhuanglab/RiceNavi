#!/usr/bin/perl

my $QTNlist;
my $diff_samples_file;

if(@ARGV < 2){
  die "Usage: perl $0 <Diff_samples> <Selected_QTNlist>";	
}else{
  ($diff_samples_file,$QTNlist) = @ARGV;
}

my $sample_type = "QTNpickLib/QTNlib.404lines_rice.type";

open TYPE, $sample_type;
open DIFF_GENO, $diff_samples_file;
open QTN, $QTNlist;

open OUT, ">Candidate.Improvement_samples";

my @selected_alleles;
while(<QTN>){chomp; push @selected_alleles, $_;}
close QTN;
my %selected_alleles = map {$_,1} @selected_alleles;
my $NO_genes_selected = scalar @selected_alleles;

my %sample2type;
my %count_samples;

while(<TYPE>){chomp; $sample2type{$1} = $2 if /(.+?)\t(.+)/;}
close TYPE;

print OUT "GeneName\tChrom\tSite\tGenotyping_Method\tNipAllele(0|0)\tAltAllele\tAllele_Effect(AltAllele)\tNo. samples with different alleles\tQTNlib samples with different alleles\n";
while(<DIFF_GENO>){
  chomp;
  my @tmp = split/\t/;
  #print $tmp[0]."\t".$tmp[1]."\t".$tmp[2]."\n";
  if(defined $selected_alleles{$tmp[0]."\t".$tmp[1]."\t".$tmp[2]}){
    my @samples = split/\|/,$tmp[-1];
    print OUT "$_\n";
    foreach my $sample (@samples){
    	if($sample =~ /(.+?)\(/){
    	  $count_samples{$1}++;	
    	}else{
     	  $count_samples{$sample}++;   	  	
    	}
    }
  }
}
close DIFF_GENO;

my $count_overlap_lines = 0;
my %overlap_samples;

foreach my $sample (sort keys%count_samples){
  if($count_samples{$sample} == $NO_genes_selected){
  	$count_overlap_lines++;
    $overlap_samples{$sample}++;
  }
}
print OUT "\n##############################\n";
print OUT "The samples which harbor the selected alternative alleles are listed below.\nNo. candidate accessions: $count_overlap_lines\n";

if($count_overlap_lines){
  foreach (sort keys%overlap_samples){
    print OUT "$_\t$sample2type{$_}\n";	
  }	
}