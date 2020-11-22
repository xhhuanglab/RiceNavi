#!/usr/bin/perl

################################################################################################
# Usage: perl step3_samples_with_diff_CausalVar.pl <Sample.RiceNavi_Causal_Var.site.geno>
################################################################################################


my $sample_geno = shift or die "Please Input <Sample.RiceNavi_Causal_Var.site.geno>";
my $popdb_SNP = "QTNpickLib/RiceNavi_QTNLib.genomatrix"; 

my %samples;
my %num2sample;
my %snp2method;
my (%refbase,%altbase);
my %posi2genename;
my %samplegeno;

open POPSNP, $popdb_SNP;
open GENO, $sample_geno;
open OUT, ">$sample_geno.samples";

my $head = <POPSNP>; chomp $head;
my @head = split/\t/,$head;
for my $i (5..$#head){
  $num2sample{$i} = $head[$i];
  $samples{$head[$i]}++;
}
while(<POPSNP>){
  chomp;
  my @tmp = split/\t/;
  my ($genename,$chrom,$snp_posi,$method,$refbase,$altbases) = @tmp[0,1,2,3,4,5];
  my $snploci = $chrom."\t".$snp_posi;
  $posi2genename{$chrom."\t".$snp_posi} = $genename;
  ($refbase{$snploci},$altbase{$snploci}) = ($refbase,$altbases);
  $snp2method{$snploci} = $method; 
  for my $i (6..$#tmp){
  	 next if $tmp[$i] eq '.|.';
     $samplegeno{$num2sample{$i}}{$snploci} = $tmp[$i];	
  }	
}
close POPSNP;


## variants in Sample.RiceNavi_Causal_Var.site.geno file ##
my $geno_header = <GENO>; chomp $geno_header;
print OUT "GeneName\t$geno_header\tNo. samples with different alleles\tQTNlib samples with different alleles\n";
while(<GENO>){
  chomp;
  next if /^\#/;
  my @tmp = split/\t/;
  my ($chrom,$site,$method,$refbaes,$altbase,$altfunc,$samplebase) = @tmp;
  my $chrom_site = "$chrom\t$site";
  my $candidate_samples = '';
  my $count_samples = 0;
  foreach my $pop_sample (sort keys%samples){
  	my $pop_sample_geno = $samplegeno{$pop_sample}{$chrom."\t".$site} || 'NA';
    if($snp2method{$chrom_site} ne 'Manta'){
  	  if($pop_sample_geno =~ /(.+)\|(.+)/){
  	  	my ($pop_sample_allele1,$pop_sample_allele2) = ($1,$2);
        if($pop_sample_allele1 eq $pop_sample_allele2
        and $pop_sample_geno ne $samplebase){
        	$count_samples++;
          $candidate_samples .= $pop_sample."\|";
        }
        elsif($pop_sample_allele1 ne $pop_sample_allele2
        and $pop_sample_geno ne $samplebase){
        	$count_samples++;
          $candidate_samples .= $pop_sample."(het)\|";	
        }
      }
    }else{
      if($pop_sample_geno ne $samplebase and $pop_sample_geno =~ /\|/){
        $count_samples++;
        $candidate_samples .= $pop_sample."\|";
      }
    }
  }
  chop $candidate_samples;
  if($samplebase eq '.|.'){
    print OUT "$posi2genename{$chrom_site}\t$_\tNA\tNA\n"; 	
  }else{
    print OUT "$posi2genename{$chrom_site}\t$_\t$count_samples\t$candidate_samples\n";    	
  }
}
close GENO;
