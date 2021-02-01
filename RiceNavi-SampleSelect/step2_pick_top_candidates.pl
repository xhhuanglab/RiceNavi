#!/usr/bin/perl 

use strict;

my $all_sampleselect_out = shift or die "Usage: perl $0 <File_Generated_by_SampleSelect>";

open IN, $all_sampleselect_out;
open OUT, ">$all_sampleselect_out.top_candidates";

my $head;
my %count;
my @match_lines;

while(<IN>){
  chomp;
  $head = $_ if /^Indiv/;
  my @tmp = split/\t/;
  my $mark = $tmp[-1];
  #print "$mark\n";
  if($mark eq 'Y'){
    push @match_lines,$_;
    last if scalar @match_lines == 5;	
  }
}
close IN;


my $num_top_indiv = scalar @match_lines;
if($num_top_indiv == 5){
  print OUT "## Information for Top 5 best candidates ##\n";
  print OUT "$head\n";
  foreach (@match_lines){
    print OUT "$_\n";	
  }
}
elsif($num_top_indiv < 5 and $num_top_indiv > 0){
  print OUT "## Only $num_top_indiv best candidates available in this population ##\n";
  print OUT "$head\n";	
  foreach (@match_lines){
    print OUT "$_\n";	
  }
}
elsif($num_top_indiv == 0){
  print OUT "## No candidate available in this population ##\n";	
}
