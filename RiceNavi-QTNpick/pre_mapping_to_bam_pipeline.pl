#/usr/bin/perl 
## Before runnning this script, make sure the reference genome has been 
## indexed by Bowtie2-build, for example: Bowtie2-build Rice_MSUv7.fa Rice_MSUv7

use strict;
use FindBin;

my ($cfg,$fq1,$fq2,$sample_prefix);
if(@ARGV < 4){
  die "Usage: perl $0 <RiceNavi-QTNpick.cfg> <fq1> <fq2> <sample_prefix(sample1)>";	
}else{
  ($cfg,$fq1,$fq2,$sample_prefix) = @ARGV;	
}

my ($bowtie2,$gatk3,$samtools,$ref_genome);

system(qq(dos2unix $cfg));
open CFG, $cfg;

while(<CFG>){
  chomp;
  s/\s+//;
  s/\r\n/\n/g;
  $bowtie2 = $1 if /bowtie2=(.+)/i;
  $gatk3 = $1 if /GATK3=(.+)/i;
  $samtools = $1 if /samtools=(.+)/i;
  $ref_genome = $1 if /Ref_genome=(.+)/i;
}
print "$ref_genome\n";
print "$samtools\n";
my $gatk_tmpfile = "$sample_prefix/GATK_tmp";

mkdir $gatk_tmpfile;
mkdir $sample_prefix;
open OUT, ">$sample_prefix/$sample_prefix.time";

my $genome_prefix = $1 if $ref_genome =~ /(.+)\.fa/;


##### bowtie2 mapping ########
#############################
print OUT "Mapping Started: ".localtime()."\n";
system(qq($bowtie2 -p 30 -x $genome_prefix -1 $fq1 -2 $fq2 --rg-id $sample_prefix --rg "PL:ILLUMINA" --rg "SM:$sample_prefix" -S $sample_prefix/$sample_prefix.sam));
system(qq($samtools view -bS $sample_prefix/$sample_prefix.sam > $sample_prefix/$sample_prefix.bam));
system(qq($samtools sort --threads 30 -o $sample_prefix/$sample_prefix.sorted.bam $sample_prefix/$sample_prefix.bam));


my $sortedbam = $sample_prefix.'.sorted.bam';
system(qq($samtools index $sample_prefix/$sortedbam));
system(qq($samtools flagstat $sample_prefix/$sortedbam > $sample_prefix/$sample_prefix.flagstat));


##  quality control #######
###########################
print OUT "remove duplicated Started: ".localtime()."\n";
system(qq($samtools rmdup -sS $sample_prefix/$sample_prefix.sorted.bam $sample_prefix/$sample_prefix.rmdup.bam));
system(qq($samtools index $sample_prefix/$sample_prefix.rmdup.bam));

print OUT "Realign Started: ".localtime()."\n";
system(qq(java -Xmx20g -Djava.io.tmpdir=$gatk_tmpfile -jar $gatk3 -R $ref_genome -T RealignerTargetCreator -o $sample_prefix/$sample_prefix.realn.intervals -I $sample_prefix/$sample_prefix.rmdup.bam));

## generating realgined bam
system(qq(java -Xmx20g -Djava.io.tmpdir=$gatk_tmpfile -jar $gatk3 -R $ref_genome -T IndelRealigner -targetIntervals $sample_prefix/$sample_prefix.realn.intervals -I $sample_prefix/$sample_prefix.rmdup.bam -o $sample_prefix/$sample_prefix.realn.bam));

## remove useless large files
&remove_file($sample_prefix,$sample_prefix);
print OUT "$sample_prefix finished: ".localtime()."\n";


###### SUB ########
###################

sub remove_file {
  my ($sample_prefix,$prefixname) = @_;
  if(-f "$sample_prefix/$prefixname.sorted.bam"){
    system(qq(rm -f $sample_prefix/$prefixname.sam));	
    system(qq(rm -f $sample_prefix/$prefixname.bam));
    system(qq(rm -f $sample_prefix/$prefixname.rmdup.bam));
  }	
}
