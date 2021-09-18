# fq1, fq2: ST fastq files
# barcode_file: barcode file provided by ST
# folder_out: folder to put output fastq.gz files
#!/usr/bin/perl
use strict;
use warnings;

my ($fq1,$fq2,$barcode_file,$folder_out)=@ARGV;
my (@items,$i,$line,$fq,$line1,%barcodes,$barcode);

# read barcodes
open(FILE,$barcode_file);
while ($line=<FILE>)
{
  @items=split("\t",$line);
  $barcodes{$items[0]}={"n"=>0,"fq"=>""};
}
close(FILE);

# clean output folder
unless (-d $folder_out) {mkdir($folder_out);}
foreach (keys %barcodes) {unlink($folder_out."/".$_.".fastq");}

# read fastq files
open(FILE1,"zcat ".$fq1." |");
open(FILE2,"zcat ".$fq2." |");

while ($line=<FILE1>)
{
  # read barcode
  $line=<FILE1>;
  $barcode=substr $line,0,18;
  $line=<FILE1>;
  $line=<FILE1>;

  # read fq
  $line1=<FILE2>;
  $fq=<FILE2>;
  $line1.=rev_complement($fq);
  $line1.=<FILE2>;
  $fq=<FILE2>;
  $fq=~s/\n//;
  $line1.=reverse($fq)."\n";

  # output
  if (exists $barcodes{$barcode})
  {
    $barcodes{$barcode}->{"n"}++;
    $barcodes{$barcode}->{"fq"}.=$line1;
    if ($barcodes{$barcode}->{"n"}>=1000)
    {
      open(FILE_OUT,">>".$folder_out."/".$barcode.".fastq");
      print FILE_OUT $barcodes{$barcode}->{"fq"};
      close(FILE_OUT);
      $barcodes{$barcode}->{"n"}=0;
      $barcodes{$barcode}->{"fq"}="";
    }
  }
}

close(FILE1);
close(FILE2);

# clean up
foreach $barcode (keys %barcodes)
{
  if ($barcodes{$barcode}->{"n"}>0)
  {
    open(FILE_OUT,">>".$folder_out."/".$barcode.".fastq");
    print FILE_OUT $barcodes{$barcode}->{"fq"};
    close(FILE_OUT);
  }
  if (-e $folder_out."/".$barcode.".fastq") {system("gzip ".$folder_out."/".$barcode.".fastq");}
}

sub rev_complement
{
  $line=$_[0];
  my ($i,$read);

  $line=~s/\n//;
  $line=uc(reverse($line));
  
  $read="";
  for $i (split //, $line)
  {
    if ($i eq "A") {
      $read.="T";
    }elsif ($i eq "T") {
      $read.="A";
    }elsif ($i eq "C") {
      $read.="G";
    }elsif ($i eq "G") {
      $read.="C";
    }else {
      $read.="N";
    }
  }
  
  $read."\n";
}
