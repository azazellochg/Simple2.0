#!/usr/bin/perl

use warnings;
use strict;
use Env;
use Cwd qw(getcwd);
use List::Util qw(max min);

# ASSUMES SIMPLEBIN ENVIRONMENT VARIABLE SET
# ASSUMES SIMPLESYS ENVIRONMENT VARIABLE SET

# check so that at least one argument
if(scalar(@ARGV) == 0){die "need at least one argument: prg=<simple program name>\n"};

# parse the command line
# $cmd_string is the command line
# %name_value is the hash with input data
my $cmd_string = '';
my %name_value;
foreach my $cmd (@ARGV){
    chomp($cmd);
    my @tmp = split('=', $cmd);
    if( $tmp[0] ne 'prg' ){
        $cmd_string = $cmd_string.$cmd.' ';
    }
    $name_value{$tmp[0]} = $tmp[1];
}
$cmd_string =~ s/\s+$//; # removes trailing whitespace
$cmd_string =~ s/queue\=//;
$cmd_string =~ s/grant//;
$cmd_string =~ s/igbmc2015//;

# check so that prg is set (program name)
if(defined $name_value{'prg'}){
}else{
    die "prg=<simple program name> need to be set\n";      
}

# check so that queue is set (queue name)
if(defined $name_value{'queue'}){
}else{
    die "queue=<grant|igbmc2015> need to be set\n";      
}

# set account name for the queue
if( $name_value{'queue'} eq 'grant' ){
    $name_value{'account'} = 'g2014a90';
}elsif( $name_value{'queue'} eq 'igbmc2015' ){
    $name_value{'account'} = 'qosigbmc2015';
}else{
    die "Unknown queue!\n";
}

# set the check_nptcls program
$name_value{'check_nptcls'} = $ENV{SIMPLEBIN}.'/simple_check_nptcls';

# make absolute path to the program & set execution directory
my $abs_prg = $ENV{SIMPLEBIN}.'/'.$name_value{'prg'};
delete $name_value{'prg'};
$name_value{'prg'} = $abs_prg;
my $execdir = getcwd();

# if only program name is given, then execute
if(scalar(@ARGV) == 1){
    system("$name_value{prg}");
    die;
}

# set number of states (by counting input volumes)
my @keys = keys %name_value;
my $nstates = 0;
foreach my $key(@keys){
    if($key =~ /vol\d/ ){$nstates++};   
}
$name_value{'nstates'} = max 1,$nstates;

# set time per image
$name_value{'time_per_image'} = 60*$name_value{'nstates'};

# set the number of particles
$name_value{'nptcls'} = exec_check_nptcls();

# check so that npart is set (number of partitions)
if(defined $name_value{'npart'}){
    if($name_value{'npart'} < 1){
        die "npart need to be > 0\n";   
    }
}else{
    die "npart need to be set\n";      
}

# check so that nthr is set (number of threads)
if(defined $name_value{'nthr'}){
}else{
    $name_value{'nthr'} = 1;   
}

# parallel execution
my $stop;
my $start;
my $ptcls_per_part = int($name_value{'nptcls'}/$name_value{'npart'});
for(my $i=1; $i<=$name_value{'npart'}; $i++){
    $stop     = min $i*$ptcls_per_part,$name_value{'nptcls'};
    if($i == $name_value{'npart'}){$stop = $name_value{'nptcls'}};
    $start    = max 1,$i*$ptcls_per_part-$ptcls_per_part+1;
    my $time  = $name_value{'time_per_image'}*($stop-$start+1);
    my $hours = int($time/3600);
    my $min   = (($time/60)%60);
    open(FHANDLE, ">script_$i");
    print FHANDLE "#!/bin/bash
#SBATCH -n $name_value{'nthr'}
#SBATCH -J distr_simple_$i
#SBATCH -p $name_value{'queue'}
#SBATCH -A $name_value{'account'}
cd $execdir
$name_value{'prg'} $cmd_string fromp=$start top=$stop part=$i outfile=algndoc_$i.txt > OUT$i\nexit\n";
    close(FHANDLE);
    chmod 0777, 'script_'.$i;
    if( $ENV{SIMPLESYS} eq 'LOCAL' ){
        system("nohup ./script_$i &\n");
    }elsif( $ENV{SIMPLESYS} eq 'HPC' ){
        system("sbatch ./script_$i");
    }else{
        die "Unrecognized system!\n";
    }
}

sub exec_check_nptcls{
    my$nptcls;
    my$check_nptcls_cmd_string = $name_value{'check_nptcls'}.' ';
    $check_nptcls_cmd_string = add2string('stk', $check_nptcls_cmd_string);
    $check_nptcls_cmd_string = add2string('box', $check_nptcls_cmd_string);
    $check_nptcls_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $check_nptcls_cmd_string, "\n";
    my$nptcls_info = `$check_nptcls_cmd_string`;
    if( $nptcls_info =~ />>>\sNPTCLS\:\s+(\d+)/ ){
        $nptcls = $1;   
    }
    return $nptcls;
}

sub add2string{
    my$var = shift;
    my$str = shift;
    if( defined($name_value{$var}) ){
        $str = $str.$var.'='.$name_value{$var}.' ';
    }
    return $str;
}
