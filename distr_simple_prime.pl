#!/usr/bin/perl

use warnings;
use strict;
use Env;
use Cwd qw(getcwd);
use List::Util qw(max min);

# ASSUMES SIMPLEBIN ENVIRONMENT VARIABLE SET
# ASSUMES SIMPLESYS ENVIRONMENT VARIABLE SET

# global variables
my%name_value;
my$folder;
$name_value{'prime'}          = $ENV{SIMPLEBIN}.'/simple_prime';
$name_value{'merge_algndocs'} = $ENV{SIMPLEBIN}.'/simple_merge_algndocs';
$name_value{'volassemble'}    = $ENV{SIMPLEBIN}.'/simple_volassemble';
$name_value{'automask'}       = $ENV{SIMPLEBIN}.'/simple_automask';
$name_value{'check_conv'}     = $ENV{SIMPLEBIN}.'/simple_check_conv';
$name_value{'check_nptcls'}   = $ENV{SIMPLEBIN}.'/simple_check_nptcls';
$name_value{'execdir'}        = getcwd();

# make command line instructions
my$instructions = `$name_value{'prime'}`;
$instructions =~ s/SIMPLE_PRIME/DISTR_SIMPLE_PRIME/;
$instructions =~ s/\[nthr\=\<nr of OpenMP threads\{1\}\>\]/npart\=\<nr of partitions\>/;
my@should_not_be_there = glob("recvol_state*");
if( scalar(@should_not_be_there) >= 1 ){
    die "ERROR: files match the pattern recvol_state*, remove these from cwd\n";
}

# if too few argument given, print instructions
if(scalar(@ARGV) < 5){
  print $instructions;
}

# parse the command line
# %name_value is the hash with input data
foreach my $cmd (@ARGV){
    chomp($cmd);
    my @tmp = split('=', $cmd);
    $name_value{$tmp[0]} = $tmp[1];
}

# check so that a volume is inputted
if( !defined($name_value{'vol1'}) ){
    die "Need at least one starting volume for prime exec!\n";   
}

# set the number of particles
$name_value{'nptcls'} = exec_check_nptcls();

# check so that npart is set (number of partitions)
if( !defined $name_value{'npart'} ){
    die "npart need to be set\n";
}

# check so that nspace is set (number of projection directions)
if( !defined $name_value{'nspace'} ){
    $name_value{'nspace'} = 1000;   
}

# check so that nthr is set (number of threads)
if( !defined $name_value{'nthr'} ){
    $name_value{'nthr'} = 1;   
}

# check so that fstep is set
if( !defined $name_value{'fstep'} ){
    $name_value{'fstep'} = 1;   
}

# check so that dynlp is set 
if( !defined($name_value{'dynlp'}) ){
    $name_value{'dynlp'} = 'no';
}

# determine initial and final Fourier index
if( $name_value{'dynlp'} eq 'yes' ){
    # determine the initial low-pass limit
    if( !defined($name_value{'find'}) ){
        $name_value{'find'} = get_find(($name_value{'box'}*0.7*$name_value{'smpd'})/8.0);
        my$lp = get_lp($name_value{'find'});
#        print "LPLOW :$lp\n";
    }
    # determine the final low-pass limit
    if( defined($name_value{'fmax'}) ){
        delete $name_value{'fmax'};
    }
    $name_value{'fmax'} = get_find(($name_value{'box'}*0.7*$name_value{'smpd'})*0.07);
    my$lp = get_lp($name_value{'fmax'});
#    print "LPHIGH: $lp\n";
}

# check so that trs is set (origin shift)
if( !defined $name_value{'trs'} ){
    $name_value{'trs'} = 0;   
}

# check so that maxits is set 
if( !defined($name_value{'maxits'}) ){
    $name_value{'maxits'} = 50;
}

# set number of states (by counting input volumes)
my @keys    = keys %name_value;
my $nstates = 0;
foreach my $key(@keys){
    if($key =~ /vol\d/ ){$nstates++};   
}
$name_value{'nstates'} = max 1,$nstates;

# set time per image
if( defined($name_value{'time_per_image'}) ){
} else {
    if( defined($name_value{'trs'}) ){
        $name_value{'time_per_image'} = 100*$name_value{'nstates'}*($name_value{'nspace'}/1000);
    }else{
        $name_value{'time_per_image'} = 90*$name_value{'nstates'}*($name_value{'nspace'}/1000);
    }
}

# do the refinement
my$converged = 0;
my$round;
my$round_stop;
if( defined($name_value{'startit'}) ){
    $round = $name_value{'startit'}-1; 
} else {
    $round = 0;
}
while( $round < $name_value{'maxits'} ){
    
    # make refinement folder
    $round++;
    $folder = 'prime_round_'.$round;
    mkdir $folder;
    
    # parallell execution
    exec_prime_para();
    
    # determine when the jobs have finished
    system("rm -f JOB_FINISHED_*");
    sleep(10) while( njobs_finished() < $name_value{'npart'} );
    system("rm -f JOB_FINISHED_* OUT* prime_script*");
    
    # merge alignment docs
    if( defined($name_value{'oritab'}) ){
        delete $name_value{'oritab'};  
    }
    $name_value{'oritab'} = 'merged.txt';
    system("$name_value{'merge_algndocs'} fbody=algndoc_ nptcls=$name_value{'nptcls'} ndocs=$name_value{'npart'} outfile=$name_value{'oritab'}");
    system("mv algndoc_* ./$folder");
    system("cp merged.* ./$folder");
    
    # assemble volumes
    exec_assemble_volumes();
    
    # automask
    if( defined($name_value{'mw'}) ){ 
        exec_automask();
    }
    
    my $conv = exec_check_conv();
    if( $round > 1 ){
        if( $conv == 2 ){
            die "**** DISTR_SIMPLE_PRIME NORMAL STOP ****\n";
        }elsif ( $conv == 1 ){
            if( !defined($name_value{'find'}) ){
                die "Want 2 do dynamic low-pass update, but find is not defined, weird!\n";   
            }
            if( defined($name_value{'deterministic'}) ){
                delete $name_value{'deterministic'};
                $name_value{'deterministic'} = 'yes';
            }
            $name_value{'find'} = $name_value{'find'}+$name_value{'fstep'};
            if( $name_value{'find'} > $name_value{'fmax'} ){
                $name_value{'find'} = $name_value{'fmax'}      
            }
        }
    }
}

sub exec_prime_para{
    my$i;
    my$vol;
    my@submitted;
    my$qsubout;
    my@args = parse_args($name_value{'prime'});
    my$prime_cmd_string = $name_value{'prime'}.' ';
    $i = 1;
    $vol = 'vol'.$i;
    while( defined($name_value{$vol}) ){
        $prime_cmd_string = add2string($vol, $prime_cmd_string);
        $i++;
        $vol = 'vol'.$i;
    }
    foreach (@args) {
        $prime_cmd_string = add2string($_, $prime_cmd_string);
    }   
    $prime_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $prime_cmd_string, "\n";
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
        open(FHANDLE, ">prime_script_$i");
        if ($ENV{SIMPLESYS} eq 'HPC'){
            print FHANDLE "#!/bin/bash
#SBATCH -n $name_value{'nthr'}
##SBATCH -N 1-1
#SBATCH -J distr_simple
#SBATCH -p grant
#SBATCH -A g2014a90
#SBATCH -t 4-00:00:00
cd $name_value{'execdir'} 
$prime_cmd_string fromp=$start top=$stop part=$i outfile=algndoc_$i.txt > OUT$i\nexit\n";
        }elsif ($ENV{SIMPLESYS} eq 'DELPHI' ){
            print FHANDLE "#!/bin/bash
#\$ \-N d_simple
#\$ \-cwd
#\$ \-pe multi $name_value{'nthr'}
#\$ \-V
#\$ \-o outfile.\$PBS_JOBID
#\$ \-e errfile.\$PBS_JOBID
cd $name_value{'execdir'} 
$prime_cmd_string fromp=$start top=$stop part=$i outfile=algndoc_$i.txt > OUT$i\nexit\n";
        } else {
            die "Unrecognized system!\n";
        }
        close(FHANDLE);
        chmod 0777, 'prime_script_'.$i;  
        if( $ENV{SIMPLESYS} eq 'LOCAL' ){
            system("nohup ./prime_script_$i &\n");
        }elsif( $ENV{SIMPLESYS} eq 'HPC' ){
            my$qsubcmd = "sbatch ./prime_script_$i";
            $qsubout = `$qsubcmd`;
            if( $qsubout =~ /Submitted batch job/ ){
                $submitted[$i] = 1;   
            }else{    
                $submitted[$i] = 0; 
            }
        }elsif( $ENV{SIMPLESYS} eq 'DELPHI' ){
            my$qsubcmd = "qsub ./prime_script_$i";
            $qsubout = `$qsubcmd`;
            if( $qsubout =~ /has been submitted/ ){
                $submitted[$i] = 1;   
            }else{    
                $submitted[$i] = 0; 
            }
        }else{
            die "Unrecognized system!\n";
        }
    }
    if( $ENV{SIMPLESYS} eq 'HPC' ){
        # resubmission if failed
        for(my $i=1; $i<=$name_value{'npart'}; $i++){
            if( $submitted[$i] == 0 ){
                my$qsubcmd = "sbatch ./prime_script_$i";
                $qsubout = `$qsubcmd`;
                if( $qsubout !~ /Submitted batch job/ ){
                    chomp($qsubout);  
                    print "ERROR, $qsubout\n";
                }
            }
        }
    } elsif( $ENV{SIMPLESYS} eq 'DELPHI' ){
        # resubmission if failed
        for(my $i=1; $i<=$name_value{'npart'}; $i++){
            if( $submitted[$i] == 0 ){
                my$qsubcmd = "qsub ./prime_script_$i";
                $qsubout = `$qsubcmd`;
                if( $qsubout !~ /has been submitted/ ){
                    chomp($qsubout);  
                    print "ERROR, $qsubout\n";
                }
            }
        }
    }
}

sub njobs_finished{
    my$nfini = 0;
    my@jobs;
    @jobs = <JOB_FINISHED_*>;
    return scalar(@jobs);     
}

sub exec_assemble_volumes{
    my$i;
    my$vol;
    my@refvols;
    my@recvols;
    my@args = parse_args($name_value{'volassemble'});
    my$assemble_cmd_string = $name_value{'volassemble'}.' ';
    foreach (@args) {
        $assemble_cmd_string = add2string($_, $assemble_cmd_string);
    }
    $assemble_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $assemble_cmd_string, "\n";
    system($assemble_cmd_string);     
    my@kernels = <kernel_state*bin>;
    my@recbins = <recvol_state*bin>;
    if( scalar(@kernels) != scalar(@recbins) ){
        die "Number of kernels not equal to number of recbins";   
    }
    foreach my$i (0 .. $#kernels){
        system("rm -f $kernels[$i] $recbins[$i]");
    }
    @recvols = <recvol_state*spi>;
    if( scalar(@recvols) != $name_value{'nstates'} ){
        die "Number of reference volumes not equal to number of states";   
    } else {
        @refvols = @recvols;
    }
    foreach $i ( 1 .. $name_value{'nstates'} ){
        $vol = 'vol'.$i;
        if( defined($name_value{$vol}) ){
            delete $name_value{$vol};
        }
        $name_value{$vol} = $refvols[$i-1];
        system("cp $name_value{$vol} ./$folder");
    }
}

sub exec_automask{
    my$i;
    my$vol;
    my@args = parse_args($name_value{'automask'});
    my$automask_cmd_string = $name_value{'automask'}.' ';
    $i = 1;
    $vol = 'vol'.$i;
    while( defined($name_value{$vol}) ){
        $automask_cmd_string = add2string($vol, $automask_cmd_string);
        $i++;
        $vol = 'vol'.$i;
    }
    foreach (@args) {
        $automask_cmd_string = add2string($_, $automask_cmd_string);
    }
    $automask_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $automask_cmd_string, "\n";
    system($automask_cmd_string);
    my$fbody      = $name_value{'vol1'};
    $fbody        =~ s/\d+\.spi//;
    my@maskedvols = <$fbody*msk.spi>;
    foreach $i ( 1 .. $name_value{'nstates'} ){
        $vol = 'vol'.$i;
        system("cp $maskedvols[$i-1] ./$folder");
        system("mv $maskedvols[$i-1] $name_value{$vol}");
    }
}

sub exec_check_conv{
    my$status;
    my@args = parse_args($name_value{'check_conv'});
    my$check_conv_cmd_string = $name_value{'check_conv'}.' ';
    foreach (@args) {
        $check_conv_cmd_string = add2string($_, $check_conv_cmd_string);
    }
    $check_conv_cmd_string =~ s/\s+$//; # removes trailing whitespace
    print $check_conv_cmd_string, "\n";
    my$conv_info = `$check_conv_cmd_string`;
    my@conv_info_lines = split(/\n/,$conv_info);
    $status = 0;
    foreach my$line (@conv_info_lines){
        print $line, "\n";
        if( $line =~ /\>\>\> UPDATE LOW-PASS LIMIT\: \.YES\./ ){ $status = 1 };
        if( $line =~ /\>\>\> CONVERGED\: \.YES\./ )            { $status = 2 };
    }
    return $status;
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

sub parse_args{
    my$prg = shift; # get program name with absolute path
    # make command line instructions
    my$instructions = `$prg`;
    # parse the instructions
    my@splitted = split(/\s/,$instructions);
    my@args;
    foreach (@splitted) {
        if( $_ =~ /(\w+\d*)\=/ ){
            if( $1 !~ /vol\d+/ ){
                push(@args,$1);
            }  
        }
    }
    return @args;
}

sub get_find{
    my$res = shift;
    return int((($name_value{'box'}-1)*$name_value{'smpd'})/$res)
}

sub get_lp{
    my$find = shift;
    return (($name_value{'box'}-1)*$name_value{'smpd'})/$find
}
