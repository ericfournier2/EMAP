#!/usr/bin/perl

use CGI;
use Archive::Zip qw( :ERROR_CODES :CONSTANTS );
use List::MoreUtils 'first_index';
use Cwd;
use File::Basename;
use POSIX;

sub doError {
    my $errMessage = shift;
    print "<p>Your request could not be processed because the following error occured:</p>";
    print '<p style="margin-left:10%">' . $errMessage . "</p>";
    print "<p>Please make sure you have read the instruction on the request page. If you find yourself unable to correct your request, contact the system administrator.</p></body>";
    exit 0;
}

my $q = new CGI;

print CGI::header;
print "<html><head><title>EMAP web interface</title></head><body>";

# Read in form parameters.
my $email = $q->param("email");
my $expName = $q->param("expName");
my $password = $q->param("passwd");
my $platform = $q->param("platform");
my $file = $q->param("file");
my $refCond = $q->param("refCond");

my $doEpi = $q->param("doEpi");
my $epiFC = $q->param("epiFC");
my $epiP = $q->param("epiP");
my $epiFDR = $q->param("epiFDR") eq "Yes";

my $doTrans = $q->param("doTrans");
my $transFC = $q->param("transFC");
my $transP = $q->param("transP");
my $transFDR = $q->param("transFDR") eq "Yes";

#my $password = "Epigenetics";
#my $expName = "Eric";
#my $refCond = "Embryo";
#my $doEpi = 1;
#my $doTrans = 0;

# Make sure this is an authorized user.
if($password ne "Epigenetics") {
    doError("The password you provided is incorrect.");
}

# Generate a directory to hold the uploaded files.
my $save_dir = "/home/EMAP-Upload/" . $expName . "/";
umask 000;
mkdir($save_dir, 0775);

# Write the uploaded zip archive to disk.
open(OUTFILE, ">$save_dir/file.zip") || die "can't create file: $!";
while (read($file, $buffer, 1024)) { print OUTFILE $buffer; }
close(OUTFILE);
chmod (0666, "$save_dir/file.zip");

# Unzip the archive.
my $file_patch = $save_dir . "/file.zip";
my $zip = Archive::Zip->new($file_patch);
my @files = $zip->members();

foreach (@files) {
    my $unzip_file = $_->fileName;
    my $unzip_target = $save_dir . "/" . $unzip_file;

    if ($zip->extractMember($unzip_file, $unzip_target) != 0) {
      doError("Extraction of $file failed: $!");
    }
}

# Is there one and only one epi.target file?
my $epiTargetFile = "";
if($doEpi) {
    my $epiGlob = $save_dir . "*.epi.target*";
    @epi_target_files = glob $epiGlob;
    if(scalar(@epi_target_files) == 0) {
        doError("No epi.target file found in archive.");
    } elsif(scalar(@epi_target_files) > 1) {
        doError("Archive must contain exactly one epi.target file.");
    }
    
    $epiTargetFile = basename($epi_target_files[0]);
}

# Is there one and only one trans.target file?
my $transTargetFile = "";
if($doTrans) {
    my $transGlob = $save_dir . "*.trans.target*";
    @trans_target_files = glob $transGlob;
    if(scalar(@trans_target_files) == 0) {
        doError("No trans.target file found in archive.");
    } elsif(scalar(@trans_target_files) > 1) {
        doError("Archive must contain exactly one trans.target file.");
    }
    
    $transTargetFile = basename($trans_target_files[0]);
}

# Validate the contents of the target files.
my $allTargetGlob = $save_dir . "*.target*";
@target_files = glob $allTargetGlob;
foreach my $targetFile (@target_files) {
    # Are all the right columns present?
    open(TARGET, "<", $targetFile);
    my $header = <TARGET>;
    $header =~ s/\r\n/\n/;
    chomp($header);
    my @colnames = split("\t", $header);
    my $filenameColumn = first_index { /Filename/ } @colnames;
    my $cy5Column = first_index { /Cy3/ } @colnames;
    my $cy3Column = first_index { /Cy5/ } @colnames;

    if(!defined($filenameColumn) || !defined($cy5Column) || !defined($cy3Column)) {
        doError($targetFile . " does not have the correct column headers.");
    }

    # Validate individual lines in the target file.
    while(my $line=<TARGET>) {
        $line =~ s/\r\n/\n/;
        chomp($line);
        my @data=split("\t", $line);

        # Are all files named in the .target file present within the archive?
        my $filename=$data[$filenameColumn];
        if(!(-e $save_dir . "/" . $filename)) {
            doError("File " . $filename . " is referenced in " . $targetFile . ", but absent from the archive.");
        }
        
        # Is the reference condition one of either the Cy3 or Cy5 channel?
        if(($data[$cy5Column] ne $refCond) && ($data[$cy3Column] ne $refCond)) {
            doError("File " . $filename . " does not have any channel with the reference condition.");
        }
    }
    close(TARGET);
}
# "Simon/Epigenetic/Combined" "Design-Combined.target" "Raw/Simon"
# Build command line for launching EMAP
chdir("/home/efournier/Projects/Epigenetics/EMAP");

my $cmdLine = 'bash Scripts/FullAnalysis.sh ';

if($doEpi) {
    my $epiResults = $expName . "/Epigenetic";
    $cmdLine = $cmdLine . '"' . $epiResults . '" "' . $epiTargetFile . '" "' . $save_dir . '" ';
} else {
    $cmdLine = $cmdLine . '"" "" "" ';
}

if($doTrans) {
    my $transResults = $expName . "/Transcriptomic";
    $cmdLine = $cmdLine . '"' . $transResults . '" "' . $transTargetFile . '" "' . $save_dir . '" ';
} else {
    $cmdLine = $cmdLine . '"" "" "" ';
}

$cmdLine = $cmdLine . $refCond . " ";
if($doTrans && $doEpi) {
    my $combinedResults = $expName . "/Combined";
    $cmdLine = $cmdLine . $combinedResults . " ";
} else {
    $cmdLine = $cmdLine . ' "" ';
}

if($doEpi) {
    $cmdLine = $cmdLine . $epiFC . " " . $epiP . " ";
    if($epiFDR) {
        $cmdLine = $cmdLine . "TRUE ";
    } else {
        $cmdLine = $cmdLine . "FALSE ";
    }
} else {
    $cmdLine = $cmdLine . "1.5 0.05 FALSE ";
}

if($doTrans) {
    $cmdLine = $cmdLine . $transFC . " " . $transP . " ";
    if($transFDR) {
        $cmdLine = $cmdLine . "TRUE ";
    } else {
        $cmdLine = $cmdLine . "FALSE ";
    }
} else {
    $cmdLine = $cmdLine . "1.5 0.05 FALSE ";
}

$cmdLine = $cmdLine . $platform . " " . $email;

# Launch the analysis process.
if(fork()) {
    # We're in the daemon's thread. Set the session ID, close file handles, 
    # and otherwise prepare for running EMAP.
    setsid();
    fork and exit;
    
    close STDIN;
    close STDOUT;
    close STDERR;
    
    # Run EMAP.
    system($cmdLine);
} else {
    print "<p>Your request has been processed successfully.</p>";
    print "<p>You should be notified by e-mail when the results become available: the process usually takes a few hours.</p>";
    print "<p>If you do not receive a notification e-mail, please check your spam filter.</p></body>";
}
