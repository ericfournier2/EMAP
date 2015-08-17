my $window = $ARGV[0];

my $lastPos = 1;
my $nextStop = $lastPos + $window;
my $curChr = "";
my @valueBuffer = ();
while(my $line = <STDIN>) {
     my @data = split("\t", $line);
    
    if(($data[1] > $nextStop) || ($data[0] ne $curChr)) {
        my $sum = 0;
        foreach my $val (@valueBuffer) {
            $sum += $val;
        }
        
        my $mean = 0;
        if(scalar(@valueBuffer) != 0) {
            $mean = $sum / scalar(@valueBuffer);
            print $curChr . "\t" . $lastPos . "\t" . $nextStop . "\t" . $mean . "\n";
        }
                
        if($data[0] ne $curChr) {
            $curChr = $data[0];
            $lastPos = $data[1];
        } else {
            $lastPos = $data[1];
        }
        $nextStop = $lastPos + $window;
        @valueBuffer = ($data[3]);
    } else {
        push(@valueBuffer, $data[3]);
    }
}

my $sum = 0;
foreach my $val (@valueBuffer) {
    $sum += $val;
}

my $mean = 0;
if(scalar(@valueBuffer) != 0) {
    $mean = $sum / scalar(@valueBuffer);
    print $curChr . "\t" . $lastPos . "\t" . $nextStop . "\t" . $mean . "\n";
}