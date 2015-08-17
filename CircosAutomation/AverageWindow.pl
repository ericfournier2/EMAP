my $window = $ARGV[0];
my $minimumProbes = $ARGV[1];

my $lastPos = 1;
my $nextStop = $lastPos + $window;
my $curChr = "";
my @valueBuffer = ();
while(my $line = <STDIN>) {
    my @data = split("\t", $line);
    
     if($curChr eq "") {
        $curChr = $data[0];
     }
    
    if(($data[1] > $nextStop) || ($data[0] ne $curChr)) {
        my $sum = 0;
        foreach my $val (@valueBuffer) {
            $sum += $val;
        }
        
        my $mean = 0;
        if(scalar(@valueBuffer) != 0) {
            $mean = $sum / scalar(@valueBuffer);
        }
        
        if(scalar(@valueBuffer) >= $minimumProbes) {
            print $curChr . "\t" . $lastPos . "\t" . $nextStop . "\t" . $mean . "\n";
        }
        
        if($data[0] ne $curChr) {
            $curChr = $data[0];
            $lastPos = 1;
        } else {
            $lastPos = $nextStop;
        }
        $nextStop = $lastPos + $window;
        @valueBuffer = ();
    } else {
        push(@valueBuffer, $data[3]);
    }
}