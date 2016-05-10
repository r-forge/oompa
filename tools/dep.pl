#! perl -w

use strict;
use warnings;
use Cwd;
my $location = getcwd;

my $oompafiles = "oompafiles.txt";
my @packList = ();
open(OOMPA, "<$oompafiles") or die "Unable to open '$oompafiles': $!\n";
while (my $line = <OOMPA>) {
    chomp($line);
    push @packList, $line;
}
close(OOMPA);

my $desc = 
chdir "../.." or die "Unable to move up two directories: $!\n";
my $home = getcwd;
my %dephash = ();

print STDERR "Working from $home\nLooking for:\n";
foreach my $p (@packList) {
    print STDERR "\t$p\n";
}

foreach my $pdir (@packList) {
    print STDERR "ready for '$pdir'\n";
    chdir $pdir or die "Unable to change to '$pdir': $!\n";

    open(DESC, "<DESCRIPTION") or die "Unable to open 'DESCRIPTION' in '$pdir': $!\n";
    while (my $line = <DESC>) {
	if ($line =~ /^Depends/ or 
	    $line =~ /^Imports/ or
	    $line =~ /^Suggests/) {
	    chomp($line);
	    my @stuff = split /: /, $line;
	    my @packages = split ',', $stuff[1];
	    foreach my $p (@packages) {
		if ($p =~ /(.*?) \((.*?)\)/) {
		    $p = $1;
		}
		$p =~ s/^\s+//;
		$p =~ s/\s+$//;
		$dephash{$p} = 1;
	    }
	}
    }
    close(DESC);
    chdir $home or die "Thomas Wolfe was right: $!\n";
}
chdir $location;

my $words = join("\", \"", keys(%dephash));
open(DEP, ">dependencies.R") or die "Unable to create 'dependencies.txt': $!\n";
print DEP "source(\"https://bioconductor.org/biocLite.R\")\nbiocLite()\n\n";
print DEP "install.packages(c(\"$words\"))\n";
close(DEP);

exit;
__END__
