#! perl -w

use strict;
use warnings;
use File::Copy;
use Cwd;
my $home = getcwd;

# get command line arguments
my $oompafiles = shift || "oompafiles.txt";
my $tools = shift || "tools.cfg";

# Assume file consists of lines of the form "KEY = VALUE"
open(TOOL, "<$tools") or die "Unable to open '$tools': $!\n";
my %toolhash = ();
while (my $line = <TOOL>) {
    chomp $line;
    my ($key, $value) = split(" = ", $line);
    $toolhash{$key} = $value;
} 
die "'$tools' must define the 'rdir' location where R is installed" unless $toolhash{rdir};
die "'$tools' must define the 'perl' location" unless $toolhash{perl};
die "'$tools' must define the 'miktex' location" unless $toolhash{miktex};
die "'$tools' must define the 'qpdf' location" unless $toolhash{qpdf};
die "'$tools' must define the 'gs' location" unless $toolhash{gs};

my $Ver  = $toolhash{major};
my $s    = $toolhash{point};
my $arch = $toolhash{arch};

# Assume file consists of one package per line.
my @packList = ();
open(OOMPA, "<$oompafiles") or die "Unable to open '$oompafiles': $!\n";
while (my $line = <OOMPA>) {
    chomp($line);
    push @packList, $line;
}

my $manDir = "$home/oompa/www/manuals";
if (-e $manDir) {
    die("'$manDir' exists and is not a directory") unless (-d $manDir);
} else {
    mkdir($manDir);
}

foreach my $packname (@packList) {
    last if $packname =~ /^#/;
    my @parts = split /\//, $packname;
    my $pack = pop @parts;
    my $pdir = join("/", @parts);
    $pdir = "$pdir/$pack.Rcheck";
    chdir $pdir or die "Unable to change to '$pdir': $!\n";
    my $man = "$pack-manual.pdf";
    if (-e $man) {
	copy($man, $manDir) or die "Copy failed for '$man': $!\n";
    } else {
	print STDERR "Missing manual: $man\n";
    }
    chdir $home or die "Unable to return to '$home': $!\n";;
}

exit;
__END__
