#! perl -w

use strict;
use warnings;
use File::Path;
use Cwd;
my $home = getcwd;

# get command line arguments
my $oompafiles = shift or die "You must supply a list of package directories!\n";
my $tools = shift or die "You must supplay a tool configuration file\n";
my $startdir = shift || '..\..'; # assumes tools lives two levels below RForge

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

if ($startdir) {
    chdir($startdir) or die "Unable to change directory to '$startdir': $!\n";
    $home = cwd;
}
print STDERR "Home is '$home'.\n";


foreach my $packname (@packList) {
    last if $packname =~ /^#/;
    my @parts = split /\//, $packname;
    my $pack = pop @parts;
    my $pdir = join("/", @parts);
    chdir $pdir or die "Unable to change to '$pdir': $!\n";
    my $checkDir = "$pdir/$pack.Rcheck";
    if (-e $checkDir) {
	my $n = rmtree($checkDir);
	print STDERR "Removed $n files from '$checkDir'\n";
    }
    chdir $home or die "Unable to return to '$home': $!\n";;
}
exit;
__END__
