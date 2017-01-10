#! perl -w

use strict;
use warnings;
use Cwd;
my $home = getcwd;

# get command line arguments
my $oompafiles = shift || "oompafiles.txt";
my $tools = shift || "tools.cfg";
my $startdir = shift || '..\..'; # assumes tools lives two levels below RForge

# TODO: allow cfg file to be passed on the command line ?
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

my $profile =<<EOP
r <- getOption("repos")
r["CRAN"] <- "https://cran.mtu.edu"
r["R-Forge"] <- "http://R-Forge.R-project.org"
options(repos = r)
rm(r)
EOP
    ;
my $profloc = $ENV{R_PROFILE_USER} || "$ENV{HOME}/.Rprofile";
my $rmProfile = 0;
if (!(-e $profloc)) {
    my $ok = open(PROF, ">$profloc");
    if ($ok) {
	print PROF $profile;
	close PROF;
	$rmProfile = 1;
    } else {
	print STDERR "Unable to create $profloc: $!\n";
    }
}
END {
    unlink($profloc) if ($rmProfile);
}


# TODO: allow file to be passed on the command line ?
# Assume file consists of one package per line.
my @packList = ();
open(OOMPA, "<$oompafiles") or die "Unable to open '$oompafiles': $!\n";
while (my $line = <OOMPA>) {
    chomp($line);
    push @packList, $line;
}


#####################################################
if ($startdir) {
    chdir($startdir) or die "Unable to change directory to '$startdir': $!\n";
    $home = cwd;
}
print STDERR "\nStarting at '$home'\n\n";

# remember the PATH where we started
my $path = $ENV{PATH};
print STDERR "start path:\n$path\n\n";

# make sure there is a directory for source tarballs
my $buildDir = "Build-$Ver.$s";
if (-e $buildDir) {
    die("'$buildDir' exists and is not a directory") unless (-d $buildDir);
} else {
    mkdir($buildDir);
}

# make sure there is a directory for binary packages
my $archDir = "Binary-$Ver.$s";
$archDir = "$archDir-$arch" if $arch; 
if (-e $archDir) {
    die("'$archDir' exists and is not a directory") unless (-d $archDir);
} else {
    mkdir($archDir);
}

#####################################################
# set the path explicitly. Note that the path to the main
# R executable varies by architecture.
my $Rpath = "$toolhash{rdir}\\R-$Ver.$s\\bin\\$arch";

my @paths = ($Rpath,
	     "C:\\Rtools\\$Ver\\bin",
	     "C:\\Rtools\\$Ver\\MinGW\\bin",
	     "C:\\Rtools\\$Ver\\gcc-4.6.3\\bin",
	     $toolhash{miktex},
	     "C:\\Windows\\system32",
	     "C:\\Windows",
	     $toolhash{perl}."\\bin",
	     $toolhash{perl}."\\site\\bin",
	     $toolhash{qpdf},
	     $toolhash{gs},
	     $toolhash{pandoc},
);

$ENV{PATH} = join(";", @paths);
$ENV{nodosfilewarning} = '1';
print STDERR "working path:\n$ENV{PATH}\n\n";

# process each package
foreach my $packname (@packList) {
    my @parts = split /\//, $packname;
    my $pack = pop @parts;
    my $pdir = join("/", @parts) || $home;
    chdir $pdir or die "Unable to change to '$pdir': $!\n";
# build the source tarball
    my @cargs = ('R', 'CMD', 'build', $pack);
    my $cmd = join(' ', @cargs);
    print STDERR "Running $cmd...\n";
    my $check = system(@cargs);
    last unless $check == 0;
    my @files = glob("$pack*tar.gz");
    die("got the wrong number of files") if ($#files);
    my $tarball = $files[0];
# build the binary version for the current architecture
    @cargs = ('R', 'CMD', 'INSTALL', '--build', $tarball);
    $cmd = join(' ', @cargs);
    print STDERR "Running $cmd...\n";
    $check = system(@cargs);
    last unless $check == 0;
# move it to the architecture binary directory
    @files = glob("$pack*zip");
    die("got the wrong number of files") if ($#files);
    my $file = $files[0];
    rename $file, "$home/$archDir/$file";
# check the package from the tarball
    @cargs = ('R', 'CMD', 'check', '--as-cran', $tarball);
    $cmd = join(' ', @cargs);
    print STDERR "\nRunning $cmd...\n";
    $check = system(@cargs);
    last unless $check == 0;
# install the package in the current R version
    @cargs = ('R', 'CMD', 'INSTALL', $tarball);
    $cmd = join(' ', @cargs);
    print STDERR "Running $cmd...\n";
    $check = system(@cargs);
    last unless $check == 0;
# move the tarball to the build directory
    rename $tarball, "$home/$buildDir/$tarball";
    chdir $home or die "Thomas Wolfe was right: $!\n";
}

exit;
__END__
