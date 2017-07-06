#! perl -w

use strict;
use warnings;
use File::Copy;
use Cwd;
my $home = getcwd;

# get command line arguments
my $oompafiles = shift || "oompafiles.txt";
my $tools = shift || "tools.cfg";
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
    my $webDir = $pdir;
    $webDir =~ s/pkg/www/; # web root for the current project
    unless (-e $webDir) {
	print STDERR "Cannot find '$webDir'. Skipping '$pack'.\n";
	next;
    } else {
	$webDir = "$home/$webDir/packages"; # where to put the auto-generated web pages
	printf STDERR " Target location is '$webDir'.\n";
    }
    unless (-e $webDir) {
	mkdir($webDir);
    }
    my $infoDir = "$webDir/$pack"; # place to put package stuff, like R CHECK logs and vignettes
    unless (-e $infoDir) {
	mkdir($infoDir);
    }
    open(TGT, ">$webDir/$pack.html") or die "Unable to create '$pack.html': $!\n";
    my $header = <<"ENDHEAD";
<htmL>
<!-- autocreated by harvestHTML.pl -->
<title>Package $pack</title><head></head>
<body><h1>Package $pack</h1>
<table>
ENDHEAD
;
    print TGT $header;
    $pdir = "$pdir/$pack.Rcheck";
    chdir $pdir or die "Unable to change to '$pdir': $!\n";
    open(DESC, "<$pack/DESCRIPTION") or die "Unable to open 'DESCRIPTION' file: $!\n";
    my $tag = '';
    my $value = '';
    while (my $line = <DESC>) {
	chomp $line;
	if ($line =~ /^([a-zA-Z]+)\: (.*)$/) {
	    print TGT "<tr><td>$tag</td>\n\t<td>$value</td></tr>\n" if $tag;
	    $tag = $1;
	    $value = $2
	} else {
	    $value .= $line;
	    $value =~ s/\s+/ /g;
    }
    }
    print TGT "<tr><td>$tag</td>\n\t<td>$value</td>\n</tr>\n" if $tag;
    close DESC;
    my $man = "$pack-manual.pdf";
    if (-e $man) {
	copy($man, $infoDir) or die "Copy failed for '$man': $!\n";
	print TGT "<tr><td>User Manual</td>",
	"<td><A HREF=\"$pack/$man\">$man</A></td></tr>";
    } else {
	print STDERR "Missing manual: $man\n";
    }
    my $check = "00check.log";
    if (-e $check) {
	copy($check, $infoDir) or die "Copy failed for '$check': $!\n";
	print TGT "<tr><td>R CHECK</td>",
	"<td><A HREF=\"$pack/$check\">$check</A></td></tr>";
    } else {
	print STDERR "Missing manual: $check\n";
    }
    chdir "$pack/doc";
    my @vigs = glob "*.pdf *.html";
    if ($#vigs) { # shortcut, since 0 means just index.html
	print TGT "<tr><td>Vignettes</td>\n  <td>";
	foreach my $f (@vigs) {
	    next if ($f eq 'index.html');
	    copy($f, $infoDir) or die "Copy failed for vignette '$f': $!\n";
	    print TGT "<A HREF=\"$pack/$f\">$f</A>, ";
	}
	print TGT "</td></ttr>\n";
    }
    print TGT "\n</table></body></html>\n";
    close(TGT);
    chdir $home or die "Unable to return to '$home': $!\n";;
}

exit;
__END__
