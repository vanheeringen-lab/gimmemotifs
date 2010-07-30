use ExtUtils::MakeMaker;
use File::Copy;

if ($^V lt v5.6) {

   die("Algorithm::Cluster is only supported on Perl 5.6.0 and later.\n" .
   "Your perl version is $].\n");
}

copy("perl/MANIFEST.perl","MANIFEST");

WriteMakefile(
	NAME         => 'Algorithm::Cluster',
	VERSION_FROM => 'perl/Cluster.pm',
        AUTHOR       => 'John Nolan and Michiel de Hoon (mdehoon "AT" gsc.riken.jp)',
        ABSTRACT     => 'Perl interface to the C Clustering Library',
	DIR          => [ 
		'src', 
		'perl',
	],
);
