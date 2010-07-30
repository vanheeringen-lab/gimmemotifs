use ExtUtils::MakeMaker;
use Config;

WriteMakefile(
	NAME		=> 'Algorithm::Cluster',
	AUTHOR		=> 'John Nolan and Michiel de Hoon',
	ABSTRACT	=> 'Perl interface to the C Clustering Library',
	VERSION_FROM	=> 'Cluster.pm',
	PM		=> 
		{
			'Cluster.pm' =>	'$(INST_LIBDIR)/Cluster.pm',
			'Record.pm' =>	'$(INST_LIBDIR)/Cluster/Record.pm',
		},
	LIBS		=> '-lm',
	INC		=> '-I../src',
	MYEXTLIB	=> '../src/libcluster$(LIB_EXT)',
);

