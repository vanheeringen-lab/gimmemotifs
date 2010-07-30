use ExtUtils::MakeMaker;
use Config;
use POSIX;

# Check if this is a 64 bit machine
my $CCFLAGS = '';
my $machine = (POSIX::uname())[4];
if ($machine =~ /64/) {
        $CCFLAGS = '-fPIC';
}

WriteMakefile(
	NAME         => 'libcluster',
	VERSION_FROM => '../perl/Cluster.pm',
	SKIP         => [qw(all static static_lib dynamic dynamic_lib)],
	C	     => ['cluster.c'],
	H	     => ['cluster.h'],
	OBJECT       => 'cluster.o',
	MYEXTLIB     => 'libcluster$(LIB_EXT)',
	CCFLAGS      => $CCFLAGS,
	clean        => {'FILES' => 'libcluster$(LIBEEXT) $(OBJECT)'},
);


sub MY::top_targets {
'
all :: static

pure_all :: static

static :: $(OBJECT) libcluster$(LIB_EXT)

cluster.o :
	$(CC) -c $(INC) $(DEFINE_VERSION) $(XS_DEFINE_VERSION) $(CCFLAGS) -o cluster.o cluster.c

libcluster$(LIB_EXT): $(O_FILES)
	$(AR) cru libcluster$(LIB_EXT) $(OBJECT)
	$(RANLIB) libcluster$(LIB_EXT)
';
}
