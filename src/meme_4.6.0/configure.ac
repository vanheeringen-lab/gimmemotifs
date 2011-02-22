## $Id: configure.ac 5293 2011-01-11 00:54:31Z cegrant $
##
## Process this file with autoconf to produce a configure script.
##
##

AC_PREREQ(2.53)
AC_INIT([meme], [4.6.0])
AM_CONFIG_HEADER([config.h])
AC_CONFIG_SRCDIR(src/alphtype.c)
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([1.9 tar-pax])
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AM_PROG_CC_C_O
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PATH_PROG(LD, ld)
AC_PATH_PROG(SED, sed)
AC_PATH_PROG(MPICC, mpicc mpcc)
AC_PATH_PROG(MAIL, mail)
AC_PATH_PROG(SENDMAIL, sendmail, [/usr/lib/sendmail], [$PATH:/usr/lib:/usr/sbin])
AC_PATH_PROG(PERL, perl)
AC_PATH_PROG(PYTHON, python)
AC_PATH_PROG(CONVERT, convert)
AC_PATH_PROG(GHOSTSCRIPT, gs)
AC_PATH_PROG(DOT, dot)
AC_PATH_PROG(XML2_CONFIG, xml2-config)
AC_PATH_PROG(XSLT_CONFIG, xslt-config)
AC_PATH_PROG(ANT, ant)

AC_DEFINE_UNQUOTED(HOSTNAME, "`hostname`", [Build host name])

# Checks for libraries.
AC_CHECK_LIB([m], 
             [isnan, sqrt, log, exp, fabs, pow, floor, ceil], 
             [], 
             [AC_MSG_WARN([Math library not found])])

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([arpa/inet.h ctype.h fcntl.h limits.h netdb.h netinet/in.h stdlib.h string.h strings.h sys/socket.h sys/time.h unistd.h ])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_CHECK_FUNCS([floor getcwd gethostbyname isascii memset pow rint socket sqrt strchr strcspn strdup strspn strstr])
AC_CHECK_FUNCS([strlcpy], [have_strlcpy=yes], [have_strlcpy=no])
AM_CONDITIONAL(WANT_STRLCPY, test ["$have_strlcpy" = no])

# Check for compiler-specific settings
LARGE_FILE_SUPPORT=`getconf LFS_CFLAGS`
if test "${CC}" = "gcc"; then
  ac_cc_opt="-O3"
  ac_cc_warn="-Wall -Wno-unused"
  ac_cc_debug="-ggdb"
  CFLAGS="$ac_cc_warn -DUNIX -D__USE_FIXED_PROTOTYPES__ $LARGE_FILE_SUPPORT"
else
  ac_cc_opt="-O"
  ac_cc_warn=""
  ac_cc_debug="-g"
  CFLAGS="$ac_cc_warn -DUNIX -D__USE_FIXED_PROTOTYPES__ $LARGE_FILE_SUPPORT"
fi

# Set up default libxml2 settings
# check for installed copy of libxml2
build_libxml2="yes"
if test -n "$XML2_CONFIG"; then
  # Require at least version 2.6
  required_libxml2_version="2006000"
  LIBXML2_VERSION="`$XML2_CONFIG --version | \
    awk 'BEGIN { FS = "."; } { printf "%d", ($1 * 1000 + $2) * 1000 + $3;}'`"
  if test -n "$LIBXML2_VERSION" && \
     test "$LIBXML2_VERSION"  -ge "$required_libxml2_version"; then 
    build_libxml2="no"
    LIBXML2_LIBS="`$XML2_CONFIG --libs`"
    LIBXML2_CFLAGS="`$XML2_CONFIG --cflags`"
  fi
fi
if test -z "$LIBXML2_LIBS"; then
  AC_MSG_WARN([libxml2 library not found])
  LIBXML2_LIBS="\${top_builddir}/src/libxml2/libxml2.la"
  LIBXML2_CFLAGS="-I\${top_srcdir}/src/libxml2/include"
fi

# Set up default libxslt settings
# check for installed copy of libxslt
build_libxslt="yes"
if test -n "$XSLT_CONFIG"; then
  build_libxslt="no"
  LIBXSLT_LIBS="`$XSLT_CONFIG --libs`"
  LIBXSLT_CFLAGS="`$XSLT_CONFIG --cflags`"
fi
if test -z "$LIBXSLT_LIBS"; then
  AC_MSG_WARN([libxslt library not found])
  LIBXSLT_LIBS="\${top_builddir}/src/libxslt/libxslt.la"
  LIBXSLT_CFLAGS="-I\${top_srcdir}/src/"
fi

# who is installing 
USER=`whoami`

# Check for platform-specific settings
case $host in
   *cygwin*)
      AC_DEFINE(Cygwin, [], [Using cygwin])
      MACHID="cygwin"
      ;;
   *darwin*)
      AC_DEFINE(MACOSX, [], [Using MacOS X])
      MACHID="darwin"
      ;;
   *linux*)
      AC_DEFINE(Linux, [], [Using Linux])
      MACHID="linux"
      ;;
   *aix*)
      ac_mpicc_flags="-O3 -qstrict -qarch=pwr3 -qtune=pwr3"
      AC_DEFINE(AIX, [], [Using AIX])
      MACHID="aix"
      if /usr/sbin/lsattr -EHl proc0 | grep POWER >/dev/null 2>&1; then
          MACHID=rs6000
      fi
      ;;
   *cray*)
      ac_mpicc_flags="-O3 -h msglevel_3"
      AC_DEFINE(cray3e, [], [Using cray3e])
      MACHID="cray"
      ;;
   *irix6*)
      ac_mpicc_flags="-O3"
      AC_DEFINE(IRIX64, [], [Using IRIX64])
      MACHID="irix"
      ;;
   *irix*)
      ac_mpicc_flags="-O3"
      MACHID="irix"
      ;;
   *solaris*|*sunos*)
      if test "${CC}" != "gcc"; then
          ac_cc_opt="-O -Xa"
      fi
      AC_DEFINE(SunOS, [], [Using SunOS])
      AC_CHECK_LIB([nsl], [gethostbyname], [], [AC_MSG_ERROR([Nsl library not found])])
      AC_CHECK_LIB([socket], [socket], [], [AC_MSG_ERROR([Socket library not found])])
      LIBS="-lsocket -lnsl -lm"
      MACHID="solaris"
      ;;
esac

####################################
# Check for configuration switches
####################################

# debugging
AC_ARG_ENABLE(
    debug,
    [  --enable-debug           Produce an executable with debugging symbols. Default is disabled.],
    [if test "$enableval" = "yes"; then
        DEBUG_FLAGS="$ac_cc_debug";
     fi],
    [DEBUG_FLAGS=""])

# optimization
AC_ARG_ENABLE(
    opt,
    [  --enable-opt             Produce an executable with optimization. Default compiler dependent.],
    [if test "$enableval" = "yes"; then
       OPT_FLAGS="$ac_cc_opt";
     fi],
    [OPT_FLAGS="$ac_cc_opt"])

# get mpi info, and check if parallel version can be compiled
AC_ARG_WITH(
    mpicc,
    [  --with-mpicc=MPICC      Set MPICC to use. Default is found mpicc.],
    [MPICC="$withval"],
    [MPICC="$MPICC"])
    
AC_ARG_WITH(
    mpidir,
    [  --with-mpidir=MPIDIR    Set MPI installation path. Default is path of found mpicc.],
    [mpidir="$withval";
     MPICC=$mpidir/bin/mpicc;
     MPILIBDIR=$mpidir/lib;
     MPIINC="-DMPI -I$mpidir/include";
     AC_CHECK_LIB([mpi], [MPI_Comm_rank])
     parallel="yes";],
    [if test "$MPICC" = ""; then
        MPILIBDIR="";
        MPIINC="";
        parallel="no";
     else
        mpidir=`echo $MPICC | sed -e 's/\/bin\/.*//'`
        MPILIBDIR="$mpidir/lib"
        MPIINC="-DMPI -I$mpidir/include"
        AC_CHECK_LIB([mpi], [MPI_Comm_rank])
        parallel="yes";
     fi])

# enable only serial installation
AC_ARG_ENABLE(
    serial,
    [  --enable-serial          Produce only serial executables. Default parallel and serial.],
    serial=$enableval, serial=no)
if test "$serial" = "yes"; then
    parallel="no";
fi

if test "$parallel" = "yes"; then
   lammpi=`strings $MPICC | grep liblam`;
   if test -n "$lammpi" ; then
      MPI_CMD="\${prefix}/bin/lrun -n";
      MPITYPE="lam"
      LAMBIN="$mpidir/bin"
   else
      MPI_CMD="$mpidir/bin/mpirun -np";
      MPITYPE="mpich"
   fi
fi
AM_CONDITIONAL(PARALLEL, test ["$parallel" != no])
AM_CONDITIONAL(PARALLEL_IS_LAM, test -n ["$lammpi"]) 

# OpenMP support
AC_MSG_CHECKING(whether to enable OpenMP multithreading support)
AC_ARG_ENABLE(openmp, [  --enable-openmp          Enable on OpenMP multithreading support (default disabled)], enable_openmp=$enableval, enable_openmp=no)
if test "x$enable_openmp" = "xyes"; then
	AC_MSG_RESULT(yes)
	OPENMP_CFLAGS="$CFLAGS -fopenmp"
else
	AC_MSG_RESULT(no)
fi

AC_ARG_ENABLE(webservice,
		[  --enable-webservice=/path/to/build.xml		Deploy the MEME service on an existing Opal installation],
		[	if test "x$enableval" = "xyes"; then
				AC_MSG_ERROR(Path to build.xml not specified with --enable-webservice option!)
			else
				if test "" = "$CATALINA_HOME"; then
					AC_MSG_ERROR(Please ensure you have the environment variable \$CATALINA_HOME set to the base jakarta-tomcat directory.  Something like /opt/jakarta-tomcat-5.0.30)
				else
					OPAL_XML=$enableval
				fi
			fi
    ],OPAL_XML=no)
AM_CONDITIONAL(WEBSERVICE, test ["$OPAL_XML" != no])

# enable website installation
AC_ARG_ENABLE(
    web,
    [  --enable-web=opal_url          Build website applications. Default is disabled.  Specify the Opal URL for the service to communicate with.],
    [	if test "x$enableval" = "xyes"; then
				OPAL_URL="http://ws.nbcr.net/opal2/services"
				AC_MSG_WARN(Opal URL not specified with --enable-web option!  Defaulting to http://ws.nbcr.net/opal2/services)
			else
				OPAL_URL=$enableval
			fi
		],OPAL_URL=no)
AM_CONDITIONAL(WEBSITE, test ["$OPAL_URL" != no])

# enable building of included libxml2
AC_ARG_ENABLE(
    build_libxml2,
    [  --enable-build-libxml2   Build included version of libxml2. Default is to use localy insstalled copy.],
    build_libxml2=$enableval; \
    LIBXML2_LIBS="\${top_builddir}/src/libxml2/libxml2.la"; \
    LIBXML2_CFLAGS="-I\${top_srcdir}/src/libxml2/include" , )
AM_CONDITIONAL(BUILD_LIBXML2, test ["$build_libxml2" = yes])

# enable building of included libxslt
AC_ARG_ENABLE(
    build_libxslt,
    [  --enable-build-libxslt   Build included version of libxslt. Default is to use localy insstalled copy.],
    build_libxslt=$enableval; \
    LIBXSLT_LIBS="\${top_builddir}/src/libxslt/libxslt.la"; \
    LIBXSLT_CFLAGS="-I\${top_srcdir}/src/" , )
AM_CONDITIONAL(BUILD_LIBXSLT, test ["$build_libxslt" = yes])

# set website url
AC_ARG_WITH(
    url,
    [  --with-url=URL	  Set url for the website. Default http://`hostname`/meme.],
    [URL=$withval],
    [URL="http://$HOSTNAME/meme"])

# set BLOCKS server url
AC_ARG_WITH(
    blocks_url,
    [  --with-blocks-url=URL	  Set url for the BLOCKMAKER website. Default http://blocks.fhcrc.org/blocks-bin/process_blocks.pl],
    [BLOCKS_URL=$withval],
    [BLOCKS_URL="http://blocks.fhcrc.org/blocks-bin/process_blocks.pl"])
    
# set go server url
AC_ARG_WITH(
    goserver_url,
    [  --with-goserver-url=URL Set url for the Gene Ontology website. Default http://amigo.geneontology.org],
    [GOSERVER_URL=$withval],
    [GOSERVER_URL="http://amigo.geneontology.org/cgi-bin/amigo/"])
    
# set contact for website
AC_ARG_WITH(
    contact,
    [  --with-contact=EMAIL    Set email of website administrator. Default `whoami`@`hostname`.],
    [CONTACT=$withval],
    [CONTACT="$MEME_USER@$HOSTNAME"])

# set website location
AC_ARG_WITH(
    webdir,
    [  --with-webdir=DIR       Set path to website part. Default $prefix/web.],
    [WEBDIR=$withval],
    [WEBDIR="\${prefix}/web"])

# set database location
AC_ARG_WITH(
    db,
    [  --with-db=DBDIR         Set path to databases. Default $prefix/db.],
    [MEME_DB=$withval],
    [MEME_DB="\${prefix}/db"])

# set logs location
AC_ARG_WITH(
    logs,
    [  --with-logs=LOGDIR      Set path to logs. Default $prefix/LOGS.],
    [MEME_LOGS=$withval],
    [MEME_LOGS="\${prefix}/LOGS"])

# set path to Perl
AC_ARG_WITH(
    perl,
    [  --with-perl=NAME     Set path to Perl. Default is to take from user's PATH.],
    [PERL=$withval],
    [AC_PATH_PROG(PERL, perl)])
AC_DEFINE_UNQUOTED(PERL, "$PERL", [Path to Perl] )

# set path to Python
AC_ARG_WITH(
    python,
    [  --with-python=NAME     Set path to Python. Default is to take from user's PATH.],
    [PYTHON=$withval],
    [AC_PATH_PROG(PYTHON, python)])
AC_DEFINE_UNQUOTED(PYTHON, "$PYTHON", [Path to Python] )

# set path to convert utility
AC_ARG_WITH(
    convert,
    [  --with-convert=NAME     Set path to convert utility. Default is to take from user's PATH.],
    [CONVERT=$withval],
    [AC_PATH_PROG(CONVERT, convert)])
AC_DEFINE_UNQUOTED(CONVERT_PATH, "$CONVERT", [Path to convert utility] )


# set path to ghostscript utility (to replace convert)
AC_ARG_WITH(
    gs,
    [  --with-gs=NAME     Set path to ghostscript utility. Default is to take from user's PATH.],
    [GHOSTSCRIPT=$withval],
    [AC_PATH_PROG(GHOSTSCRIPT, gs)])
AC_DEFINE_UNQUOTED(GHOSTSCRIPT_PATH, "$GHOSTSCRIPT", [Path to ghostscript utility] )

# set path to dot utility
AC_ARG_WITH(
    dot,
    [  --with-dot=NAME	  Set path to dot utility. Default is to take from user's PATH.],
    [DOT=$withval],
    [AC_PATH_PROG(DOT, dot)])
AC_DEFINE_UNQUOTED(DOT_PATH, "$dot", [Path to dot utility] )

# set wall time limit 
AC_ARG_WITH(
    maxtime,
    [  --with-maxtime=TIME     Set wall time limit (sec) a web-submitted job may use. Default 7200.],
    [MAXTIME=$withval],
    [MAXTIME="7200"])

AM_CONDITIONAL(METAMEME_EXTERNAL, test ["$CONVERT" != "" -a "$DOT" != ""])
CFLAGS="$CFLAGS $DEBUG_FLAGS $OPT_FLAGS"
MPIFLAGS="-DPARALLEL $ac_mpicc_flags"
MEME_DIRECTORY="\${prefix}"

AC_SUBST(MEME_DIRECTORY)
AC_SUBST(MACHID)
AC_SUBST(MEME_PORT)
AC_SUBST(MAST_PORT)
AC_SUBST(DEBUG_FLAGS)
AC_SUBST(OPT_FLAGS)
AC_SUBST(LIBS)
AC_SUBST(MPICC)
AC_SUBST(MPIINC)
AC_SUBST(MPI_CMD)
AC_SUBST(MPILIBDIR)
AC_SUBST(MPIFLAGS)
AC_SUBST(MPITYPE)
AC_SUBST(VERSION)
AC_SUBST(LD)
AC_SUBST(LAMBIN)
AC_SUBST(URL)
AC_SUBST(BLOCKS_URL)
AC_SUBST(GOSERVER_URL)
AC_SUBST(CONTACT)
AC_SUBST(MAIL)
AC_SUBST(SENDMAIL)
AC_SUBST(HOSTNAME)
AC_SUBST(WEBSERVICE)
AC_SUBST(WEBDIR)
AC_SUBST(MEME_DB)
AC_SUBST(MEME_LOGS)
AC_SUBST(MAXTIME)
AC_SUBST(LIBXML2_CFLAGS)
AC_SUBST(LIBXML2_LIBS)
AC_SUBST(LIBXSLT_CFLAGS)
AC_SUBST(LIBXSLT_LIBS)
AC_SUBST(OPAL_URL)
AC_SUBST(ANT)
AC_SUBST(OPAL_XML)

AC_CONFIG_FILES([Makefile \
    doc/Makefile \
    doc/examples/Makefile \
    doc/examples/compute_prior_dist_example_output_files/Makefile \
    doc/examples/fimo_example_output_files/Makefile \
    doc/examples/glam2_example_output_files/Makefile \
    doc/examples/glam2scan_example_output_files/Makefile \
    doc/examples/mast_example_output_files/Makefile \
    doc/examples/mcast_example_output_files/Makefile \
    doc/examples/meme_example_output_files/Makefile \
    doc/examples/sample_opal_scripts/Makefile \
    doc/examples/tomtom_example_output_files/Makefile \
    doc/images/Makefile \
    etc/Makefile \
    src/filters/Makefile \
    src/filters/purge/Makefile \
    src/filters/dust/Makefile \
    scripts/Makefile \
    src/Makefile \
    src/libxml2/Makefile \
    src/libxml2/include/Makefile \
    src/libxml2/include/libxml/Makefile \
    src/libxslt/Makefile \
    src/parallel/Makefile \
    tests/Makefile \
    tests/common/Makefile \
    tests/mast/Makefile \
    tests/mcast/Makefile \
    tests/meme/Makefile \
    tests/clustalw2fasta/Makefile \
    tests/draw-mhmm/Makefile \
    tests/fasta-center/Makefile \
    tests/fimo/Makefile \
    tests/mhmm/Makefile \
    tests/mhmms/Makefile \
    tests/motiph/Makefile \
    tests/psp-gen/Makefile \
    tests/gomo/Makefile \
    tests/tomtom/Makefile \
    tests/scaffold/Makefile \
    tests/qvalue/Makefile \
    tests/web/Makefile \
    website/Makefile \
    website/cgi-bin/Makefile \
    website/html/Makefile \
    website/scripts/Makefile \
    website/images/Makefile])
AC_OUTPUT

# Raise warning if convert utility can't be found.
if test -z "$CONVERT" -o ! -e "$CONVERT"
then
  AC_MSG_WARN( \
    [The 'convert' utility can't be found.] \
  )
fi

# Raise warning if dot utility can't be found.
if test "$DOT" = ""
then
  AC_MSG_WARN( \
    [The 'dot' utility can't be found.] \
  )
fi

echo \
"
        ========================
        Configuration parameters
        ========================

  Install path:               ${prefix}
  Install UID:                ${MEME_USER}
  Version:                    ${VERSION}
  C compiler:                 ${CC}
  C compiler flags:           ${CFLAGS}
  Linker:                     ${LD}
  Special Libs:               ${LIBS}
  MPICC                       ${MPICC}
  MPIINC                      ${MPIINC}
  MPILIBDIR                   ${MPILIBDIR}
  MPIFLAGS                    ${MPIFLAGS}
  MPIRUN                      ${MPI_CMD}
  OpenMP Support:             ${enable_openmp}
  MEME URL:                   ${URL}
  MEME WEB DIRECTORY:         ${WEBDIR}
  MEME CONTACT:               ${CONTACT}
  MEME_LOGS:                  ${MEME_LOGS}
  BLOCKMAKER URL:             ${BLOCKS_URL}
  GO-SERVER URL:              ${GOSERVER_URL}
  PERL:                       ${PERL}
  PYTHON:                     ${PYTHON}
  CONVERT:                    ${CONVERT}
  DOT:                        ${DOT}
  BUILD LIBXML2:              ${build_libxml2}
  LIBXML2 compiler flags:     ${LIBXML2_CFLAGS}
  LIBXML2 libs:               ${LIBXML2_LIBS}
  BUILD LIBXSLT:              ${build_libxslt}
  LIBXSLT compiler flags:     ${LIBXSLT_CFLAGS}
  LIBXSLT libs:               ${LIBXSLT_LIBS}
  MEME_DB:                    ${MEME_DB}
  MEME_LOGS:                  ${MEME_LOGS}
  OPAL URL:                   ${OPAL_URL}
  OPAL XML:                   ${OPAL_XML}

  Run the following commands to compile, test and install meme:
        make   
        make test
        make install
"

