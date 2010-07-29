#!/bin/sh
#
# $Id: pre-install_trawler.sh,v 1.1 2009/04/28 07:34:00 haudry Exp $
#
# Trawler pre-installation script
# @author Yannick Haudry (haudry@embl.de)   #
#
# Run with --help to see more information
#

#uncomment to enable debug
#set -x

#############
# CONSTANTS #
#############

PERL_REQUIRED_VERSION="5.006.000"
JAVA_REQUIRED_VERSION="1.5"

TRAWLER_URL="http://ani.embl.de/laurence/blog/trawler_download/trawler_standalone-1.0.tar.gz"
TRAWLER_PACKAGE="trawler_standalone-1.0.tar.gz"

#################
# DOCUMENTATION #
#################

# command line usage
usage () {
    echo "`basename $0` version $VERSION"
    cat <<USAGE
$0 [-h|--help]|[-u|--usage]|[-v|--version]
USAGE
}

help () {
    echo "install_trawler.sh version $VERSION"
    usage
    cat <<HELP

  DESCRIPTION:

      This script checks the availability of required dependencies
      to run Trawler on your local machine.

  OPTIONS:
      -h|--help    - print this help and exit
      -u|--usage   - print usage and exit
      -v|--version - print version and exit

  AUTHOR:
      yannick.haudry@embl.de

HELP
}

#############
# FUNCTIONS #
#############

ACTION=""
fail () {
    echo "$ACTION failed: aborting!" 1>&2
    exit 3;
}

success () {
    echo "OK" 1>&2
}

###############################
# CHECK REQUIRED DEPENDENCIES #
###############################

# 1. Check for BASH environment
check_bash () {
    if [ -z "$BASH" ]; then
        echo "Please run this script using bash interpreter," 1>&2
        echo "using the following command: > bash install_trawler.sh" 1>&2
        exit 1;
    fi
}

# 2. Check for Perl installation / version
check_perl () {
    if ! perl -M$PERL_REQUIRED_VERSION -e1; then
        echo
        echo "FAILED" 1>&2
        echo "Aborting: Trawler requires Perl >= v5.6.0" 1>&2
        echo "You can use your package manager or go to http://www.perl.org/ for more information." 1>&2
        # exit 2;
    else
        success
    fi
}

# 3. Check for Java installation / version
check_java () {
    # Transform the required version string into a number that can be used in comparisons
    JAVA_REQUIRED_VERSION=`echo $JAVA_REQUIRED_VERSION | sed -e 's;\.;0;g'`

    # Check JAVA_HOME directory to see if Java version is adequate
    java -version 2> tmp.ver
    VERSION=`cat tmp.ver | grep "java version" | awk '{ print substr($3, 2, length($3)-2); }'`
    rm tmp.ver
    VERSION=`echo $VERSION | awk '{ print substr($1, 1, 3); }' | sed -e 's;\.;0;g'`
    if [ $VERSION ]; then
        if [ $VERSION -ge $JAVA_REQUIRED_VERSION ]; then
            JAVA_HOME=`echo java | awk '{ print substr($1, 1, length($1)-9); }'`
            success
        else
            echo "FAILED" 1>&2
            echo "Aborting: Trawler requires Java >= v1.5" 1>&2
            echo "Please refer to the Java website: http://java.sun.com/javase/downloads/index.jsp" 1>&2
            # exit 2;
        fi
    fi
}

# 4. Check for Ghostscript
check_gs () {
    gs -version >/dev/null 2>&1
    if [ $? != 0 ]; then
        echo "FAILED" 1>&2
        echo "Aborting: Trawler requires Ghostscript" 1>&2
        echo "You can use you package manager or go to http://pages.cs.wisc.edu/~ghost/ for more information" 1>&2
        # exit 2;
    else
        success
    fi
}

# 5. Algorithm::Cluster
check_cluster () {
    perl -e 'use Algorithm::Cluster'
    if [ $? != 0 ]; then
        echo "FAILED" 1>&2
        echo "Aborting: Trawler requires Algorithm::Cluster" 1>&2
        echo "read the trawler install instructions or the install notes from http://search.cpan.org/~mdehoon/"
        # exit 2;
    else
        success
    fi
}

# Run all dependencies checking
check_all () {
    check_bash
    echo "--- Checking Perl installation ---"
    check_perl
    echo
    echo "--- Checking Java installation ---"
    check_java
    echo
    echo "--- Checking Ghostscript installation ---"
    check_gs
    echo
    echo "--- Checking Algorithm::Cluster installation ---"
    check_cluster
    echo
}


######################
# ARGUMENTS HANDLING #
######################

VERSION=0.1
PRINT_USAGE=0
PRINT_HELP=0
PRINT_VERSION=0

while [ $# -gt 0 ]; do
    case "$1" in
        -u|--usage) PRINT_USAGE=1; shift ;;
        -h|--help) PRINT_HELP=1; shift ;;
        -v|--version) PRINT_VERSION=1; shift ;;
        --) shift ; break ;;
        -*) echo "Unknown command-line option: $1" ; exit 1 ;;
        *)  echo "breaking"; break;; # terminate while loop
    esac
    shift
done

if [ "$PRINT_USAGE" = 1 ]; then usage; exit; fi
if [ "$PRINT_VERSION" = 1 ]; then echo Version: $VERSION; exit; fi
if [ "$PRINT_HELP" = 1 ]; then help; exit; fi

check_all

exit 0;
