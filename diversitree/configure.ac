AC_INIT([diversitree], 1.0, [rich.fitzjohn@gmail.com])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`

PKG_CPPFLAGS="-I."

dnl FFTW library location
AC_ARG_WITH([fftw],
	AS_HELP_STRING([--with-fftw=PATH],
		[the location of fftw libraries and headers]),
	[fftw_path=$withval])
if test [ -n $fftw_path ]; then
  if test [ -d "${fftw_path}/lib" ]; then
     PKG_LIBS="-L$fftw_path/lib ${PKG_LIBS}"
  fi
  if test [ -d "${fftw_path}/include" ]; then
     PKG_CPPFLAGS="-I$fftw_path/include ${PKG_CPPFLAGS}"
  fi
fi

# GSL configuration based on package 'gsl'
## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
   # Use gsl-config for header and linker arguments
   GSL_CFLAGS=`${GSL_CONFIG} --cflags`
   GSL_LIBS=`${GSL_CONFIG} --libs`
else
   AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi

AC_MSG_CHECKING([if GSL version >= 1.12])
AC_RUN_IFELSE([AC_LANG_SOURCE([[
#include <string.h>
#include <gsl/gsl_version.h>
int main() {
#ifdef GSL_VERSION
   int major, minor;
   char *gslv = GSL_VERSION;
   if ((sscanf(gslv, "%d.%d", &major, &minor)) != 2) {
     exit (1);
   }
   exit (minor < 12);
#else
  exit(1);
#endif
}
]])],
[gsl_version_ok=yes],
[gsl_version_ok=no],
[gsl_version_ok=yes])
if test "${gsl_version_ok}" = no; then
  AC_MSG_ERROR([Need GSL version >= 1.12])
else
  AC_MSG_RESULT([yes])
fi

CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"
dnl AC_CHECK_HEADERS([fftw3.h], AC_DEFINE(HAVE_FFTW3_H, 1, Define to 1 if you have FFTW))
AC_CHECK_HEADERS([fftw3.h])
if test "${ac_cv_header_fftw3_h}" = no; then
   AC_MSG_WARN(No fftw found - QuaSSE/fftC will not be available)
else
   echo "Using fftw for QuaSSE FFT calculations"
   PKG_LIBS="${PKG_LIBS} -lfftw3"
fi

AC_ARG_ENABLE([gcc-strict],
	AS_HELP_STRING([--enable-gcc-strict],
		[use strict gcc error checking]),
	[gcc_strict=$enableval])
if [ test x"$gcc_strict" = xyes ]; then
   echo "Using strict C checking"
   PKG_CFLAGS="-Wall -pedantic ${PKG_CFLAGS}"
else
   echo "Using normal C checking"
fi

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LIBS)
AC_SUBST(GSL_CFLAGS) 
AC_SUBST(GSL_LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_CONFIG_HEADERS([src/config.h])
AC_OUTPUT
