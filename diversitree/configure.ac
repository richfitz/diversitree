AC_INIT([diversitree], 0.6, [fitzjohn@zoology.ubc.ca])

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

AC_CONFIG_HEADERS([src/config.h])

CPPFLAGS="${CPPFLAGS} ${PKG_CPPFLAGS}"
dnl AC_CHECK_HEADERS([fftw3.h], AC_DEFINE(HAVE_FFTW3_H, 1, Define to 1 if you have FFTW))
AC_CHECK_HEADERS([fftw3.h])
if test "${ac_cv_header_fftw3_h}" = no; then
dnl AC_MSG_ERROR("FFTW header fftw3.h not found")
   AC_MSG_WARN(No fftw found - QuaSSE/fftC will not be available)
else
   echo "Using fftw for QuaSSE FFT calculations"
   PKG_LIBS="${PKG_LIBS} -lfftw3"
fi

dnl Faster gcc flags
AC_ARG_ENABLE([gcc-speedy],
	AS_HELP_STRING([--enable-gcc-speedy],
		[use agressive optimisation with gcc]),
	[gcc_speedy=$enableval])
if [ test x"$gcc_speedy" = xyes ]; then
	PKG_CFLAGS="-O3 -ffast-math -funroll-all-loops -fomit-frame-pointer"
	echo "Using fast options"
else
	PKG_CFLAGS=""
	echo "Using default options"
fi

AC_ARG_ENABLE([gcc-strict],
	AS_HELP_STRING([--enable-gcc-strict],
		[use strict gcc error checking]),
	[gcc_strict=$enableval])
if [ test x"$gcc_strict" = xyes ]; then
   echo "Using strict C checking"
   PKG_CFLAGS="-Wall -ansi -pedantic ${PKG_CFLAGS}"
else
   echo "Using normal C checking"
fi

AC_SUBST(PKG_CPPFLAGS)
AC_SUBST(PKG_CFLAGS)
AC_SUBST(PKG_LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
