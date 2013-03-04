MESSAGE("\n")
MESSAGE("Make sure you have libgfortran in your LD_LIBRARY_PATH")
MESSAGE("Otherwise you get many of \"undefined reference\" errors \n")
enable_language (Fortran)

ADD_SUBDIRECTORY( ctq )
ADD_SUBDIRECTORY( gen )

