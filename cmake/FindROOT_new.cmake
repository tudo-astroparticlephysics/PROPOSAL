# copied from https://github.com/FairRootGroup/FairRoot

################################################################################
#    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
#                                                                              #
#              This software is distributed under the terms of the             #
#         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #
#                  copied verbatim in the file "LICENSE"                       #
################################################################################
# - Find ROOT instalation
# This module tries to find the ROOT installation on your system.
# It tries to find the root-config script which gives you all the needed 
# information.
# If the system variable ROOTSYS is set this is straight forward.
# If not the module uses the pathes given in ROOT_CONFIG_SEARCHPATH.
# If you need an other path you should add this path to this varaible.
# The root-config script is then used to detect basically everything else.
# This module defines a number of key variables and macros.
#
# Variables defined by this module:
#
#   ROOT_FOUND               System has ROOT, this means the root-config 
#                            executable was found.
#
#   ROOT_INCLUDE_DIR         ROOT include directories: not cached
#
#   ROOT_INCLUDES            Same as above,
#
#   ROOT_LIBRARIES           Link to these to use the ROOT libraries, not cached
#
#   ROOT_LIBRARY_DIR         The path to where the ROOT library files are.
#
#   ROOT_VERSION_STRING      The version string of the ROOT libraries which
#                            is reported by root-config
#
#   ROOT_VERSION_MAJOR       Major version number of ROOT
#   ROOT_VERSION_MINOR       Minor version number of ROOT
#   ROOT_VERSION_PATCH       Patch version number of ROOT
#
#   ROOT_VERSION_NUMBER      A unique version number which is calculated from 
#                            major, minor and patch version found
#
#   ROOT_CINT_EXECUTABLE     The rootcint executable.
#
#   RLIBMAP_EXECUTABLE       The rlibmap executable.

Message(STATUS "Looking for Root...")

Set(ROOT_CONFIG_SEARCHPATH
  ${ROOT_DIR}/bin
  ${SIMPATH}/bin
  ${SIMPATH}/tools/root/bin
  $ENV{ROOTSYS}/bin
)

Set(ROOT_FOUND FALSE)
Set(ROOT_DEFINITIONS "")
Set(ROOT_INSTALLED_VERSION_TOO_OLD FALSE)
Set(ROOT_CONFIG_EXECUTABLE ROOT_CONFIG_EXECUTABLE-NOTFOUND)

Find_Program(ROOT_CONFIG_EXECUTABLE
            NAMES root-config 
            PATHS ${ROOT_CONFIG_SEARCHPATH}
            NO_DEFAULT_PATH
            )

If(ROOT_CONFIG_EXECUTABLE)
   
  String(REGEX REPLACE "(^.*)/bin/root-config" "\\1" test ${ROOT_CONFIG_EXECUTABLE}) 
  Set(ENV{ROOTSYS} ${test})
  Set(ROOTSYS ${test})

  Execute_Process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --version 
                  OUTPUT_VARIABLE ROOT_VERSION_STRING
                 )
  Execute_Process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --prefix
                  OUTPUT_VARIABLE ROOT_INSTALL_DIR
                 )
  String(STRIP ${ROOT_VERSION_STRING} ROOT_VERSION_STRING)
  String(STRIP ${ROOT_INSTALL_DIR} ROOT_INSTALL_DIR)


  MESSAGE(STATUS "Looking for Root... - Found ${ROOT_INSTALL_DIR}/bin/root")
  MESSAGE(STATUS "Looking for Root... - Found version is ${ROOT_VERSION_STRING} ")   
   
  # extract major, minor, and patch versions from
  # the version string given by root-config
  String(REGEX REPLACE "^([0-9]+)\\.[0-9][0-9]+\\/[0-9][0-9]+.*" "\\1" ROOT_VERSION_MAJOR "${ROOT_VERSION_STRING}")
  String(REGEX REPLACE "^[0-9]+\\.([0-9][0-9])+\\/[0-9][0-9]+.*" "\\1" ROOT_VERSION_MINOR "${ROOT_VERSION_STRING}")
  String(REGEX REPLACE "^[0-9]+\\.[0-9][0-9]+\\/([0-9][0-9]+).*" "\\1" ROOT_VERSION_PATCH "${ROOT_VERSION_STRING}")

  # compute overall version numbers which can be compared at once
  Math(EXPR req_vers "${ROOT_VERSION_MAJOR}*10000 + ${ROOT_VERSION_MINOR}*100 + ${ROOT_VERSION_PATCH}")
  Math(EXPR found_vers "${ROOT_VERSION_MAJOR}*10000 + ${ROOT_VERSION_MINOR}*100 + ${ROOT_VERSION_PATCH}")
  Math(EXPR ROOT_FOUND_VERSION "${ROOT_VERSION_MAJOR}*10000 + ${ROOT_VERSION_MINOR}*100 + ${ROOT_VERSION_PATCH}")

  Set(ROOT_Version ${found_vers})
  Set(ROOT_VERSION_NUMBER ${found_vers})

  If(found_vers LESS req_vers)
    Set(ROOT_FOUND FALSE)
    Set(ROOT_INSTALLED_VERSION_TOO_OLD TRUE)
  Else(found_vers LESS req_vers)
    Set(ROOT_FOUND TRUE)
  EndIf(found_vers LESS req_vers)

Else(ROOT_CONFIG_EXECUTABLE)
  Message(STATUS "Looking for Root... - Not found")
  Message(FATAL_ERROR "ROOT not installed in the searchpath and ROOTSYS is not set. Please set ROOTSYS or add the path to your ROOT installation in the Macro FindROOT.cmake in the subdirectory cmake/modules.")
Endif(ROOT_CONFIG_EXECUTABLE)


If(ROOT_FOUND)

  # ask root-config for the library dir
  # Set ROOT_LIBRARY_DIR
  Execute_Process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --libdir
                  OUTPUT_VARIABLE ROOT_LIBRARY_DIR
                 )
  String(STRIP ${ROOT_LIBRARY_DIR} ROOT_LIBRARY_DIR)

  # ask root-config for the binary dir
  Execute_Process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --bindir
                  OUTPUT_VARIABLE ROOT_BINARY_DIR
                 )
  String(STRIP ${ROOT_BINARY_DIR} ROOT_BINARY_DIR)

  # ask root-config for the include dir
  Execute_Process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --incdir
                  OUTPUT_VARIABLE ROOT_INCLUDE_DIR
                 )
  String(STRIP ${ROOT_INCLUDE_DIR} ROOT_INCLUDE_DIR)

  # ask root-config for the library varaibles
  Execute_Process(COMMAND ${ROOT_CONFIG_EXECUTABLE} --glibs
                  OUTPUT_VARIABLE ROOT_LIBRARIES
                 )
  String(STRIP ${ROOT_LIBRARIES} ROOT_LIBRARIES)

  # Make variables changeble to the advanced user
  Mark_As_Advanced(ROOT_LIBRARY_DIR ROOT_INCLUDE_DIR ROOT_DEFINITIONS)

  # Set ROOT_INCLUDES
  Set(ROOT_INCLUDES ${ROOT_INCLUDE_DIR})

  Set(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${ROOT_LIBRARY_DIR})

  #######################################
  #
  #       Check the executables of ROOT 
  #          ( rootcint ) 
  #
  #######################################

  Find_Program(ROOT_CINT_EXECUTABLE
    NAMES rootcint
    PATHS ${ROOT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  Find_Program(RLIBMAP_EXECUTABLE
    NAMES rlibmap
    PATHS ${ROOT_BINARY_DIR}
    NO_DEFAULT_PATH
    )

  Function(Format _output input prefix suffix)

  # DevNotes - input should be put in quotes or the complete list does not get passed to the function
    set(format)
    foreach(arg ${input})
      set(item ${arg})
      if(prefix)
        string(REGEX MATCH "^${prefix}" pre ${arg})
      endif(prefix)
      if(suffix)
        string(REGEX MATCH "${suffix}$" suf ${arg})
      endif(suffix)
      if(NOT pre)
        set(item "${prefix}${item}")
      endif(NOT pre)
      if(NOT suf)
        set(item "${item}${suffix}")
      endif(NOT suf)
      list(APPEND format ${item})
    endforeach(arg)
    set(${_output} ${format} PARENT_SCOPE)

  endfunction(Format)


    ###########################################
    #
    #       Macros for building ROOT dictionary
    #
    ###########################################
  Macro(ROOT_GENERATE_DICTIONARY)

    # Macro to switch between the old implementation with parameters
    # and the new implementation without parameters.
    # For the new implementation some CMake variables has to be defined
    # before calling the macro.

    If(${ARGC} EQUAL 0)
  #    Message("New Version")
      ROOT_GENERATE_DICTIONARY_NEW()
    Else(${ARGC} EQUAL 0)
      If(${ARGC} EQUAL 4)
  #      Message("Old Version")
        ROOT_GENERATE_DICTIONARY_OLD("${ARGV0}" "${ARGV1}" "${ARGV2}" "${ARGV3}")
      Else(${ARGC} EQUAL 4)
        Message(FATAL_ERROR "Has to be implemented")
      EndIf(${ARGC} EQUAL 4)
    EndIf(${ARGC} EQUAL 0)

  EndMacro(ROOT_GENERATE_DICTIONARY)

  Macro(ROOT_GENERATE_DICTIONARY_NEW)

    # All Arguments needed for this new version of the macro are defined
    # in the parent scope, namely in the CMakeLists.txt of the submodule
    set(Int_LINKDEF ${LINKDEF})
    set(Int_DICTIONARY ${DICTIONARY})
    set(Int_LIB ${LIBRARY_NAME})

    set(Int_INC ${INCLUDE_DIRECTORIES} ${SYSTEM_INCLUDE_DIRECTORIES})
    set(Int_HDRS ${HDRS})
    set(Int_DEF ${DEFINITIONS})

    # Convert the values of the variable to a semi-colon separated list
    separate_arguments(Int_INC)
    separate_arguments(Int_HDRS)
    separate_arguments(Int_DEF)

    # Format neccesary arguments
    # Add -I and -D to include directories and definitions
    Format(Int_INC "${Int_INC}" "-I" "")
    Format(Int_DEF "${Int_DEF}" "-D" "")

    #---call rootcint / cling --------------------------------
    set(OUTPUT_FILES ${Int_DICTIONARY})

    If (CMAKE_SYSTEM_NAME MATCHES Linux)
      Set(MY_LD_LIBRARY_PATH ${ROOT_LIBRARY_DIR}:${_intel_lib_dirs}:$ENV{LD_LIBRARY_PATH})
    ElseIf(CMAKE_SYSTEM_NAME MATCHES Darwin)
      Set(MY_LD_LIBRARY_PATH ${ROOT_LIBRARY_DIR}:$ENV{DYLD_LIBRARY_PATH})
    EndIf()

    get_filename_component(script_name ${Int_DICTIONARY} NAME_WE)
    String(REPLACE ";" " " Int_DEF_STR "${Int_DEF}")
    String(REPLACE ";" " " Int_INC_STR "${Int_INC}")
    String(REPLACE ";" " " Int_HDRS_STR "${Int_HDRS}")

    Set(EXTRA_DICT_PARAMETERS "")
    If (ROOT_FOUND_VERSION GREATER 59999)

      Set(Int_ROOTMAPFILE ${LIBRARY_OUTPUT_PATH}/lib${Int_LIB}.rootmap)
      Set(Int_PCMFILE G__${Int_LIB}Dict_rdict.pcm)
      Set(OUTPUT_FILES ${OUTPUT_FILES} ${Int_PCMFILE} ${Int_ROOTMAPFILE})
      Set(EXTRA_DICT_PARAMETERS ${EXTRA_DICT_PARAMETERS}
          -inlineInputHeader -rmf ${Int_ROOTMAPFILE} 
          -rml ${Int_LIB}${CMAKE_SHARED_LIBRARY_SUFFIX})
      Set_Source_Files_Properties(${OUTPUT_FILES} PROPERTIES GENERATED TRUE)
      String(REPLACE ";" " " EXTRA_DICT_PARAMETERS_STR "${EXTRA_DICT_PARAMETERS}")

    EndIf()

    # We need some environment variables which are present when running cmake at the
    # time we run make. To pass the variables a script is created containing the
    # correct values for the needed variables

    IF(FAIRROOTPATH)
      Configure_File(${FAIRROOTPATH}/share/fairbase/cmake/scripts/generate_dictionary_root.sh.in
                     ${CMAKE_CURRENT_BINARY_DIR}/generate_dictionary_${script_name}.sh
                    )
      EXEC_PROGRAM(/bin/chmod ARGS "u+x ${CMAKE_CURRENT_BINARY_DIR}/generate_dictionary_${script_name}.sh")
    ELSE(FAIRROOTPATH)
      Configure_File(${PROJECT_SOURCE_DIR}/cmake/scripts/generate_dictionary_root.sh.in
                     ${CMAKE_CURRENT_BINARY_DIR}/generate_dictionary_${script_name}.sh
                    )
    ENDIF(FAIRROOTPATH)


    If (ROOT_FOUND_VERSION GREATER 59999)
      Add_Custom_Command(OUTPUT  ${OUTPUT_FILES}
                         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/generate_dictionary_${script_name}.sh
                         COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/${Int_PCMFILE} ${LIBRARY_OUTPUT_PATH}/${Int_PCMFILE}
                         DEPENDS ${Int_HDRS} ${Int_LINKDEF}
                        )
      Install(FILES ${LIBRARY_OUTPUT_PATH}/${Int_PCMFILE} ${Int_ROOTMAPFILE} DESTINATION lib)
    Else()
      Add_Custom_Command(OUTPUT  ${OUTPUT_FILES}
                         COMMAND ${CMAKE_CURRENT_BINARY_DIR}/generate_dictionary_${script_name}.sh
                         DEPENDS ${Int_HDRS} ${Int_LINKDEF}
                        )
    EndIf()

  endmacro(ROOT_GENERATE_DICTIONARY_NEW)


  MACRO (ROOT_GENERATE_DICTIONARY_OLD INFILES LINKDEF_FILE OUTFILE INCLUDE_DIRS_IN)

    set(INCLUDE_DIRS)

    foreach (_current_FILE ${INCLUDE_DIRS_IN})
      set(INCLUDE_DIRS ${INCLUDE_DIRS} -I${_current_FILE})
    endforeach (_current_FILE ${INCLUDE_DIRS_IN})

  #  Message("Definitions: ${DEFINITIONS}")
  #  MESSAGE("INFILES: ${INFILES}")
  #  MESSAGE("OutFILE: ${OUTFILE}")
  #  MESSAGE("LINKDEF_FILE: ${LINKDEF_FILE}")
  #  MESSAGE("INCLUDE_DIRS: ${INCLUDE_DIRS}")

    STRING(REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" bla "${OUTFILE}")
  #  MESSAGE("BLA: ${bla}")
    SET (OUTFILES ${OUTFILE} ${bla})


    if (CMAKE_SYSTEM_NAME MATCHES Linux)
      ADD_CUSTOM_COMMAND(OUTPUT ${OUTFILES}
         COMMAND LD_LIBRARY_PATH=${ROOT_LIBRARY_DIR}:${_intel_lib_dirs} ROOTSYS=${ROOTSYS} ${ROOT_CINT_EXECUTABLE}
         ARGS -f ${OUTFILE} -c -DHAVE_CONFIG_H ${INCLUDE_DIRS} ${INFILES} ${LINKDEF_FILE} DEPENDS ${INFILES} ${LINKDEF_FILE})
    else (CMAKE_SYSTEM_NAME MATCHES Linux)
      if (CMAKE_SYSTEM_NAME MATCHES Darwin)
        ADD_CUSTOM_COMMAND(OUTPUT ${OUTFILES}
         COMMAND DYLD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS} ${ROOT_CINT_EXECUTABLE}
         ARGS -f ${OUTFILE} -c -DHAVE_CONFIG_H ${INCLUDE_DIRS} ${INFILES} ${LINKDEF_FILE} DEPENDS ${INFILES} ${LINKDEF_FILE})
      endif (CMAKE_SYSTEM_NAME MATCHES Darwin)
    endif (CMAKE_SYSTEM_NAME MATCHES Linux)

  ENDMACRO (ROOT_GENERATE_DICTIONARY_OLD)

  MACRO (GENERATE_ROOT_TEST_SCRIPT SCRIPT_FULL_NAME)

    get_filename_component(path_name ${SCRIPT_FULL_NAME} PATH)
    get_filename_component(file_extension ${SCRIPT_FULL_NAME} EXT)
    get_filename_component(file_name ${SCRIPT_FULL_NAME} NAME_WE)
    set(shell_script_name "${file_name}.sh")

    #MESSAGE("PATH: ${path_name}")
    #MESSAGE("Ext: ${file_extension}")
    #MESSAGE("Name: ${file_name}")
    #MESSAGE("Shell Name: ${shell_script_name}")

    string(REPLACE ${PROJECT_SOURCE_DIR}
           ${PROJECT_BINARY_DIR} new_path ${path_name}
          )

    #MESSAGE("New PATH: ${new_path}")

    file(MAKE_DIRECTORY ${new_path}/data)

    CONVERT_LIST_TO_STRING(${LD_LIBRARY_PATH})
    set(MY_LD_LIBRARY_PATH ${output})

    CONVERT_LIST_TO_STRING(${ROOT_INCLUDE_PATH})
    set(MY_ROOT_INCLUDE_PATH ${output})

    set(my_script_name ${SCRIPT_FULL_NAME})

    IF(FAIRROOTPATH)
      configure_file(${FAIRROOTPATH}/share/fairbase/cmake/scripts/root_macro.sh.in
                     ${new_path}/${shell_script_name}
                    )
    ELSE(FAIRROOTPATH)
      configure_file(${PROJECT_SOURCE_DIR}/cmake/scripts/root_macro.sh.in
                     ${new_path}/${shell_script_name}
                    )
    ENDIF(FAIRROOTPATH)

    EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${new_path}/${shell_script_name}")

  ENDMACRO (GENERATE_ROOT_TEST_SCRIPT)

Else(ROOT_FOUND)

  If(ROOT_FIND_REQUIRED)
    Message(STATUS "Looking for ROOT... - Found version to old.")
    Message(STATUS "Looking for ROOT... - Minimum required version is ${ROOT_FIND_VERSION}")
    Message(FATAL_ERROR "Stop here because of a wrong Root version.")
  EndIf(ROOT_FIND_REQUIRED)

Endif(ROOT_FOUND)
