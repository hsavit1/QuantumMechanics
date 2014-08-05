################################################################
# 'make'        build executables files - one for each source
# 'make clean'  removes all object, dependency, and executable files
################################################################

###################### # C compiler # ##########################
# define the C compiler to use (Supported: Intel)
CC = Intel

################### # Compile-time flags # #####################
# define any compile-time flags

# Use this for common flags (gets - in linux and /Q in windows)
CFLAGS = std=c++11

# Use these for platform specific flags (write verbose e.g. including /Q if needed!)
WIN_CFLAGS = 

# Use these for platform specific flags (write verbose e.g. including - if needed!)
LINUX_CFLAGS =

################# # Compiled Library Paths # ###################
# define library paths in addition to standard

# Use this for common paths (gets -L in linux and /L in windows)
LPATHS = 

# Use these for platform specific paths (write verbose e.g. including /L if needed!)
WIN_LPATHS = 

# Use these for platform specific paths (write verbose e.g. including -L if needed!)
LINUX_LPATHS =

################### # Compiled Libraries # #####################
# define any libraries to link into executable

# Use this for common libs (gets - in linux and /Q in windows)
LIBS = mkl tbb

# Use these for platform specific libs (write verbose e.g. including /Q if needed!)
WIN_LIBS = 

# Use these for platform specific libs (write verbose e.g. including - if needed!)
LINUX_LIBS =

################# # Header Files/Libraries # ###################
# define any directories containing header files other than standard


# Use this for common libs (gets - in linux and /Q in windows)
INCLUDES = M:/Code/Libraries/ M:/Code/Includes

# Use these for platform specific libs (write verbose e.g. including /Q if needed!)
WIN_INCLUDES = 

# Use these for platform specific libs (write verbose e.g. including - if needed!)
LINUX_INCLUDES =

################################################################
###################### # Source Files # ########################
################################################################

# define the main source files
SRCS = GreensFormalismUnittesting.cpp \
	EigensystemUnittesting.cpp
	
################################################################
####################### # Main Files # #########################
################################################################

# define the executable file 
MAIN = $(SRCS:%.cpp=%)
	
################################################################
################################################################
#
#	The following part of the makefile is generic; it can be used to 
#	build any executable just by changing the definitions above.
#	Do not change anything below with out testing first.
#
################################################################
################################################################

# find the C compiler to use
ifeq ($(OS),Windows_NT)
	ifeq ($(CC),Intel)
		-include makefile_intel_win64
	endif
else
	ifeq ($(CC),Intel)
		-include makefile_intel_linux64
	endif
endif