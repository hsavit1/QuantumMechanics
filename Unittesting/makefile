#
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C source files
SRCS = GreensFormalismUnittesting.cpp \
	EigensystemUnittesting.cpp \
	CircularBufferUnitTesting.cpp

# define the C compiler to use
ifeq ($(OS),Windows_NT)
	-include makefile_win64
else
	-include makefile_linux64
endif
