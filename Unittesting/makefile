#
# 'make'        build executable file 'mycc'
# 'make clean'  removes all .o and executable files
#

# define the C compiler to use
CC = icpc

# define any compile-time flags
CFLAGS = -std=c++11 -DMKL_Complex8="std::complex<float>" -DMKL_Complex16="std::complex<double>" -DEIGEN_USE_MKL_ALL -O3 -ip

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -mkl -tbb

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = 

# define any directories containing header files other than /usr/include
INCLUDES = 

# define the C source files
SRCS = GreensFormalismUnittesting.cpp

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

# define the executable file 
MAIN = $(SRCS:%.cpp=%Executable)

DEPDIR = .deps
DEPEXT = dep
DEPS = $(SRCS:%.cpp=$(DEPDIR)/%.$(DEPEXT))
OBJDIR = .obj
OBJS = $(SRCS:%.cpp=$(OBJDIR)/%.o)

.PHONY: clean

all:    $(MAIN)
	@echo  Executable \'$(MAIN)\' has been compiled.

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LFLAGS) -o $@ $<


$(OBJS): $(SRCS) $(OBJDIR) $(DEPDIR) makefile
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LFLAGS) -MMD -MP -MF$(<:%.cpp=$(DEPDIR)/%.$(DEPEXT)) -c -o $@ $<

$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(DEPDIR):
	@mkdir -p $(DEPDIR)

clean:
	$(RM) *~ *.bak $(OBJS) $(DEPS) $(MAIN)
	@rmdir $(OBJDIR) $(DEPDIR)


-include $(SRCS:%.cpp=$(DEPDIR)/%.$(DEPEXT))