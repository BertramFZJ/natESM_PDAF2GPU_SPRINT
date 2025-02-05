FC = nvfortran

SMARCH = -gpu=sm_86
LTOFLAG =
CUDAVER = -gpu=cuda12.3

COMPILER_OPTIONS = -mp -acc=gpu -gpu=ptxinfo $(CUDAVER) $(SMARCH) $(LTOFLAG) -Minfo=accel -cpp -O2 -Mnoautoinline -Mstandard -Minform=inform

LINKER_OPTIONS = -mp -acc=gpu $(SMARCH) $(CUDAVER) $(LTOFLAG)

TARGET = fortranMemoryHeap

SOURCES = mo_mem_workspaces.F90 main.F90

OBJECTS = $(patsubst %.F90, %.o, $(SOURCES))

%.o : %.F90
	@echo "compiling $^ ..."
	@$(FC) $(COMPILER_OPTIONS) -c $^ -o $@

all : $(OBJECTS)
	@rm -f $(TARGET).exe
	@echo "building $(TARGET) ..."
	@$(FC) $(LINKER_OPTIONS) $(OBJECTS) -o $(TARGET).exe

clean :
	@rm -f $(TARGET).exe $(OBJECTS) *.mod
