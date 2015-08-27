CFLAGS = -ansi -pedantic -w -g -O3 -fPIC -I$(NETCDF_INCDIR) #-DDEBUG $(EXTRA_FLAGS)
#CXX = OMPI_CXX=/opt/intel/bin/icpc mpic++ -DOMPI_SKIP_MPICXX
CXX = mpicxx -DOMPI_SKIP_MPICXX
CC = mpicc

vpath  %c src
vpath  %cpp src

SRC_C = gpc.c 
SRC_CPP = circle.cpp cputime.cpp grid.cpp inside.cpp intersect.cpp intersection_ym.cpp\
	mapper.cpp meshutil.cpp mpi_routing.cpp mpi_cascade.cpp \
	node.cpp parallel_tree.cpp polyg.cpp \
	timer.cpp tree.cpp triple.cpp \
	libmapper.cpp clipper.cpp

OBJ_C =   $(patsubst %.c,obj/%.o, $(SRC_C))
OBJ_CPP = $(patsubst %.cpp,obj/%.o,$(SRC_CPP))

lib/libmapper.so: $(OBJ_CPP) $(OBJ_C)
	mkdir -p lib
	echo $(OBJ_C)
	$(CXX) -shared -o $@ $^ -lc

test: $(OBJ_CPP) $(OBJ_C) obj/test-main.o
	$(CXX) $(CFLAGS) -o $@  $^ -L$(NETCDF_LIBDIR) -lnetcdf

try: test
	cd run && mpirun -n 32 ./test

obj/%.o: %.cpp
	@mkdir -p obj
	$(CXX) -c $< -o $@ $(CFLAGS)

obj/%.o: %.c
	@mkdir -p obj
	$(CXX) -c $< -o $@ $(CFLAGS)

clean:
	rm -rf obj/* lib/*
