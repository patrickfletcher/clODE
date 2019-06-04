OBJS1 = testTransient.o CLODE.o OpenCLResource.o
OBJS2 = testTrajectory.o CLODE.o CLODEtrajectory.o OpenCLResource.o
OBJS3 = testFeatures.o CLODE.o CLODEfeatures.o OpenCLResource.o
CXX = g++
DEBUG = 
CPPFLAGS = -Wall -c -std=c++0x $(DEBUG)
LFLAGS = -Wall $(DEBUG)
LDLIBS=-lOpenCL

CLODEDIR := ${CURDIR}

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	CPPFLAGS += -I/usr/local/cuda-7.5/targets/x86_64-linux/include
endif
ifeq ($(UNAME_S),Darwin)
	CPPFLAGS += -framework OpenCL
endif

all: testTrans testTraj testFeat

testTrans : $(OBJS1)
	$(CXX) $(LFLAGS) -o testTrans $(OBJS1) $(LDLIBS)
	
testTraj : $(OBJS2)
	$(CXX) $(LFLAGS) -o testTraj $(OBJS2) $(LDLIBS)
	
testFeat : $(OBJS3)
	$(CXX) $(LFLAGS) -o testFeat $(OBJS3) $(LDLIBS)

testTransient.o: testTransient.cpp OpenCLResource.hpp CLODE.hpp
	$(CXX) $(CPPFLAGS) testTransient.cpp 
	
testTrajectory.o: testTrajectory.cpp OpenCLResource.hpp CLODE.hpp CLODEtrajectory.hpp
	$(CXX) $(CPPFLAGS) testTrajectory.cpp 
	
testFeatures.o: testFeatures.cpp OpenCLResource.hpp CLODE.hpp CLODEfeatures.hpp
	$(CXX) $(CPPFLAGS) testFeatures.cpp 

CLODE.o : CLODE.cpp CLODE.hpp
	$(CXX) $(CPPFLAGS) CLODE.cpp -DCLODE_ROOT=\"$(CLODEDIR)/\"
	
CLODEtrajectory.o : CLODEtrajectory.cpp CLODEtrajectory.hpp
	$(CXX) $(CPPFLAGS) CLODEtrajectory.cpp -DCLODE_ROOT=\"$(CLODEDIR)/\"
	
CLODEfeatures.o : CLODEfeatures.cpp CLODEfeatures.hpp
	$(CXX) $(CPPFLAGS) CLODEfeatures.cpp -DCLODE_ROOT=\"$(CLODEDIR)/\"

OpenCLResource.o : OpenCLResource.cpp  OpenCLResource.hpp
	$(CXX) $(CPPFLAGS) OpenCLResource.cpp
	
.PHONY: clean
clean:
	\rm *.o testTrans testTraj testFeat
