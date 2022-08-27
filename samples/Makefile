OBJS1 = testTransient.o CLODE.o OpenCLResource.o
OBJS2 = testTrajectory.o CLODE.o CLODEtrajectory.o OpenCLResource.o
OBJS3 = testFeatures.o CLODE.o CLODEfeatures.o OpenCLResource.o
CXX = g++
DEBUG = 
CPPFLAGS = -Wall -c -std=c++0x $(DEBUG) -I../src
LFLAGS = -Wall $(DEBUG)
LDLIBS=

CLODEDIR := ${CURDIR}

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	CPPFLAGS += -I/usr/local/cuda-7.5/targets/x86_64-linux/include
	LDFLAGS += -DCLODE_ROOT=\"$(CLODEDIR)/\"
	LDLIBS += -lOpenCL
endif
ifeq ($(UNAME_S),Darwin)
	CPPFLAGS += -framework OpenCL
	LDFLAGS += -DCLODE_ROOT=\"$(CLODEDIR)/\"
	LDLIBS += -framework OpenCL
endif

all: testTrans testTraj testFeat

testTrans : $(OBJS1)
	$(CXX) $(LFLAGS) -o testTrans $(OBJS1) $(LDLIBS) #-framework OpenCL -DCLODE_ROOT=\"$(CLODEDIR)/\"
	
testTraj : $(OBJS2)
	$(CXX) $(LFLAGS) -o testTraj $(OBJS2) $(LDLIBS) # -framework OpenCL -DCLODE_ROOT=\"$(CLODEDIR)/\" 
	
testFeat : $(OBJS3)
	$(CXX) $(LFLAGS) -o testFeat $(OBJS3) $(LDLIBS) # -framework OpenCL -DCLODE_ROOT=\"$(CLODEDIR)/\"

testTransient.o: #testTransient.cpp OpenCLResource.hpp CLODE.hpp
	$(CXX) $(CPPFLAGS) testTransient.cpp -DCLODE_ROOT=\"$(CLODEDIR)/\"
	
testTrajectory.o: #testTrajectory.cpp OpenCLResource.hpp CLODE.hpp CLODEtrajectory.hpp
	$(CXX) $(CPPFLAGS) testTrajectory.cpp -DCLODE_ROOT=\"$(CLODEDIR)/\" 
	
testFeatures.o: #testFeatures.cpp OpenCLResource.hpp CLODE.hpp CLODEfeatures.hpp
	$(CXX) $(CPPFLAGS) testFeatures.cpp -DCLODE_ROOT=\"$(CLODEDIR)/\"

CLODE.o : #CLODE.cpp CLODE.hpp
	$(CXX) $(CPPFLAGS) ../src/CLODE.cpp #-DCLODE_ROOT=\"$(CLODEDIR)/\"
	
CLODEtrajectory.o : #CLODEtrajectory.cpp CLODEtrajectory.hpp
	$(CXX) $(CPPFLAGS) ../src/CLODEtrajectory.cpp #-DCLODE_ROOT=\"$(CLODEDIR)/\"
	
CLODEfeatures.o : #CLODEfeatures.cpp CLODEfeatures.hpp
	$(CXX) $(CPPFLAGS) ../src/CLODEfeatures.cpp #-DCLODE_ROOT=\"$(CLODEDIR)/\"

OpenCLResource.o : #OpenCLResource.cpp  OpenCLResource.hpp
	$(CXX) $(CPPFLAGS) ../src/OpenCLResource.cpp
	
.PHONY: clean
clean:
	\rm -f *.o testTrans testTraj testFeat
