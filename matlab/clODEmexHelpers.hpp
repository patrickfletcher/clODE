//helper functions for MEX interface to clODE

std::string getMatlabString(const mxArray *strData) {
    int strLength = (int)mxGetN(strData) + 1;
    std::vector<char> buf( strLength, 0 );

    mxGetString(strData, &buf[0], strLength);
    return std::string(&buf[0]);
}


cl_deviceType getDeviceTypeEnum(int deviceInt) {
	switch (deviceInt) {
		case 0:
			return CL_DEVICE_TYPE_DEFAULT;
		case 1:
			return CL_DEVICE_TYPE_CPU;
		case 2:
			return CL_DEVICE_TYPE_GPU;
		case 3:
			return CL_DEVICE_TYPE_ACCELERATOR;
		case 4:
			return CL_DEVICE_TYPE_ALL;
		default:
			return CL_DEVICE_TYPE_DEFAULT;
	}
}


SolverParams<double> getMatlabSPstruct(const mxArray *spptr) {
	SolverParams<double> sp;
	sp.dt=mxGetScalar( mxGetField(spptr,0,"dt") );
	sp.dtmax=mxGetScalar( mxGetField(spptr,0,"dtmax") );
	sp.abstol=mxGetScalar( mxGetField(spptr,0,"abstol") );
	sp.reltol=mxGetScalar( mxGetField(spptr,0,"reltol") );
	sp.max_steps=(int)mxGetScalar( mxGetField(spptr,0,"max_steps") );
	sp.max_store=(int)mxGetScalar( mxGetField(spptr,0,"max_store") );
	sp.nout=(int)mxGetScalar( mxGetField(spptr,0,"nout") );
	return sp;
}


ProblemInfo getMatlabProblemStruct(const mxArray *probptr) {
	std::string clFilename=getMatlabString( mxGetField(probptr,0,"clRHSfilename") );
	int nVar = mxGetScalar( mxGetField(probptr,0,"nVar") );
	int nPar = mxGetScalar( mxGetField(probptr,0,"nPar") );
	int nAux = mxGetScalar( mxGetField(probptr,0,"nAux") );
	int nWiener = mxGetScalar( mxGetField(probptr,0,"nWiener") );
	ProblemInfo newProblem;
	newProblem.clRHSfilename=clFilename;
	newProblem.nVar=nVar;
	newProblem.nPar=nPar;
	newProblem.nAux=nAux;
	newProblem.nWiener=nWiener;
	return newProblem;
}
