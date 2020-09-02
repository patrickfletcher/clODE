//helper functions for MEX interface to clODE

#ifndef CLODE_MEX_HELPERS_H
#define CLODE_MEX_HELPERS_H

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


SolverParams<cl_double> getMatlabSPstruct(const mxArray *spptr) {
	SolverParams<cl_double> sp;
	sp.dt=mxGetScalar( mxGetField(spptr,0,"dt") );
	sp.dtmax=mxGetScalar( mxGetField(spptr,0,"dtmax") );
	sp.abstol=mxGetScalar( mxGetField(spptr,0,"abstol") );
	sp.reltol=mxGetScalar( mxGetField(spptr,0,"reltol") );
	sp.max_steps=(cl_int)mxGetScalar( mxGetField(spptr,0,"max_steps") );
	sp.max_store=(cl_int)mxGetScalar( mxGetField(spptr,0,"max_store") );
	sp.nout=(cl_int)mxGetScalar( mxGetField(spptr,0,"nout") );
	return sp;
}


ProblemInfo getMatlabProblemStruct(const mxArray *probptr) {
	ProblemInfo newProblem;
	newProblem.clRHSfilename=getMatlabString( mxGetField(probptr,0,"clRHSfilename") );
	newProblem.nVar=(cl_int)mxGetScalar( mxGetField(probptr,0,"nVar") );
	newProblem.nPar=(cl_int)mxGetScalar( mxGetField(probptr,0,"nPar") );
	newProblem.nAux=(cl_int)mxGetScalar( mxGetField(probptr,0,"nAux") );
	newProblem.nWiener=(cl_int)mxGetScalar( mxGetField(probptr,0,"nWiener") );

	const mxArray *namesPtr;
	mwSize nNames;
	namesPtr = mxGetField(probptr,0,"varNames");
	nNames = mxGetNumberOfElements(namesPtr); 
	for (mwIndex i=0; i<nNames; i++)
		newProblem.varNames.push_back(mxArrayToString(mxGetCell(namesPtr, i)));

	namesPtr = mxGetField(probptr,0,"parNames");
	nNames = mxGetNumberOfElements(namesPtr); 
	for (mwIndex i=0; i<nNames; i++)
		newProblem.parNames.push_back(mxArrayToString(mxGetCell(namesPtr, i)));

	namesPtr = mxGetField(probptr,0,"auxNames");
	nNames = mxGetNumberOfElements(namesPtr); 
	for (mwIndex i=0; i<nNames; i++)
		newProblem.auxNames.push_back(mxArrayToString(mxGetCell(namesPtr, i)));

	// for (auto n : newProblem.varNames)
	// 	printf("%s\n",n);


	return newProblem;
}

#endif // CLODE_MEX_HELPERS_H