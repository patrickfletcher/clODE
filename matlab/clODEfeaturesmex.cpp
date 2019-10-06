// mex interface to clODE class
// (c) Patrick Fletcher 2017
//
// based on: 
// class_wrapper_template.cpp
// Example of using a C++ class via a MEX-file
// by Jonathan Chappelow (chappjc)
  
#include "mex.h"

#include <vector>
#include <memory> //shared_ptr
#include <map>
#include <sstream>

////////////////////////  BEGIN Step 1: Configuration  ////////////////////////

#include "OpenCLResource.hpp"
#include "CLODE.hpp"
#include "CLODEfeatures.hpp"
#include "clODEmexHelpers.hpp"

typedef CLODEfeatures class_type;

// List actions
enum class Action
{
    New,
    Delete,
    SetNewProblem,
    SetStepper,
    SetPrecision,
    SetOpenCL,
    Initialize, //overridden for clODEfeatures
    SetProblemData,
    SetTspan,
    SetX0,
    SetPars,
    SetSolverPars,
    SeedRNG,
    Transient,
    UpdateTspan,
    UpdateX0,
    GetTspan,
    GetX0,
    GetXf,
    GetAuxf,
//from here are clODEfeatures derived actions
    SetObserverPars,
    SetObserver,
    Features,
    GetNFeatures,
    GetF
};

// Map string (first input argument to mexFunction) to an Action
const std::map<std::string, Action> actionTypeMap =
{
    { "new",            Action::New },
    { "delete",         Action::Delete },
    { "setnewproblem",  Action::SetNewProblem },
    { "setstepper",     Action::SetStepper },
    { "setprecision",   Action::SetPrecision },
    { "setopencl",      Action::SetOpenCL },
    { "initialize",     Action::Initialize },
    { "setproblemdata", Action::SetProblemData },
    { "settspan",       Action::SetTspan },
    { "setx0",          Action::SetX0 },
    { "setpars",        Action::SetPars },
    { "setsolverpars",  Action::SetSolverPars },
    { "seedrng",        Action::SeedRNG },
    { "transient",      Action::Transient },
    { "updatetspan",    Action::UpdateTspan },
    { "updatex0",       Action::UpdateX0 },
    { "gettspan",       Action::GetTspan },
    { "getx0",          Action::GetX0 },
    { "getxf",          Action::GetXf },
    { "getauxf",        Action::GetAuxf },
//from here are clODEfeatures derived actions
    { "setobserverpars", Action::SetObserverPars },
    { "setobserver",    Action::SetObserver },
    { "features",       Action::Features },
    { "getnfeatures",   Action::GetNFeatures },
    { "getf",        	Action::GetF }
}; 

//custom helper prototype
ObserverParams<double> getMatlabOPstruct(const mxArray *opptr);

/////////////////////////  END Step 1: Configuration  /////////////////////////

typedef unsigned int handle_type;
typedef std::pair<handle_type, std::shared_ptr<class_type>> indPtrPair_type; // or boost::shared_ptr
typedef std::map<indPtrPair_type::first_type, indPtrPair_type::second_type> instanceMap_type;
typedef indPtrPair_type::second_type instPtr_t;

// getHandle pulls the integer handle out of prhs[1]
handle_type getHandle(int nrhs, const mxArray *prhs[]);
// checkHandle gets the position in the instance table
instanceMap_type::const_iterator checkHandle(const instanceMap_type&, handle_type);


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // static storage duration object for table mapping handles to instances
    static instanceMap_type instanceTab;

    if (nrhs < 1 || !mxIsChar(prhs[0]))
        mexErrMsgTxt("First input must be an action string ('new', 'delete', or a method name).");

    //~ char *actionCstr = mxArrayToString(prhs[0]); // convert char16_t to char
    //~ std::string actionStr(actionCstr); mxFree(actionCstr);
	std::string actionStr=getMatlabString(prhs[0]);

    for (auto & c : actionStr) c = ::tolower(c); // remove this for case sensitivity

    if (actionTypeMap.count(actionStr) == 0)
        mexErrMsgTxt(("Unrecognized action (not in actionTypeMap): " + actionStr).c_str());

    // If action is not "new" or "delete" try to locate an existing instance based on input handle
    instPtr_t instance;
    if (actionTypeMap.at(actionStr) != Action::New && actionTypeMap.at(actionStr) != Action::Delete) {
        handle_type h = getHandle(nrhs, prhs);
        instanceMap_type::const_iterator instIt = checkHandle(instanceTab, h);
        instance = instIt->second;
    }

    //NOTE: if not 'new', the first two RHS args are actionStr, instanceHandle. User inputs start at 3rd arg
    //TODO: could wrap this somehow to make the code here more readable...
    
	//////// Step 2: customize the each action in the switch in mexFuction ////////
    switch (actionTypeMap.at(actionStr))
    {
    case Action::New:
    {
		
 	#if defined(WIN32)||defined(_WIN64)
 		_putenv_s("CUDA_CACHE_DISABLE", "1");
 	#else
 		setenv("CUDA_CACHE_DISABLE", "1", 1);
 	#endif
 	
		//sig: clODEobjective(nPar,nVar,nObjTimes,nObjVars,devicetype=all,vendor=any)
		
        handle_type newHandle = instanceTab.size() ? (instanceTab.rbegin())->first + 1 : 1;

		//create a new object
        std::pair<instanceMap_type::iterator, bool> insResult;
        
		//PARSE INPUT arguments ('new', problemInfoStruct, stepperInt, clSinglePrecisionBool, openclVendor=ANY, openclDeviceType=DEFAULT)
		
        if (nrhs < 4) {
			mexErrMsgTxt("Incorrect number of input arguments for clODEobjective object constructor");
		}
		
        std::string clFilename=getMatlabString( mxGetField(prhs[1],0,"clRHSfilename") );
        int nVar = mxGetScalar( mxGetField(prhs[1],0,"nVar") );
        int nPar = mxGetScalar( mxGetField(prhs[1],0,"nPar") );
        int nAux = mxGetScalar( mxGetField(prhs[1],0,"nAux") );
        int nWiener = mxGetScalar( mxGetField(prhs[1],0,"nWiener") );
        
        ProblemInfo newProblem;
        newProblem.clRHSfilename=clFilename;
        newProblem.nVar=nVar;
        newProblem.nPar=nPar;
        newProblem.nAux=nAux;
        newProblem.nWiener=nWiener;
        
        //correct integer must be supplied from MatLab
		int stepint=(int) mxGetScalar(prhs[2]);
        StepperType steppertype=(StepperType) stepint;
        
		bool clSinglePrecision=(bool) mxGetScalar(prhs[3]);
        
        //opencl device selection: force matlab user to use "vendor" and/or "devicetype", always pass in args for this
		// cl_vendor vendor = static_cast<cl_vendor>((int)mxGetScalar(prhs[4]));
		// cl_deviceType devicetype = getDeviceTypeEnum(static_cast<int>(mxGetScalar(prhs[5])) );	
		// OpenCLResource opencl(devicetype,vendor);

        //opencl device selection: assume the matlab caller selects by plaformID and deviceID, as returned by queryOpenCL 
        unsigned int platformID = (unsigned int)mxGetScalar(prhs[4]); 
        unsigned int deviceID = (unsigned int)mxGetScalar(prhs[5]);
		OpenCLResource opencl(platformID, deviceID);
		
        //correct integer must be supplied from MatLab
		int obsint=(int) mxGetScalar(prhs[6]);
        ObserverType observertype=(ObserverType) obsint;
        
			
		insResult = instanceTab.insert(indPtrPair_type(newHandle, std::make_shared<class_type>(newProblem,steppertype,observertype,clSinglePrecision, opencl)));

        if (!insResult.second) // sanity check
            mexPrintf("Oh, bad news.  Tried to add an existing handle."); // shouldn't ever happen
        else
            mexLock(); // add to the lock count

		// return the handle
        plhs[0] = mxCreateDoubleScalar(insResult.first->first); // == newHandle

        break;
    }
    case Action::Delete:
    { //rhs='delete',instanceID
        instanceMap_type::const_iterator instIt = checkHandle(instanceTab, getHandle(nrhs, prhs));
        instanceTab.erase(instIt);
        mexUnlock();
        plhs[0] = mxCreateLogicalScalar(instanceTab.empty()); // info
        break;
    }
    //base class CLODE methods (rhs 0='methodName', 1=instanceID, ...)
    case Action::SetNewProblem:
	{ //inputs: prob
        instance->setNewProblem(getMatlabProblemStruct(prhs[2]));
        break;
	}
    case Action::SetStepper:
	{ //inputs: steppertype int
		int stepint=(int) mxGetScalar(prhs[2]);
        StepperType steppertype=(StepperType) stepint;
        instance->setStepper(steppertype);
        break;
	}
    case Action::SetPrecision:
	{ //inputs: clSinglePrecision
        instance->setPrecision((bool) mxGetScalar(prhs[2]));
        break;
	}
    case Action::SetOpenCL:
	{ //inputs: vendor/devicetype
		cl_vendor vendor = static_cast<cl_vendor>((int)mxGetScalar(prhs[2]));
		cl_deviceType devicetype = getDeviceTypeEnum(static_cast<int>(mxGetScalar(prhs[3])) );	
		OpenCLResource opencl(devicetype,vendor);
        instance->setOpenCL(opencl);
        break;
	}
    case Action::SetProblemData:
	{ //inputs: x0, pars
        std::vector<double> x0( static_cast<double *>(mxGetData(prhs[2])),  static_cast<double *>(mxGetData(prhs[2])) + mxGetNumberOfElements(prhs[2]) ); 
        std::vector<double> pars(static_cast<double *>(mxGetData(prhs[3])),  static_cast<double *>(mxGetData(prhs[3])) + mxGetNumberOfElements(prhs[3]) );  		
        instance->setProblemData(x0, pars);
        break;
	}
    case Action::SetTspan:
	{ //inputs: tspan
        std::vector<double> tspan ( static_cast<double *>(mxGetData(prhs[2])),  static_cast<double *>(mxGetData(prhs[2])) + mxGetNumberOfElements(prhs[2]) ); 
        instance->setTspan(tspan);
        break;
	}
    case Action::SetX0:
	{ //inputs: x0
        std::vector<double> x0( static_cast<double *>(mxGetData(prhs[2])),  static_cast<double *>(mxGetData(prhs[2])) + mxGetNumberOfElements(prhs[2]) );  
        instance->setX0(x0);
        break;
	}
    case Action::SetPars:
	{ //inputs: pars
        std::vector<double> pars( static_cast<double *>(mxGetData(prhs[2])),  static_cast<double *>(mxGetData(prhs[2])) + mxGetNumberOfElements(prhs[2]) );  
        instance->setPars(pars);
        break;
	}
    case Action::SetSolverPars:
	{	//inputs: sp 
		SolverParams<double> sp = getMatlabSPstruct(prhs[2]);
        instance->setSolverParams(sp);
        break;
	}
    case Action::SeedRNG:
	{ //inputs: none, or mySeedInt
		if (nrhs==2) 
			instance->seedRNG();
		else if (nrhs==3)
			instance->seedRNG((int)mxGetScalar(prhs[2]));

        break;
	}
    case Action::Transient:
	{
        instance->transient();
        break;
	}
    case Action::UpdateTspan:
	{
        instance->updateTspan();

        break;
	}
    case Action::UpdateX0:
	{
        instance->updateX0();

        break;
	}
    case Action::GetTspan:
    {
        std::vector<double> tspan=instance->getTspan();
		plhs[0]=mxCreateDoubleMatrix(tspan.size(), 1, mxREAL);
        std::copy(tspan.begin(), tspan.end(), (double *)mxGetData(plhs[0]));
        break;
	}
    case Action::GetX0:
    {
        std::vector<double> x0=instance->getX0();
		plhs[0]=mxCreateDoubleMatrix(x0.size(), 1, mxREAL);
        std::copy(x0.begin(), x0.end(), (double *)mxGetData(plhs[0]));
        break;
	}
    case Action::GetXf:
    {
        std::vector<double> xf=instance->getXf();
		plhs[0]=mxCreateDoubleMatrix(1, xf.size(), mxREAL);
        std::copy(xf.begin(), xf.end(), (double *)mxGetData(plhs[0]));
        break;
	}
    case Action::GetAuxf:
    {
        std::vector<double> auxf=instance->getAuxf();
		plhs[0]=mxCreateDoubleMatrix(auxf.size(), 1, mxREAL);
        std::copy(auxf.begin(), auxf.end(), (double *)mxGetData(plhs[0]));
        break;
	}
	//CLODEfeatures methods:
    case Action::Initialize:
	{ //inputs: tspan, x0, pars, sp, op //OVERRIDE
        std::vector<double> tspan ( static_cast<double *>(mxGetData(prhs[2])),  static_cast<double *>(mxGetData(prhs[2])) + mxGetNumberOfElements(prhs[2]) ); 
        std::vector<double> x0 ( static_cast<double *>(mxGetData(prhs[3])),  static_cast<double *>(mxGetData(prhs[3])) + mxGetNumberOfElements(prhs[3]) ); 
        std::vector<double> pars (static_cast<double *>(mxGetData(prhs[4])),  static_cast<double *>(mxGetData(prhs[4])) + mxGetNumberOfElements(prhs[4]) );  
		SolverParams<double> sp = getMatlabSPstruct(prhs[5]);
		ObserverParams<double> op = getMatlabOPstruct(prhs[6]);
        instance->initialize(tspan, x0, pars, sp, op);       
        break;
	}
    case Action::SetObserverPars:
	{	//inputs: op 
		ObserverParams<double> op = getMatlabOPstruct(prhs[2]);
        instance->setObserverParams(op);
        break;
	}
    case Action::SetObserver:
    {   //input int corresponding to an observer type
		int obsint=(int) mxGetScalar(prhs[2]);
        ObserverType observertype=(ObserverType) obsint;
        instance->setObserver(observertype);
        break;
	}
    case Action::Features:
    {   //inputs: none, or int - 1=do init, 0=don't init
		if (nrhs==2) { 
			instance->features(); }
		else if (nrhs==3) {
			int doInit=(int)mxGetScalar(prhs[2]);
			if (doInit<0 || doInit>1) {mexErrMsgTxt("Argument must be 0 or 1 for features with observer data initialization control");}
			instance->features(doInit);
		}
        break;
	}
    case Action::GetNFeatures:
    {  
		plhs[0]=mxCreateDoubleScalar(instance->getNFeatures());
        break;
	}
    case Action::GetF:
    {  
        std::vector<cl_double> F=instance->getF();
		plhs[0]=mxCreateDoubleMatrix(F.size(), 1, mxREAL);
        std::copy(F.begin(), F.end(), (double *)mxGetData(plhs[0]));
        break;
	}
    default:
        mexErrMsgTxt(("Unhandled action: " + actionStr).c_str());
        break;
    }
    ////////////////////////////////  DONE!  ////////////////////////////////
}

handle_type getHandle(int nrhs, const mxArray *prhs[])
{
    if (nrhs < 2 || mxGetNumberOfElements(prhs[1]) != 1) // mxIsScalar in R2015a+
        mexErrMsgTxt("Specify an instance with an integer handle.");
    return static_cast<handle_type>(mxGetScalar(prhs[1]));
}

instanceMap_type::const_iterator checkHandle(const instanceMap_type& m, handle_type h)
{
    auto it = m.find(h);

    if (it == m.end()) {
        std::stringstream ss; ss << "No instance corresponding to handle " << h << " found.";
        mexErrMsgTxt(ss.str().c_str());
    }

    return it;
}
	


ObserverParams<double> getMatlabOPstruct(const mxArray *opptr) {
	ObserverParams<double> op;
	op.eVarIx=mxGetScalar( mxGetField(opptr,0,"eVarIx") ) -1; //convert matlab 1-based index to C 0-based.
	op.fVarIx=mxGetScalar( mxGetField(opptr,0,"fVarIx") ) -1;
	op.maxEventCount=mxGetScalar( mxGetField(opptr,0,"maxEventCount") );
	op.minXamp=mxGetScalar( mxGetField(opptr,0,"minXamp") );
	op.minIMI=mxGetScalar( mxGetField(opptr,0,"minIMI") );
	op.nHoodRadius=mxGetScalar( mxGetField(opptr,0,"nHoodRadius") );
	op.xUpThresh=mxGetScalar( mxGetField(opptr,0,"xUpThresh") );
	op.xDownThresh=mxGetScalar( mxGetField(opptr,0,"xDownThresh") );
	op.dxUpThresh=mxGetScalar( mxGetField(opptr,0,"dxUpThresh") );
	op.dxDownThresh=mxGetScalar( mxGetField(opptr,0,"dxDownThresh") );
	op.eps_dx=mxGetScalar( mxGetField(opptr,0,"eps_dx") );
	return op;
}
