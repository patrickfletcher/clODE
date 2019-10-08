// MEX function to wrap the OpenCLResource function of the same name. 
// Either print to console (no lhs) or return array of platform info structs
// Add rhs args to filter?

#include "mex.h"
#include "OpenCLResource.hpp"

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
{
    if (nlhs==0)
    {
        printOpenCL();
        return;
    }
    else if (nlhs>1)
    {
        mexErrMsgTxt("More than one output argument not supported");
    }

    //get the platform/device info
    std::vector<platformInfo> pinfo = queryOpenCL(); 

    //convert this info to a Matlab struct array, one element per device
    unsigned int nElem=0;
    for(unsigned int i = 0; i < pinfo.size(); ++i) {
        nElem+=pinfo[i].nDevices;
    }
    mwSize dims[2]={1, nElem};

    //TODO: add version and platform stuff too
    const char *field_names[] = {"platformID","deviceID","name", "type","vendor","computeUnits","maxClock","memSize","doubleSupport","available"};
    int nFields = sizeof(field_names) / sizeof(*field_names);

    plhs[0] = mxCreateStructArray(2,dims,nFields,field_names);

    int ix=0;
    for(unsigned int i = 0; i < pinfo.size(); ++i) {
        for (unsigned int j=0; j<pinfo[i].nDevices; j++){
            deviceInfo dinfo=pinfo[i].device_info[j];

            mxSetField(plhs[0],ix,"platformID",mxCreateDoubleScalar(i));
            mxSetField(plhs[0],ix,"deviceID",mxCreateDoubleScalar(j));
            mxSetField(plhs[0],ix,"name",mxCreateString(dinfo.name.c_str()));
            mxSetField(plhs[0],ix,"type",mxCreateString(dinfo.devTypeStr.c_str()));
            mxSetField(plhs[0],ix,"vendor",mxCreateString(dinfo.vendor.c_str()));
            mxSetField(plhs[0],ix,"computeUnits",mxCreateDoubleScalar(dinfo.computeUnits));
            mxSetField(plhs[0],ix,"maxClock",mxCreateDoubleScalar(dinfo.maxClock));
            mxSetField(plhs[0],ix,"memSize",mxCreateDoubleScalar(dinfo.deviceMemSize/1024/1024));
            mxSetField(plhs[0],ix,"doubleSupport",mxCreateDoubleScalar(dinfo.doubleSupport));
            mxSetField(plhs[0],ix,"available",mxCreateDoubleScalar(dinfo.deviceAvailable));

            ix++;
        }
    }

    return;
}