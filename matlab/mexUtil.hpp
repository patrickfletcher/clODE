
//some utils to help mex file construction

std::string getMatlabString(const mxArray *strData) {
    int strLength = mxGetN(strData) + 1;
    std::vector<char> buf( strLength, 0 );

    mxGetString(strData, &buf[0], strLength);
    return std::string(&buf[0]);
}
