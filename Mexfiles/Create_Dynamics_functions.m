% Greybox function
mex CXXFLAGS='$CXXFLAGS $COMPFLAGS -std=c++17' CRMDYN_c.cpp "../src/*.cpp" "../numerical/*.cpp" -I../src/ -I../numerical/ -I../.. -I/usr/local/include/eigen3/

% Test function DEBUGGING
mex CXXFLAGS='$CXXFLAGS $COMPFLAGS -std=c++17' CRMDYN_c_mex.cpp "../src/*.cpp" "../numerical/*.cpp" -I../src/ -I../numerical/ -I../.. -I/usr/local/include/eigen3/