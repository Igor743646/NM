/*
    Modified Fractional Differential Equations Solver 
*/

#pragma once

#include <functional>

#include "matrix2D.h"
#include "FDESbase.h"

class MFDES : public FDESbase {

public: 

    MFDES(const FDESbase&);

    void solve();
};