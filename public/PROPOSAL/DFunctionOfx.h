/**
 * function and its derivative must be included in a class that implements this interface
 */

#ifndef DFunctionOfx_H_
#define DFunctionOfx_H_

#include "PROPOSAL/FunctionOfx.h"

class DFunctionOfx : public  FunctionOfx
{
public:

    virtual double dFunction(double x){return x;}

//----------------------------------------------------------------------------//

    virtual ~DFunctionOfx()
    {

    }

};

#endif // DFunctionOfx_H

