/**
 * function to be integrated must be included in a class that implements this interface.


  Every class which might be integrated needs a function like this. This is the function which will be derived by using the Rombergmethod. Be careful which function will be integrated because each class got its own.
 */

#ifndef FunctionOfx_H
#define FunctionOfx_H

#include <stdio.h>
#include <cmath>
#include <iostream>

class FunctionOfx
{
	
public:

    virtual double function(double x)
    {
        return x;
    }
//----------------------------------------------------------------------------//
    virtual ~FunctionOfx(){}

};

#endif //FunctionOfx_H
