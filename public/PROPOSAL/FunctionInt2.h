/**
 * function to be interpolated must be included in a class that implements this interface
 */

#ifndef FunctionInt2_H
#define FunctionInt2_H

class FunctionInt2
{
public:
//----------------------------------------------------------------------------//
    virtual double functionInt(double x1, double x2)
    {
        return x1;
    }
//----------------------------------------------------------------------------//

    virtual ~FunctionInt2(){}

    		
};


#endif // FunctionInt2_H
