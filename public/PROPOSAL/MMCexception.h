#ifndef MMCEXCEPTION_H
#define MMCEXCEPTION_H

#include <exception>
#include <string>
#include <sstream>

class MMCexception : public std::exception {

public:

    MMCexception(int errorNr, std::string object, std::string filename);

//----------------------------------------------------------------------------//

    virtual ~MMCexception() throw() {}

//----------------------------------------------------------------------------//

    virtual const char* what() const throw();

private:

    int errorNr_;
    std::string object_;
    std::string filename_;
    std::string info_;

};


#endif // MMCEXCEPTION_H
