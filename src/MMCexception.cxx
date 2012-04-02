#include "MMCexception.h"

using namespace std;

MMCexception::MMCexception(int errorNr, string object, string filename){

    errorNr_ = errorNr;
    object_ = object;
    filename_ = filename;

    stringstream infostream;

    infostream<<"In class "<<object_<<":\t Error Number "<<errorNr_<<":\t";
    info_=infostream.str();
}


const char* MMCexception::what() const throw(){

    stringstream error_information;

    switch(errorNr_){

    case 0:
        error_information<<info_<<"Cannot redirect stderr to the file "<<filename_;
        break;
    default:
        error_information<<info_<<"Unkown error";
        break;

    }
    return error_information.str().c_str();
}


