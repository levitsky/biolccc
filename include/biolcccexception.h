#ifndef BIOLCCCEXCEPTION_H
#define BIOLCCCEXCEPTION_H

#include <string>

namespace BioLCCC
{

class BioLCCCException : public std::exception
{
public:
    BioLCCCException(std::string message);
    ~BioLCCCException() throw();
    virtual const char* what() const throw();
private:
    std::string mMessage;
};
}

#endif
