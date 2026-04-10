#include "biolcccexception.h"
#include <string.h>

namespace BioLCCC
{

BioLCCCException::BioLCCCException(std::string message)
{
    mMessage = message;
};

BioLCCCException::~BioLCCCException() noexcept {};

const char* BioLCCCException::what() const noexcept
{
    return mMessage.c_str();
};
}

