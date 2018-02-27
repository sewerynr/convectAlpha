#ifndef BASICEXCEPTIONS_H
#define BASICEXCEPTIONS_H

#include <iostream>
#include <string>

class BasicException
{
public:
    std::string wyjatek;
    BasicException(const std::string& text):wyjatek(text) {}

};

class FileException : public BasicException
{
public:
    FileException() : BasicException("FERR: File can not be opened !!!")
    {}
};


#endif // BASICEXCEPTIONS_H
