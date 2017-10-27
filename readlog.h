#ifndef READLOG_H
#define READLOG_H
#include <string>
#include "traiettoria.h"
#include <vector>


class ReadLog
{
public:
    ReadLog(std::string filename,Traiettoria * t=0,unsigned int skip=1);
    ~ReadLog();
    double * line(unsigned int index);
    unsigned int n_timestep(){return data.size();}
    unsigned int n_data(){return headers.size()-1;}
    unsigned int timestep(unsigned int index);
private:
    Traiettoria * traiettoria;
    std::vector<std::string> headers;
    std::vector<std::pair<unsigned int ,double*> > data;
    unsigned int skip,size,step_index;
    bool if_only_numbers(std::string str);

};

#endif // READLOG_H
