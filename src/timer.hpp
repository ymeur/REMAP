#ifndef __TIMER_HPP__
#define __TIMER_HPP__

#include "mpi.h" // some compiler fails if mpi.h is not the first header to be included
#include <string>
#include <map>

  using namespace std;
  
  
  class CTimer
  {
    public :
    
    double cumulatedTime;
    double lastTime;
    bool suspended;
    string name;
    
    CTimer(const string& name);
    void suspend(void);
    void resume(void);
    void reset(void);
    double getCumulatedTime(void);
    void print(void);
    static map<string,CTimer*> allTimer;
    static double getTime(void);
    static CTimer& get(string name);
  };




#endif
