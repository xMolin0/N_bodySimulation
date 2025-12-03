#ifndef SEQTASKING__HPP
#define SEQTASKING__HPP


#include <functional>


namespace tasking {

  void taskwait() {

  }
  
  void taskstart(std::function<void()> f)
  {
    f();
  }

  void doinparallel(std::function<void()> f, int nbthreads = 1)
  {
    f();
  }
  
};

#endif
