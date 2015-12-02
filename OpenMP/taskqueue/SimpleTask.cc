#include "TaskQueue.h"
#include <iostream>
using namespace std;

class SimpleTask: public Task {
private:
  int N;
public:
  SimpleTask(int n) {
    N=n;
  }
  void execute() {
    cout << "Task n = " << N << endl;
  }
};

int main(int argc, char **argv)
{ 
  TaskQueue queue;
  SimpleTask ST1(1);
  SimpleTask ST2(2);
  SimpleTask ST3(3);
  
  queue.AddTask(&ST1);
  queue.AddTask(&ST2);
  queue.AddTask(&ST3);

  Task* T;
  while ((T = queue.Extract()) != NULL) {
    T->execute();
  }

  return 0;
}
