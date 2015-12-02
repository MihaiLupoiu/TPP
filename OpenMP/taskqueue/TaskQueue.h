#include <cstdio>

class Task {
public:
  virtual void execute() = 0;
};

struct Node {
  Task* Data;
  Node* Next;
  Node(Task* T): Data(T) { Next=NULL; }
};
typedef struct Node* NodePtr;

class TaskQueue
{
private:
  NodePtr Head, Tail;
public:
  TaskQueue();
  ~TaskQueue();
  void AddTask(Task*);
  Task* Extract();
  bool isEmpty() { return Head==NULL; }
};

