#include "TaskQueue.h"

TaskQueue::TaskQueue()
{
  Tail = NULL;
  Head = NULL;
}

TaskQueue::~TaskQueue()
{
  NodePtr temp = Head;
  NodePtr CurrentPtr = Head;

  while (CurrentPtr != NULL) {
    CurrentPtr = CurrentPtr->Next;
    delete temp;
    temp = CurrentPtr;
  }
}

void TaskQueue::AddTask(Task* T)
{
  if (Tail==NULL) {
    Tail = new Node(T);
    Head = Tail;
  }
  else {
    Tail->Next = new Node(T);
    Tail = Tail->Next;
  }
}

Task* TaskQueue::Extract()
{
  Task* T;
  if (isEmpty()) T=NULL;
  else {
    NodePtr temp = Head;
    T = temp->Data;
    Head = temp->Next;
    delete temp;
    if (Head==NULL) Tail=NULL;
  }
  return T;
}

