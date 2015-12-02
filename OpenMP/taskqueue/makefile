CXX      = c++
CXXFLAGS = -O3
LDFLAGS  =
LIBS     = -lm

OBJ_SIMPLE = TaskQueue.o SimpleTask.o
tq_simple: $(OBJ_SIMPLE)
	$(CXX) -o $@ $(OBJ_SIMPLE) $(LDFLAGS)

OBJ_MANDEL = TaskQueue.o MandelTask.o misc.o
tq_mandel: $(OBJ_MANDEL)
	$(CXX) -o $@ $(OBJ_MANDEL) $(LDFLAGS) $(LIBS)

clean:
	@rm -f $(OBJ_SIMPLE) $(OBJ_MANDEL)
