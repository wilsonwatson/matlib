

all :
	g++ ${shell find . -name '*.cpp'} -o out -std=c++17 -lpthread -g
