CPP=g++
CC=gcc
EXEC= chaostest

all: ${EXEC}
chaostest: chaostest.cc
	${CPP} -g -Wall $^ -o $@

clean:
	-rm ${EXEC}
