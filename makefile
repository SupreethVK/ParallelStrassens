all: parastrasens.out strassens.out

parastrassens.out: parallelstrassens.c
	gcc -o parastrassens.out parallelstrassens.c -pthread
	
strassens.out: strassens.c
	gcc -o strassens.out strassens.c 
