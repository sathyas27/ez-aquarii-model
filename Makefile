all: ac-system ez-aquarii

ac-system:
	gcc -o ac_system ac_system.c -lm

ez-aquarii:
	gcc -o ez_aquarii ez_aquarii.c -lm

clean:
	rm *.gif
	rm *.txt
	rm ac_system
	rm ez_aquarii