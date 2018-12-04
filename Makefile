all: main.o input.o initial_curves.o curve_similarity.o grid.o hash.o range.o initialization.o traversal_computation.o assignment.o pam.o mean_frechet_tree.o ending.o silhouette.o output_file.o
	gcc -o cluster main.o input.o initial_curves.o curve_similarity.o grid.o hash.o range.o initialization.o traversal_computation.o assignment.o pam.o  mean_frechet_tree.o ending.o silhouette.o output_file.o -g -lm

main.o: main.c
	gcc -c -g main.c

input.o: input.c
	gcc -c -g input.c

initial_curves.o: initial_curves.c
	gcc -c -g initial_curves.c

curve_similarity.o: curve_similarity.c
	gcc -c -g -lm curve_similarity.c

grid.o: grid.c
	gcc -c -g -lm grid.c

hash.o: hash.c
	gcc -c -g -lm hash.c

range.o: range.c
	gcc -c -g -lm range.c

initialization.o: initialization.c
	gcc -c -g -lm initialization.c

traversal_computation.o: traversal_computation.c
	gcc -c -g -lm traversal_computation.c

assignment.o: assignment.c
	gcc -c -g -lm assignment.c

pam.o: pam.c
	gcc -c -g pam.c

mean_frechet_tree.o: mean_frechet_tree.c
	gcc -c -g mean_frechet_tree.c

ending.o: ending.c
	gcc -c -g ending.c

silhouette.o: silhouette.c
	gcc -c -g silhouette.c

output_file.o: output_file.c
	gcc -c -g output_file.c

clean:
	rm cluster main.o input.o initial_curves.o curve_similarity.o initialization.o traversal_computation.o grid.o hash.o range.o assignment.o pam.o mean_frechet_tree.o ending.o silhouette.o output_file.o
