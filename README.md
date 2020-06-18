# GB-500_Trabalho
Codigos do Trabalho da disciplina GB-500 Programação em Criptografia


Como compilar o arquivo initial_condition.c

gcc -c initial_condition.c -lgmp -lm -o initial_condition.o && gcc -o initial_condition initial_condition.o utils.o Pailler.o -lgmp -lm

Executar initial_condition

./initial_condition

Como compilar o arquivo principal

gcc -c trabalho.c -o trabalho.o -lgmp -lm && gcc -c simulation.c -lgmp -lm -o simulation.o && gcc -o trabalho simulation.o trabalho.o utils.o Pailler.o -lgmp -lm

Executar o arquivo principal

./trabalho
