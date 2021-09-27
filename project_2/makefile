all: compile run
compile:
	g++ main.cpp src/matrix.cpp src/jacobi.cpp src/problems.cpp -larmadillo -o main.exe 

run:
	mkdir -p ./Data
	mkdir -p ./Data/problem_5
	mkdir -p ./Data/problem_6
	mkdir -p ./Data/problem_7
	mkdir -p ./Figures
	pip install -r python/requirements.txt
	./main.exe 
	runipy python/problem_6.ipynb
	runipy python/problem_7.ipynb