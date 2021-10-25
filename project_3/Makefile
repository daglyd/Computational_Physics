all: compile run python_ 

compile:  
	mkdir -p ./bin 

	g++ ./main_one_particle_exp.cpp src/Particle.cpp src/PenningTrap.cpp -I ./include -larmadillo -o bin/main_one_particle_exp.exe
	g++ ./main_two_particle_exp.cpp src/Particle.cpp src/PenningTrap.cpp -I ./include -larmadillo -o bin/main_two_particle_exp.exe
	g++ ./main_two_particle_exp_wo_interact.cpp src/Particle.cpp src/PenningTrap.cpp -I ./include -larmadillo -o bin/main_two_particle_exp_wo_interact.exe
run:
	mkdir -p ./Data
	mkdir -p ./Data/main_one_particle_exp
	./bin/main_one_particle_exp.exe 

	mkdir -p ./Data/main_two_particle_exp
	./bin/main_two_particle_exp.exe 
	./bin/main_two_particle_exp_wo_interact.exe

 
python_: 
	mkdir -p ./python 
	mkdir -p ./Figures/main_two_particle_exp
	mkdir -p ./Figures/main_one_particle_exp
	pip install -r requirements.txt 
	runipy One_particle.ipynb
	runipy Two_particle.ipynb 