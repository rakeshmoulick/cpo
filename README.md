# cpo
Cold Plasma Oscillation (Normalized and Un-Normalized Codes along with Python Plotting Files)

# Run Instructions
1. Go to the folder in your local computer where you wish to clone the code from GitHub.
2. run > git clone https://github.com/rakeshmoulick/cpo.git
3. run > cd cpo
4. run > make veryclean
5. run > make all
6. run > ./epic input.ini
7. This should run the code in C++ (src/main.cpp).
8. Once it is completed, there will be a data folder created. Some of the files will be under ../data/files. Particle data will be saved in /data folder itself. 

# Plotting Data with Python
1. run > cd python_scripts
2. running instructions are given in every python file (beginning of the file). 
3. An example run > python3 pe_semilogy.py ../input.ini for the energy conservation plot.
