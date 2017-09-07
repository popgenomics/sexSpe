# sexSpe  
## Compilation  
g++ sexSpe.cpp -std=c++17 -O3 -o sexSpe  
  
## execution  
./sexSpe 1000 10 10 500 50  

### arguments   
**arg1**: number of diploid individuals in the population (__integer__)  
**arg2**: number of loci under **S**ex **A**ntagonistic selection (__integer__)  
**arg3**: number of loci under **S**exually **C**oncordant selection (__integer__)  
**arg4**: number of generations (__integer__)  
**arg5**: number of reproductive males (__integer__)   
  
### TODO  
At that time, mutation rate is hard-coded and equal to **nLocus x 0.00001**  
where __nLocus__ is equal to: 1 (modifier of recombination) + 1 (modifier of expression) + 2.(1 + SA + SC)  
Inversion rate is also hard-coded and fixed to 0.00001    

