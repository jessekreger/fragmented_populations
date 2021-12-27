// Example simulation code for Figure 8 in main text
// Written in C++ for speed (to perform many simulations)
// Compiled and run using Xcode version 12 on a Mac

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <random>
#include <math.h>
#include <iomanip>

using namespace std;

double choose(double Nint, double Kint) { // choose function
    double result = 1;
    for (int iii = 1; iii <= Kint; iii++)
    {
        result = result * (Nint - (Kint - iii));
        result = result / iii;
    }
    return result;
}

int main(int argc, const char * argv[]) { // main function
    
    FILE *file_handle;
    
    srand(time(0)); // seed random number generators
    srand48(time(0));
    std::mt19937_64 mt_rand(time(0));
    std::default_random_engine generator(time(0));
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    
    double Ntime = pow(10,9);
    int Ntime_int = Ntime;
    double mutationrate = pow(10,-4);
    double mutationrateback = mutationrate;
    double fitness = 0.95;
    double NC = log((((-1 + fitness)*fitness)/(-1 + fitness + mutationrate) - fitness*mutationrateback/mutationrate))/log(fitness); // calculate NC
    //printf("NC = %lf \n", NC);
    
    std::string filename = "file_name.csv";
    
    double Nmig;
    double Migset;
    double mig;
    double rr;
    double rr2;
    double Pdivm;
    double Pdiem;
    double Pup;
    double Pdown;
    int i;
    int j;
    int ii;
    int jj;
    int t;
    int k;
    int i1;
    int i2;
    double mut1;
    double mut2;
    double dif;
    double standarddev;
    double average;
    double standarderr;
    int h;
    double randomprob;
    int h2;
    double yo1;
    double yo2;
    double totalmutants;
    double Ncells;
    
    for (i=1; i<51; i++) { // define which Ncells to run and define other parameters
        Ncells = 5*i+15;
        
        double Npatches = 20;
        int Npatches_int = Npatches;
        
        int Ncells_int = Ncells;
        
        Migset = Ncells/5;
        Nmig = Npatches/5;
        Ncells_int = Ncells;
        Npatches_int = Npatches;
        
        double hg[Ncells_int+1];
        double hgsum[Ncells_int+1];
        double hg2[Ncells_int+1];
        double hgsum2[Ncells_int+1];
        
        double xold[Npatches_int];
        double xnew[Npatches_int];
        
        double nom;
        double jjj;
        
        double smbalance = (1/(2* (1 - fitness))*(1 - fitness + mutationrate +fitness* mutationrate -sqrt(pow((-1 + mutationrate),2.0) + pow(fitness,2.0)* pow((-1 + mutationrate),2.0) +2 *fitness *(-1 + 2 *mutationrate + pow(mutationrate,2.0)))));
        //printf("noc = %d, smbalance = %lf\n", Ncells_int, smbalance);
        
        double threshold = pow(10.0,-5.0); // used for adaptive migration rate (to approximate where the average number of mutants = 2*smbalance)
        double average2 = 0;
        jjj = 0;
        while (abs(average2-2*smbalance) > threshold) {
            jjj = jjj+1;
            if (jjj == 1) {
                mig = 0.01;
            }
            else {
                if (average2-2*smbalance > 0) {
                    mig = mig + pow(10,-0.1*jjj);
                    if (mig > 1)
                    {
                        mig = 0.99;
                    }
                }
                else {
                    mig = mig - pow(10,-0.1*jjj);
                    if (mig < 0)
                    {
                        mig = 0.000001;
                    }
                }
            }
            jjj = jjj+1;
            
            
            
            for (ii=0; ii<Npatches; ii++) {
                xold[ii] = 0;
            }
             

            totalmutants = 0;
            standarddev = 0;
            
            for (t=0; t<Ntime; t++) { // perform simulations
                for (ii=0; ii<Npatches; ii++) {
                    xnew[ii] = xold[ii];
                }
                for (ii=0; ii<Npatches; ii++) { // Moran step
                    Pdivm = ((1-mutationrate)*fitness*xold[ii]+mutationrate*(Ncells-xold[ii]))/(fitness*xold[ii]+Ncells-xold[ii]);
                    Pdiem = xold[ii]/Ncells;
                    Pup = Pdivm*(1-Pdiem);
                    Pdown = (1-Pdivm)*Pdiem;
            
                    rr = distribution(generator);
                    
                    if (rr < Pup) {
                        xnew[ii] = xnew[ii]+1;
                    }
                    else if (rr < Pup + Pdown) {
                        xnew[ii] = xnew[ii]-1;
                    }
                } // end of Npatches loop (ii)
            
                rr2 = distribution(generator);
                if (rr2 < mig) { // migration step
                    for (k=0; k<Nmig; k++) {
                    
                        i1 = mt_rand() % Npatches_int;
                
                        
                        for (h=0; h<Ncells+1; h++) {
                            hg[h] = 0;
                        }
                        for (h=0; h<Migset+1; h++) {
                            hg[h] = choose(xnew[i1],h)*choose((Ncells-xnew[i1]),Migset-h)/choose(Ncells,Migset);
                        }
                        hgsum[0] = hg[0];
                        for (h=1; h<Ncells+1; h++) {
                            hgsum[h] = hgsum[h-1]+hg[h];
                        }
                        randomprob = distribution(generator);
                        for (h=0; h<Ncells+1; h++) {
                            if (randomprob < hgsum[h]) {
                                mut1 = h;
                                break;
                            }
                        }
                        
                        i2 = mt_rand() % Npatches_int;
                       
                        
                        for (h2=0; h2<Ncells+1; h2++) {
                            hg2[h2] = 0;
                        }
                        for (h2=0; h2<Migset+1; h2++) {
                            hg2[h2] = choose(xnew[i2],h2)*choose((Ncells-xnew[i2]),Migset-h2)/choose(Ncells,Migset);
                        }
                        hgsum2[0] = hg2[0];
                        for (h2=1; h2<Ncells+1; h2++) {
                            hgsum2[h2] = hgsum2[h2-1]+hg2[h2];
                        }
                        randomprob = distribution(generator);
                        for (h2=0; h2<Ncells+1; h2++) {
                            if (randomprob < hgsum2[h2]) {
                                mut2 = h2;
                                break;
                            }
                        }
                        
                        yo1 = xnew[i1];
                        yo2 = xnew[i2];
                        mut1 = min(mut1,yo1);
                        mut2 = min(mut2,yo2);
                        dif = mut1 - mut2;
                        xnew[i1] = xnew[i1] - dif;
                        xnew[i2] = xnew[i2] + dif;
                    } // end of Nmig loop (k)
                
                        
                }
                nom = 0;
                for (ii=0; ii<Npatches; ii++) {
                    xold[ii] = xnew[ii];
                    totalmutants = totalmutants+xold[ii];
                    nom = nom+xold[ii];
                }
                standarddev = standarddev+pow(nom,2.0);
            } // end of time loop (t)
            
        
            average = totalmutants/Ntime;
            average2 = totalmutants/Ntime/Npatches/Ncells;
            standarddev = standarddev-(Ntime*pow(average,2.0));
            standarddev = standarddev/(Ntime-1);
            standarddev = sqrt(standarddev);
                
            standarderr = standarddev/sqrt(Ntime);
            
                
            file_handle = fopen(filename.c_str(), "a"); // if not working then change your permissions or add full file path
            fprintf(file_handle, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", Npatches, Ncells, mig, average, standarderr, Ntime, NC, mutationrate, fitness, standarddev);
            fclose(file_handle);
            
        } // end of mig loop (j)
    } // end of Ncells loop (i)
    

    return 0;
}

