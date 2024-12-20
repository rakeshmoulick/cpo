/*
 ****************************************************
 Developer:
 DR. RAKESH MOULICK, CPP-IPR, ASSAM, INDIA
 ****************************************************
 */

# include <iostream>
# include <cmath>
# include <cstdlib>
# include <vector>
# include <list>
# include <ctime>
# include <random>
# include <cstring>
# include <fstream>
# include <sstream>
# include <chrono>
extern "C"
{
	# include "iniparser.h"

}
using namespace std;

/***************** INIPARSER ***************/
int  parse_ini_file(char * ini_name);

/* Random Number Generator */
std::mt19937 mt_gen(0);
std::uniform_real_distribution<double> rnd_dist(0,1.0);
double rnd()
{
    return rnd_dist(mt_gen);
}

/* Define universal constants */
double EPS;
double K;
double massE;
double chargeE;
double eV;
double e;
double AMU;
double massI;
double EV_TO_K;
double pi;

/* Define Simulation Parameters*/
double density;       // Plasma Density
double stepSize;      // Cell Spacing
double DT;            // Time steps
double DT_coeff;

double tempE;         // Temperature of the cold electrons in eV units
double tempI;  		  // Temperature of the ion species in eV


/*Simulation Normalizing Parameters*/
double n0, ne0, ni0;
double L, Lx;
double wp;
double CS;
double A; 

/* CHANGED TYPE FROM CONST TO VAR FOR INPUT DATA CONTROL  */
int nParticlesE;     // Number of simulation electrons
int nParticlesI;     // Number of simulation ions

int NC;              // Total number of cells
int NUM_TS;          // Total Time steps (default)
int write_interval, write_interval_phase;

string output; // Open up the output file

/* Class Domain: Hold the domain paramassEters*/
class Domain
{
public:
    int ni;      // Number of nodes
    double x0;   // initial position
    double dx;   // cell spacing
    double xl;   // domain length
    double xmax; // domain maximum position

    /* Field Data structures */
    double *phi; // Electric Potential
    double *ef;  // Electric field
    double *rho; // Charge Density
};

/* Class Particle: Hold particle position, velocity and particle identity*/
class Particle
{
public:
    double pos;  // particle position
    double vel; // particle velocity
    int id;  // hold particle identity

    // Add a constructor
    Particle(double x, double v):pos(x), vel(v){};
};

/* Class Species: Hold species data*/
class Species
{
public:
	int nc = 0;
    // Use linked list for the particles
    list<Particle> part_list;
    double mass;
    double charge;
    double spwt;
    string name;

    int NUM;
    double Temp;
    double *den;
    double *vel;

    void add(Particle part)
    {
        part.id=part_id++;
        part_list.push_back(part);
    }

    // Add a constructor
    Species(string name, double mass, double charge, double spwt, int NUM, double Temp)
    {
        setName(name);
        setMass(mass);
        setCharge(charge);
        setSpwt(spwt);
        setNum(NUM);
        setTemp(Temp);
    }

    // Define the constructor functions
    void setName(string name){this->name = name;}
    void setMass(double mass){this->mass = mass;}
    void setCharge(double charge){this->charge = charge;}
    void setSpwt(double spwt){this->spwt = spwt;}
    void setNum (int NUM){this->NUM = NUM;}
    void setTemp(double Temp){this->Temp = Temp;}

private:
    int part_id = 0;
};

// Define Domain and File as the global variable
Domain domain;

// Open files to hold numerical data
FILE *file_res; // File to hold the data of overall data
FILE *file_ke; 	// File to hold the Kinetic energy data
FILE *file_pe; // File to hold the Potential energy data
FILE *f1; 	// File to hold the ion particle data
FILE *f2; 	// File to hold the cold electron particle data
FILE *f3; 	// File to hold the hot electron particle data
FILE *f4;	// File to hold the beam electron particle data
FILE *file_loc; // File to hold data at a specific location of simulation domain

FILE *file_denflucI; 	// File to hold the density fluctuations for ions
FILE *file_denflucE;	// File to hold the density fluctuations for cold electrons

// Define Helper functions
void Init(Species *species, string flag);
void ScatterSpecies(Species *species);
void ScatterSpeciesVel(Species *species);
void ComputeRho(Species *ions, Species *electrons);
void ComputeEF(double *phi, double *ef);
void PushSpecies(Species *species, double *ef);
void RewindSpecies(Species *species, double *ef);
void SanityCheck(double new_pos, double old_pos, double cell_len);

// [Write Outputs]
void Write_ts(int ts, Species *ions,Species *electrons);
void Write_Particle(FILE *file, int ts, Species *species);
void WriteKE(double Time, Species *ions, Species *electrons);
void WriteLocation(double Time, double pos);
void WriteDenFluc(FILE *file, double Time, double pos, Species *species);
void Write_Particle(FILE *file, int ts, Species *species);

double ComputeKE(Species *species, Species *electrons);
double XtoL(double pos);
double gather(double lc, double *field);
void scatter(double lc, double value, double *field); 
void scatter_qspline(double lc, double value, double *field);

// [Sample Velocities]
double SampleVelIon(double T, double mass);
double SampleVelElectron(double T, double mass);

//[ Potential Solvers]
bool SolvePotential(double *phi, double *rho);
bool SolvePotentialDirect(double *phi, double *rho);

// [Ini Parser File]
int parse_ini_file(char * ini_name)
{
    dictionary  *   ini;

    ini = iniparser_load(ini_name);
    if (ini==NULL) {
        fprintf(stderr, "cannot parse file: %s\n", ini_name);
        return -1 ;
    }
    //iniparser_dump(ini, stderr); // Comment out to fix issues with iniparser

    /* [File] */
    output		  = iniparser_getstring(ini,"file:output",NULL);

	/* [Universal Constants] */
	EPS = iniparser_getdouble(ini,"constants:EPS",-1.0);
	K   = iniparser_getdouble(ini,"constants:K",-1.0);
	eV  = iniparser_getdouble(ini,"constants:eV",-1.0);
	e   = iniparser_getdouble(ini,"constants:e",-1.0);
	AMU = iniparser_getdouble(ini,"constants:AMU",-1.0);
	EV_TO_K  = iniparser_getdouble(ini,"constants:EV_TO_K",-1.0);
	pi  = iniparser_getdouble(ini,"constants:pi",-1.0);

    /* [Time] */
    NUM_TS     = iniparser_getint(ini,"time:NUM_TS",-1);
    DT_coeff   = iniparser_getdouble(ini,"time:DT_coeff",-1.0);
    
    /* [Diagnostics] */
    write_interval = iniparser_getint(ini,"diagnostics:write_interval",-1);
    write_interval_phase = iniparser_getint(ini,"diagnostics:write_interval_phase",-1);
    
    /* [Grid] */
    NC         = iniparser_getint(ini,"grid:NC",-1);

    /* [Population] */
    nParticlesI    = iniparser_getint(ini,"population:nParticlesI",-1);
    nParticlesE    = iniparser_getint(ini,"population:nParticlesE",-1);
    A              = iniparser_getdouble(ini,"population:A",-1);
    massI          = iniparser_getdouble(ini,"population:massI",-1.0);
    massI		   = massI*AMU;
    massE          = iniparser_getdouble(ini,"population:massE",-1.0);

    chargeE        = iniparser_getdouble(ini,"population:chargeE",-1.0);
    density        = iniparser_getdouble(ini,"population:density",-1.0);

    tempE	      = iniparser_getdouble(ini,"population:tempE",-1.0);
    tempI		  = iniparser_getdouble(ini,"population:tempI",-1.0);
         
    // ------------------------------------------------------------------------ 
    /* Normalizing Parameters */
    n0 = density;
    ni0 = n0;
    ne0 = n0;     
    wp = sqrt((ne0*chargeE*chargeE)/(massE*EPS)); // Electron Plasma Frequency
    
    DT		   = DT_coeff*(1.0/wp); // This is the unnormalized time interval, typically take 1% of wp^-1. 
    Lx = 2*pi; // System Length 
    stepSize   = Lx/NC; //LD can be used only when there is temperature, but in cold plasma oscillation there is no temperature;
      
    cout << "*************** Input Sanity Check ***************" << '\n';
    bool SFLAG = true;
    if (stepSize >= 10) {
      cout<<"ERROR: stepSize is bigger than Debye length."<<endl;
      SFLAG = false;
    }
    if (DT > 0.01/wp) {
      cout<<"ERROR: timeStep is too big. The recommended value: <"<<(0.01/wp)<<" s"<<endl;
      SFLAG = false;
    }

    if (SFLAG==true) {
      cout<<"STATUS: Input parameters are compatible."<<endl;
    }
    else {
      cout<<"ERROR: Input parameters are incompatible."<<endl;
      exit (EXIT_FAILURE);
    }

    iniparser_freedict(ini);
    return 0;
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* MAIN FUNCTION ***************************/
int main(int argc, char *argv[])
{
    /************* INIPARSER ********************/
    if(argc<2) {
      cout<<"ERROR, at least one argument expected (the input file)."<<endl;
      exit (EXIT_FAILURE);
    }
    parse_ini_file(argv[1]);
    /********************************************/    
    double Time = 0;
    /*Construct the domain parameters*/
    domain.ni = NC+1;
    domain.dx = stepSize;
    domain.x0 = 0;
    domain.xl = (domain.ni-1)*domain.dx;
    domain.xmax = domain.x0 + domain.xl;
    
    cout << "System Length: " << Lx << '\t' << "Plasma Frequency: " << wp << '\t' << endl;
    cout << "DX: " << stepSize << '\t' << "DT: " << DT << '\t' << "System Length:" << domain.xl << endl;

    /*Allocate massEmory to the domain data structures (Field variables)*/
    domain.phi = new double[domain.ni];
    domain.ef  = new double[domain.ni];
    domain.rho = new double[domain.ni];

    /*Redifine the field variables */
    double *phi = domain.phi;
    double *ef = domain.ef;
    double *rho = domain.rho;

    /* Clear the domain fields*/
    memset(phi,0,sizeof(double)*domain.ni);
    memset(ef, 0,sizeof(double)*domain.ni);
    memset(rho,0,sizeof(double)*domain.ni);

    /**************************************************/

    /*Species Info: Create vector to hold the data*/
    vector <Species> species_list;

    /*Calculate the specific weights of the ions and electrons*/
    double ion_spwt = (ni0*domain.xl)/(nParticlesI);    // Normalized with L
    double electron_spwt = (ne0*domain.xl)/(nParticlesE); // Normalized with L

    /* Add singly charged Positive ions, electrons and Background Neutrals */
    /************************************************************************/
    /* Create the species lists*/
    species_list.emplace_back("Ion", massI, chargeE, ion_spwt, nParticlesI, tempI);
    species_list.emplace_back("Electrons", massE, -chargeE, electron_spwt, nParticlesE, tempE);
    
    /*Assign the species list as ions and electrons*/
    Species &ions = species_list[0];
    Species &electrons = species_list[1]; 

    /*Initiate the species density and velocity fields*/
    ions.den = new double[domain.ni];
    electrons.den = new double[domain.ni];    

    ions.vel = new double[domain.ni];
    electrons.vel = new double[domain.ni];
   
    /*Initialize electrons and ions */
    Init(&ions,"ion");
    Init(&electrons,"electron");    

    for(auto &p:species_list)
        cout<< p.name << '\n' << "mass: " << p.mass<< '\n' << "Charge: " << p.charge << '\n' << "spwt: " << p.spwt << '\n' << "Num: " << p.NUM << endl;
    /***************************************************************************/

    /*Compute Number Density*/
    ScatterSpecies(&ions);
    ScatterSpecies(&electrons);
    
    /*Compute charge density, solve for potential
     and compute the electric field*/

    ComputeRho(&ions, &electrons);
    SolvePotential(phi, rho);
    ComputeEF(phi,ef);

    RewindSpecies(&ions,ef);
    RewindSpecies(&electrons,ef);
       
    /*------------- Print Output ---------------*/

    /*create a folder named output and
	delete the previous output folder: print statement is just to show*/
    printf("Deleting the output folder ... \n");
    system(("rm -rf "+ output).c_str());

    /*create an output folder*/
    //system("mkdir output");
    system(("mkdir "+ output).c_str());
    system(("mkdir "+ output + "/files").c_str());

    char NAME[50];

    sprintf(NAME,"%s/files/results.txt",output.c_str());
    file_res = fopen(NAME,"w");

    sprintf(NAME,"%s/files/ke.txt",output.c_str());
    file_ke = fopen(NAME,"w");
    
    sprintf(NAME,"%s/files/potloc.txt",output.c_str());
    file_loc = fopen(NAME,"w");

    sprintf(NAME,"%s/files/denFlucI.txt",output.c_str());
    file_denflucI = fopen(NAME,"w");
	
    sprintf(NAME,"%s/files/denFlucE.txt",output.c_str());
    file_denflucE = fopen(NAME,"w");


    /*MAIN LOOP*/
    clock_t start = clock();
    for (int ts=0; ts<NUM_TS+1; ts++)
    {
        //Compute number density
        ScatterSpecies(&ions);
        ScatterSpecies(&electrons);	               

        //Compute charge density
        ComputeRho(&ions, &electrons);

        //SolvePotential(phi, rho);
        SolvePotentialDirect(phi, rho);
        ComputeEF(phi, ef);

	    //move particles
        //PushSpecies(&ions, ef);
        PushSpecies(&electrons, ef);
	    
        //Write diagnostics
        if(ts%write_interval == 0)
        {            
            double max_phi = phi[0];
            for(int i=0; i<domain.ni; i++)
                if (phi[i]>max_phi) max_phi=phi[i];			

            /*print diagnostics to the screen*/
	        printf("TS: %i \t delta_phi: %.3g \t nI:%ld \t nEC:%ld\n",
				   ts, max_phi-phi[0],ions.part_list.size(),electrons.part_list.size());

	        /*Write time evolution of plasma profiles and Kinetic energy*/
    	    Write_ts(ts, &ions, &electrons);
	        WriteKE(Time, &ions, &electrons);        
	    	
        }
    if(ts%write_interval_phase == 0)
        {
            sprintf(NAME,"%s/i%d.txt",output.c_str(),ts);
            f1 = fopen(NAME,"w"); // open this file should you intend to write ion particle data

            sprintf(NAME,"%s/e%d.txt",output.c_str(),ts);
            f2 = fopen(NAME,"w"); // open this file should you intend to write electron particle data
            /*Write individual particle data to  the file*/
    	    Write_Particle(f1,ts, &ions);
    	    Write_Particle(f2,ts, &electrons);
        }

        Time += DT;
    }
    clock_t end = clock();

    /*free up memory*/
    delete phi;
    delete rho;
    delete ef;
    cout << "Time = " << ((end-start)/(double)CLOCKS_PER_SEC)/60 << " minutes" << endl;
    return 0;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/********************* HELPER FUNCTIONS ***************************/

/*Initialize the particle data : initial positions and velocities of each particle of each species*/
void Init(Species *species, string flag)
{
    // sample particle positions and velocities
    for(int p=0; p<species->NUM; p++)
    {
        double x, v;
        /*Random Loading*/ 
        // double x = domain.x0 + rnd()*(domain.ni-1)*domain.dx;
        /*Uniform Loading*/
        // double x = domain.x0 +  p * (domain.xl/species->NUM)    
        if	(flag == "electron")
        {
            //x = domain.x0 + rnd()*(domain.ni-1)*domain.dx; // Random loading
            x = domain.x0 +  p * (domain.xl/species->NUM); // uniform loading
            //----------------------------------------------------------------            
            double k = 2*pi/Lx;        
            //----------------------------------------------------------------
            double dx = A*sin(k*x);
            //----------------------------------------------------------------
            x = x + dx;
            v = 0;           
        }        
        else if (flag == "ion")
        {
            //x = domain.x0 + rnd()*(domain.ni-1)*domain.dx;
            x = domain.x0 +  p * (domain.xl/species->NUM); // uniform loading
            v = 0;             
        }
	    // Add to the particle list
        species->add(Particle(x,v));
    }
}

/*Covert the physical coordinate to the logical coordinate*/
double XtoL(double pos)
{
    double li = (pos-domain.x0)/domain.dx;
    return li;
}

/*scatter the particle data to the massEsh and collect the densities at the massEsh */
void scatter(double lc, double value, double *field)
{
    int i = (int)lc;
    double di = lc-i;
    field[i] += value*(1-di);
    field[i+1] += value*(di);
}

/* Gather field values at logical coordinates*/
double gather(double lc, double *field)
{
    int i=(int)lc;
    double di = lc-i;
    double val = field[i]*(1-di) + field[i+1]*(di);
    return val;
}

/*Scatter the particles to the massEsh for evaluating densities*/
void ScatterSpecies(Species *species)
{
    /*grab a pointer to the species density data and change
     the density field using the pointer*/

    double *field = species->den;

    /*clear the field*/
    memset(field,0,sizeof(double)*domain.ni);

    /*scatter particles to the mesh*/
    for(auto &p:species->part_list)
    {
        double lc = XtoL(p.pos);
        scatter(lc,species->spwt,field);
        //scatter_qspline(lc,species->spwt,field);
    }

    /*divide by cell volume*/
	/*Hint: we divide by the cell volume because we have scattered
	the spwt to the grid, which is just a number. Again number per
	unit volume is density. Hence we further divide the field by the cell volume*/
    for(int i=0; i<domain.ni; i++){
    	field[i] /=(domain.dx);}

    field[0] += field[domain.ni-1];
    field[domain.ni-1] = field[0];
}


//*******************************************************
void PushSpecies(Species *species, double *ef)
{
    // compute charge to mass ratio
    double qm = species->charge/species->mass;
    list<Particle>::iterator it = species->part_list.begin();

    // loop over particles
    while (it!=species->part_list.end())
    {
        // grab a reference to the pointer
        Particle &part = *it;

        // compute particle node position
        double lc = XtoL(part.pos);

        // gather electric field onto particle position
        double part_ef = gather(lc,ef);
                       
        // Grab the old positions of the particles
        double old_pos = part.pos;

        // advance velocity
        part.vel += qm*part_ef*DT;

        // Advance particle position
        part.pos += part.vel*DT;

	    // Grab the new positions of the particles
	    double new_pos = part.pos;

	    // Check whether a particle is crossing one full cell length
	    SanityCheck(new_pos, old_pos, domain.dx);

        // Take care of the particle that left the Domain (PBC)
		if (part.pos < domain.x0)
		{
			part.pos = part.pos + domain.xl;
		}
		else if(part.pos>=(domain.x0+domain.xl))
		{
			part.pos = part.pos - domain.xl;
		}
			it++;
    }
}
//*********************************************************
void SanityCheck(double new_pos, double old_pos, double cell_len)
{
	if (abs(new_pos-old_pos) > cell_len)
	{
		printf("Alert! Particles are crossing one full cell!\n");
		//exit(-1);
	}
}
// ********************************************************
/*Rewind particle velocities by -0.5*DT */
void RewindSpecies(Species *species, double *ef)
{
    // compute charge to mass ratio
    double qm = species->charge/species->mass;
    for(auto &p:species->part_list)
    {
        // compute particle node position
        double lc = XtoL(p.pos);

        // gather electric field onto the particle position        
        double part_ef = gather(lc,ef);
        
        //Rewind velocity		
        p.vel -= 0.5*qm*part_ef*DT;
    }
}

/* Compute the charge densities */
void ComputeRho(Species *ions, Species *electrons)
{
    double *rho = domain.rho;
    memset(rho,0,sizeof(double)*domain.ni);

    for(int i=0; i<domain.ni; i++)
        rho[i] = ions->charge*ions->den[i] + electrons->charge*electrons->den[i];
}

/* Potential Solver: 1. Gauss-Seidel 2. Direct-Solver*/
bool SolvePotential(double *phi, double *rho)
{
    double L2;
    double dx2 = domain.dx*domain.dx;

    // Initialize boundaries
    phi[0]=phi[domain.ni-1]=0;

    // Main Solver
    for(int it=0; it<2000000; it++)
    {
        for(int i=1; i<domain.ni-1; i++)
        {
            double g = 0.5*(phi[i-1] + phi[i+1] + dx2*rho[i]/EPS);
            phi[i]=phi[i] + 1.4*(g-phi[i]);
        }
        // Check for convergence
        if(it%25==0)
        {
            double sum = 0;
            for(int i=1; i<domain.ni-1; i++)
            {
                double R = - rho[i]/EPS - (phi[i-1]-2*phi[i]+phi[i+1])/dx2;
                sum += R*R;
            }
            L2 = sqrt(sum)/domain.ni;
            if(L2<1e-4){return true;}

        }
        //printf("GS-Converged! L2=%g\n",L2);
    }
    printf("Gauss-Siedel solver failed to converge, L2=%g\n",L2);
    return false;
}

/* Potential Direct Solver */

bool SolvePotentialDirect(double *x, double *rho)
{
    /* Set coefficients, precompute them*/
    int ni = domain.ni;
    double dx2 = domain.dx*domain.dx;
    double *a = new double[ni];
    double *b = new double[ni];
    double *c = new double[ni];

    /*Centtral difference on internal nodes*/
    for(int i=1; i<ni-1; i++)
    {
        a[i] = 1; b[i] = -2; c[i] = 1;
    }

    /*Apply dirichlet boundary conditions on boundaries*/
    a[0]=0; b[0]=1; c[0]=0;
    a[ni-1]=0; b[ni-1]=1; c[ni-1]=0;

    /*multiply R.H.S.*/
    for (int i=1; i<ni-1; i++)
        x[i]=-rho[i]*dx2/EPS;

    x[0] = 0;
    x[ni-1] = 0;

    /*Modify the coefficients*/
    c[0] /=b[0];
    x[0] /=b[0];

    for(int i=1; i<ni; i++)
    {
        double id = (b[i]-c[i-1]*a[i]);
        c[i] /= id;
        x[i] = (x[i]-x[i-1]*a[i])/id;
    }

    /* Now back substitute */
    for(int i=ni-2; i>=0; i--)
        x[i] = x[i] - c[i]*x[i+1];

    return true;
}

/*Compute electric field (differentiating potential)*/
void ComputeEF(double *phi, double *ef)
{
    /*Apply central difference to the inner nodes*/
    for(int i=1; i<domain.ni-1; i++)
        ef[i] = -(phi[i+1]-phi[i-1])/(2*domain.dx);

    // Apply central difference to the boundary nodes in periodic case 
    ef[0] = -(phi[1]-phi[domain.ni-2])/(2*domain.dx);
    ef[domain.ni-1] = -(phi[1]-phi[domain.ni-2])/(2*domain.dx);

    /*Apply one sided difference at the boundary nodes*/
    // ef[0] = -(phi[1]-phi[0])/domain.dx;
    // ef[domain.ni-1] = -(phi[domain.ni-1]-phi[domain.ni-2])/domain.dx;
}


/*Write the output with Time*/
void Write_ts(int ts, Species *ions, Species *electrons)
{
    for(int i=0; i<domain.ni; i++)
    {
        fprintf(file_res,"%g \t %g \t %g \t %g \t %g \t %g\n", i*domain.dx, ions->den[i], electrons->den[i],
		domain.rho[i], domain.phi[i], domain.ef[i]);

    }
    fflush(file_res);
}

void WriteKE(double Time, Species *ions, Species *electrons)
{
    double ke_ions = ComputeKE(ions, electrons);
    double ke_electrons = ComputeKE(electrons, electrons);    
    
    fprintf(file_ke,"%g \t %g \t %g\n",Time, ke_ions, ke_electrons);

    fflush(file_ke);
}

/* Write the output particle data to respective files */
void Write_Particle(FILE *file, int ts, Species *species)
{
    for(auto& p: species->part_list)
    {
        fprintf(file,"%g \t %g\n",p.pos, p.vel);
    }
    fflush(file_res);
}

double ComputeKE(Species *species, Species *electrons_cold)
{
    double ke = 0;
    for (auto &p:species->part_list)
    {
        // un-normalize the velocity by multiplying with the cold thermal velocity
        ke += (p.vel*p.vel);
    }
    /*Multiply 0.5*mass for all particles*/
    ke *= 0.5*(species->spwt*species->mass);
    return ke;
}
