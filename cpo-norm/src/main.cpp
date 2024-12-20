/*
 ****************************************************
 CODE: epic (electrostatic particle-in-cell)
 
 DEVELOPER:
 ----------------------------------------
 DR. RAKESH MOULICK, CPP-IPR, ASSAM, INDIA
 
 ACKNOWLEDGEMENTS:
 ----------------------------------------
 DR. SAYAN ADHIKARI, IFE, NORWAY
 DR. LUBOS BRIEDA, PIC LLC, USA
 MR. KAUSHIK KALITA, CPP-IPR

 OBJ: TO OBSERVE COLD PLASMA OSCILLATION
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
double density;       // plasma Density
double stepSize;      // cell Spacing
double DT;            // time step
double DT_coeff;      // percentage of wp^-1 to be used e.g. for 1%(wp^-1) it is 0.01

double tempE;         // temperature of electrons in eV units
double tempI;  		  // temperature of ions  in eV units


/*Define individual species density and other Parameters*/
double n0, ne0, ni0;
double wp, wpe, wpi; 

double A;  // perturbation amplitude 
double Lx; // domain length

/*Normalization parameters*/
double Delta, L, W, phi0, vT, posCoff, velCoff, rhoCoff, ECoff;

/* CHANGED TYPE FROM CONST TO VAR FOR INPUT DATA CONTROL  */
int nParticlesE;     // number of simulation electrons
int nParticlesI;     // number of simulation ions

int NC;              // total number of cells
int NUM_TS;          // total Time steps (default)
int write_interval, write_interval_phase;  // data writing interval 

string output;       // output file

/* Class Domain: Hold the domain paramassEters*/
class Domain
{
public:
    int ni;      // number of domain nodes
    double x0;   // initial position
    double dx;   // cell spacing
    double xl;   // domain length
    double xmax; // domain maximum position

    /* Field Data structures */
    double *phi; // electric Potential
    double *ef;  // electric field
    double *rho; // charge Density
};

/* Class Particle: Hold particle position, velocity and particle identity*/
class Particle
{
public:
    double pos;  // particle position
    double vel;  // particle velocity
    int id;      // particle identity

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
FILE *file_pe;  // File to hold the Potential energy data
FILE *f1; 	    // File to hold the ion particle data
FILE *f2; 	    // File to hold the electron particle data
FILE *file_loc;      // File to hold data at a specific location of simulation domain
FILE *file_denflucI; // File to hold the density fluctuations for ions
FILE *file_denflucE; // File to hold the density fluctuations for cold electrons

// Define helper functions
void Init(Species *species, string flag);
void ScatterSpecies(Species *species);
void ScatterSpeciesVel(Species *species);
void ComputeRho(Species *ions, Species *electrons);
void ComputeEF(double *phi, double *ef);
double ComputeKE(Species *species, Species *electrons);
void PushSpecies(Species *species, double *ef);
void RewindSpecies(Species *species, double *ef);
void SanityCheck(double new_pos, double old_pos, double cell_len);

// [Write output functions]
void Write_ts(int ts, Species *ions,Species *electrons);
void Write_Particle(FILE *file, int ts, Species *species);
void WriteKE(double Time, Species *ions, Species *electrons);
void WriteLocation(double Time, double pos);
void WriteDenFluc(FILE *file, double Time, double pos, Species *species);

// [Compute functions]
double XtoL(double pos);
double gather(double lc, double *field);
void scatter(double lc, double value, double *field); 
void scatter_qspline(double lc, double value, double *field);

// [Sample Velocities]
double SampleVelIon(double T, double mass);
double SampleVelElectron(double T, double mass);

// [Potential Solvers]
bool SolvePotential(double *phi, double *rho);
bool SolvePotentialDirect(double *phi, double *rho);

// [Ini Parser File: Call the values from the input.ini file]
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

    /* [DIAGNOSTICS] */
    write_interval = iniparser_getint(ini,"diagnostics:write_interval",-1);
    write_interval_phase = iniparser_getint(ini,"diagnostics:write_interval_phase",-1);
    /* [Grid] */
    NC            = iniparser_getint(ini,"grid:NC",-1);

    /* [Population] */
    nParticlesI    = iniparser_getint(ini,"population:nParticlesI",-1);
    nParticlesE    = iniparser_getint(ini,"population:nParticlesE",-1);
    A              = iniparser_getdouble(ini,"population:A",-1.0);
    massI          = iniparser_getdouble(ini,"population:massI",-1.0);
    massI		   = massI*AMU;
    massE          = iniparser_getdouble(ini,"population:massE",-1.0);
    chargeE        = iniparser_getdouble(ini,"population:chargeE",-1.0);
    density        = iniparser_getdouble(ini,"population:density",-1.0);
    tempE	       = iniparser_getdouble(ini,"population:tempE",-1.0);
    tempI		   = iniparser_getdouble(ini,"population:tempI",-1.0);    
    
    /* Normalizing Parameters */
    n0 = density; // equilibrium plasma density
    ni0 = n0;     // equilibrium ion density
    ne0 = n0;     // equilibrium electron density

    wpe = sqrt((ne0*chargeE*chargeE)/(massE*EPS)); // Electron Plasma Frequency
    wpi = sqrt((ni0*chargeE*chargeE)/(massI*EPS)); // Ion Plasma Frequency
    wp = wpe + wpi;     // Total Plasma Frequency

    Lx = 2*pi;          // domain length (1D)
    stepSize   = Lx/NC; // individual cell sizes

    L = stepSize;       // spatial normalization parameter (~ stepSize)
    W = wpe;            // temporal normalization parameter (~inverse of W)
    
    /* Time step */
    DT = DT_coeff*(1.0/W); // unnormalized time interval
    DT = W*DT; // normalized time interval: temporal normalization is accomplished here

    /* Normalizing Coefficients */
    vT = L*W;                        // velocity normalization by ion thermal velocity
    phi0 = 1;
    posCoff = vT/(L*W);              // normalization coefficient for position 
    velCoff = (n0*e*L)/(EPS*W*vT);   // normalization coefficient for velocity
    rhoCoff = (e*n0*L*L)/(phi0*EPS); // normalization coefficient for charge density
    ECoff = (EPS*phi0)/(n0*e*L*L);   // normalization coefficient for electric field
   	
    cout << "*************** Input Sanity Check ***************" << '\n';
    bool SFLAG = true;
    if (stepSize > Lx/NC) {
      cout<<"ERROR: stepSize is bigger than cell sizes "<<endl;
      SFLAG = false;
    }
    if (DT > DT_coeff) {
      cout<<"ERROR: timeStep is too big. The recommended value is 1 or 2 percent of wpe^-1 [s]"<<endl;
      SFLAG = false;
    }

    if (SFLAG==true) {
      cout<<"STATUS: Input parameters are compatible "<<endl;
    }
    else {
      cout<<"ERROR: Input parameters are incompatible "<<endl;
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
    /* Print the basic normalizing parameters */
    cout << "Normalizing Length: " << L << "[m]" << endl;
    cout << "Normalizing Frequency: " << W << "[Hz]" << endl;
    cout << "DX: " << stepSize << "[m]" << endl;
    cout << "DT: " << DT << "[normalized]" << endl;

    /* Begin Time */
    double Time = 0;

    /*Construct the domain parameters*/
    domain.ni = NC+1; // number of grid points 
    domain.dx = stepSize/L; // normalized cell size (spatial normalization is accomplished here)
    domain.x0 = 0; // initial value of x (in normalized scale) 
    domain.xl = (domain.ni-1)*domain.dx; // length of the domain (in normalized scale) 
    domain.xmax = domain.x0 + domain.xl; // maximum domain length (in normalized scale)
    
    /*Allocate massEmory to the domain data structures (Field variables)*/
    domain.phi = new double[domain.ni]; // electric potential vector
    domain.ef  = new double[domain.ni]; // electric field vector
    domain.rho = new double[domain.ni]; // charge density vector

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

    /*Calculate the specific weights of ions and electrons*/
    // Since, the domain length is already normalized, we need to un-normalize it for the calculation
    double ion_spwt = (ni0*domain.xl*L)/(nParticlesI);    // particle weight factor of ions
    double electron_spwt = (ne0*domain.xl*L)/(nParticlesE); // particle weight factor of electrons
    
    /* Add singly charged Positive ions, electrons and Background Neutrals */
    /************************************************************************/
    /* Create the species lists*/
    species_list.emplace_back("Ion", massI,chargeE,ion_spwt, nParticlesI, tempI);
    species_list.emplace_back("Electrons", massE,-chargeE,electron_spwt, nParticlesE, tempE);
    
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
    
    /*Compute charge density, solve for electric potential
     and compute the electric field*/

    ComputeRho(&ions, &electrons);
    SolvePotential(phi, rho); 
    ComputeEF(phi,ef);

    /*Rewind the velcities to account for leap-froging*/
    RewindSpecies(&ions,ef);
    RewindSpecies(&electrons,ef);
       
    /*------------- Print Output ---------------*/

    /* Delete any existing output folder */
    printf("Deleting the output folder ... \n");
    system(("rm -rf "+ output).c_str());

    /*create a new output folder*/
    system(("mkdir "+ output).c_str());

    /*create another folder named 'files' within the output folder */
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

        //Compute velocity field
        //ScatterSpeciesVel(&ions);
        //ScatterSpeciesVel(&electrons);    

        //Compute charge density
        ComputeRho(&ions, &electrons);

        //SolvePotential(phi, rho);
        SolvePotentialDirect(phi, rho);
        ComputeEF(phi, ef);

	    //move particles: move only electrons for cold plasma oscillation
        // PushSpecies(&ions, ef);
        PushSpecies(&electrons, ef);
	    
        //Write diagnostics
        if(ts%write_interval == 0)
        {      
            double max_phi = phi[0];
            for(int i=0; i<domain.ni; i++)
                if (phi[i]>max_phi) max_phi=phi[i];

            /*print diagnostics to the screen*/
	        printf("TS: %i \t delta_phi: %.3g \t nI:%ld \t nEC:%ld\n",
				   ts, max_phi-phi[0], ions.part_list.size(),electrons.part_list.size());

	        /*Write time evolution of plasma profiles and Kinetic energy*/
    	    Write_ts(ts, &ions, &electrons);
	        WriteKE(Time, &ions, &electrons);
            
	        /*Write Electric Field Data at the Mid Plane*/
	        WriteLocation(Time, domain.xl/2);

	        /*Write the Density Fluctuation at the Mid Plane*/
	        WriteDenFluc(file_denflucI, Time, domain.xl/2, &ions);
	        WriteDenFluc(file_denflucE, Time, domain.xl/2, &electrons);
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
            double k = 2*pi/Lx; // mode to be excited            
            //----------------------------------------------------------------
            double dx = (1/L)*A*sin(k*x*L);
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

/*Sample Velocity (According to Birdsall)*/
double SampleVelIon(double T, double mass)
{
    double v_th = sqrt(2*K*T/mass);
	double vt = v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5);
    // Normalize particle velocity by the thermal velocity of cold electrons & return
	return vt/(vT); 
}

double SampleVelElectron(double T, double mass)
{
	double v_th = sqrt(2*K*T/mass);
    double vt = v_th*sqrt(2)*(rnd()+rnd()+rnd()-1.5);
    // Normalize particle velocity by the thermal velocity of cold electrons & return
    return vt/(vT);
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

void scatter_qspline(double lc, double value, double *field)
{
    /*compute nearest grid point*/    
    // int i = (int) (lc + 0.5);
    int i = round(lc);
    double di = lc-i;
                
    /*deposit value to center node if the left/right neighbor is out of domain*/
    double v_L = value*0.5*(0.5-di)*(0.5-di);
    double v_C = value*(0.75-di*di);
    double v_R = value*0.5*(0.5+di)*(0.5+di);

    field[i] += v_C;
    
    if (i>0) field[i-1] += v_L;
    if (i<domain.ni-1) field[i+1] += v_R;

    // if (i>0) field[i-1] += v_L; else field[i] +=v_L;
    // if (i<domain.ni-1) field[i+1] += v_R; else field[i] +=v_R;

}

double gather_qspline(double lc, double *field)
{
    int i = (int)lc;
    double di = lc-i;

    double v_C = (0.75-di*di);
    double v_L = 0.5*(0.5-di)*(0.5-di);
    double v_R = 0.5*(0.5+di)*(0.5+di);	

    double val = field[i-1]*v_L + field[i]*v_C + field[i+1]*v_R; 
    return val;
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
        // scatter_qspline(lc,species->spwt,field);
    }

    /*divide by cell volume*/
	/*Hint: we divide by the cell volume because we have scattered
	the spwt to the grid, which is just a number. Again number per
	unit volume is density. Hence we further divide the field by the cell volume*/
    for(int i=0; i<domain.ni; i++){
    	field[i] /=(domain.dx*L);}

	/*Normalize the field value*/
	for(int i=0; i<domain.ni; i++){
		field[i] /=density;}

    // Adjust the end nodes for periodic boundary condition
    field[0] += field[domain.ni-1];
    field[domain.ni-1] = field[0];
    
    // For Qudratic Spline Interpolation
    // field[1] = field[domain.ni-3];
    // field[domain.ni-2] = field[2];
}

/*Scatter the particles to the massEsh for evaluating velocities*/
void ScatterSpeciesVel(Species *species)
{
    /*grab a pointer to the species velocity field and change
     the velocity field using the pointer*/
    double *field = species->vel;

    /*clear the field*/
    memset(field,0,sizeof(double)*domain.ni);

    /*scatter particles to the mesh*/
    for(auto &p:species->part_list)
    {
        double lc = XtoL(p.pos);
        scatter(lc,species->spwt*p.vel,field);
    }

    /*divide by cell volume*/
    for(int i=0; i<domain.ni; i++){
        field[i] /=(species->den[i]*density*domain.dx*L);}

    /*adjusting values at the boundaries*/ 
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
        // double part_ef = gather_qspline(lc,ef);
        
        // Grab the old positions of the particles
        //double old_pos = part.pos;

        part.vel += qm*velCoff*part_ef*DT;
        
        // Advance particle position
        part.pos += posCoff*part.vel*DT;
	    
        // Grab the new positions of the particles
	    //double new_pos = part.pos;

	    // Check whether a particle is crossing one full cell length
	    //SanityCheck(new_pos, old_pos, domain.dx);

        // Take care of the particle that left the Domain (PBC)
		if (part.pos < domain.x0)
		{
			part.pos = part.pos + domain.xl;
		}
		else if(part.pos>=(domain.x0 + domain.xl))
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
        // double part_ef = gather_qspline(lc,ef);
        
        //advance velocity
        p.vel -= 0.5*qm*velCoff*part_ef*DT;
    }
}

/* Compute the charge densities */
void ComputeRho(Species *ions, Species *electrons)
{
    double *rho = domain.rho;
    memset(rho,0,sizeof(double)*domain.ni);

    for(int i=0; i<domain.ni; i++)
        rho[i] = rhoCoff*(ions->den[i] - electrons->den[i]);
}

/* Potential Solvers: 1. Gauss-Seidel 2. Direct-Solver*/
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
            double g = 0.5*(phi[i-1] + phi[i+1] + dx2*rho[i]);
            phi[i]=phi[i] + 1.4*(g-phi[i]);
        }
        // Check for convergence
        if(it%25==0)
        {
            double sum = 0;
            for(int i=1; i<domain.ni-1; i++)
            {
                double R = - rho[i] - (phi[i-1]-2*phi[i]+phi[i+1])/dx2;
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
        x[i]=-rho[i]*dx2;

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
        ef[i] = -ECoff*((phi[i+1]-phi[i-1])/(2*domain.dx));
    
    // Apply central difference to the boundary nodes in periodic case 
    ef[0] = -ECoff*(phi[1]-phi[domain.ni-2])/(2*domain.dx);
    ef[domain.ni-1] = -ECoff*(phi[1]-phi[domain.ni-2])/(2*domain.dx);
}


/*Write the output data with Time*/
void Write_ts(int ts, Species *ions, Species *electrons)
{
    for(int i=0; i<domain.ni; i++)
    {
        fprintf(file_res,"%g \t %g \t %g \t %g \t %g \t %g\n", i*domain.dx, ions->den[i], electrons->den[i],
		domain.rho[i], domain.phi[i], domain.ef[i]);

    }
    fflush(file_res);
}

/* Write the electric field data at a particular location */
void WriteLocation(double Time, double pos)
{
	double lc = XtoL(pos);
	int i = (int) lc;
	fprintf(file_loc,"%g \t %g\n",Time, domain.ef[i]);
}

/* Write density fluctuation data*/
void WriteDenFluc(FILE *file, double Time, double pos, Species *species)
{
    double normalized_equilibrium_density = 1.0;
	double lc = XtoL(pos);
	int i = (int) lc;
	fprintf(file,"%g \t %g\n",Time, (species->den[i] - normalized_equilibrium_density));
	fflush(file);
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

/* Write kinetic energy data*/
void WriteKE(double Time, Species *ions, Species *electrons)
{
    double ke_ions = ComputeKE(ions, electrons);
    double ke_electrons = ComputeKE(electrons, electrons);    
    
    fprintf(file_ke,"%g \t %g \t %g\n",Time, ke_ions, ke_electrons);

    fflush(file_ke);
}

/* Compute the kinetic energies*/
double ComputeKE(Species *species, Species *electrons_cold)
{
    double ke = 0;
    for (auto &p:species->part_list)
    {
        // un-normalize the velocity by multiplying with the cold thermal velocity
        ke += (p.vel*p.vel)*(vT)*(vT);
    }
    /*Multiply 0.5*mass for all particles*/
    ke *= 0.5*(species->spwt*species->mass);

    return ke;
}
