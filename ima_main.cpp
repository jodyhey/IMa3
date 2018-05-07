/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */

#define GLOBVARS
#include "ima.hpp"
#include <stack>

/*  some compilation commands
mpicxx *.cpp -D MPI_ENABLED -D NDEBUG -O3 -o IMa3
mpicxx *.cpp -D MPI_ENABLED -D NDEBUG -O3 -o IMa3 -fpermissive
 */

/********* variables local to ima_main.cpp.***********/   /*   global variables are declared in ima.hpp using gextern/extern */

static char *infilename;
static char *loadfilebase;
static char *outfilename;
static char *burnfilename;
static char *runfilename;
static char command_line[COMMANDLINESTRINGLENGTHMAX];
static char defaultpriorfilename[14]= "imapriors.txt";
static char fpstr[FPSTRIMAXLENGTH]; // probably long enough to hold entire output file text
static char genealogyinfosavefilename[FNSIZE];
static char heatingterm_str[50], modeloptions_str[50], calcoptions_str[50], outputoptions_str[50],runoptions_str[50], priors_str[50];
static char infile_name[FNSIZE];
static char mcfreadfilename[FNSIZE];
static char mcfwritefilename[FNSIZE];
static char migplotfilename[FNSIZE];
static char migrationnamefilename[FNSIZE];
static char nestedmodelfilename[FNSIZE] = {'\0'};
static char oldoutfilename[FNSIZE];
static char priorfilename[FNSIZE];
static char timeinfostring[80];
static char topologypriorinfostring[COMMANDLINESTRINGLENGTHMAX] = {'\0'};
static double generationtime;
static double hilike, hiprob;
static double currlike, currprob;
//static double hilike_rec,hiprob_rec;
static double hval1, hval2;
static double scaleumeaninput = 0;
static FILE *genealogyinfosavefile;
static FILE *migplotfile;
static FILE *migrationnamefile;
static FILE *outfile;
static int **migcount,*migcountdata, *migfrom, *migto;
static int *fpstri;
static int burndone;
static int burndurationmode, cdurationmode;
static int burntrendstartdelay;
static int continuerun = 1;
static int usegenerationtimedefault;
static int genealogiestosave;
static int gsampinflength;
static int heatmode;
static int lastgenealogysaved = 0;
static int maxedoutgenealogysave = 0;
static int memforgenealogiessaved = 0;
static int migrationnamefrom,migrationnameto;
static int numgenealogyfiles;
static int phylogeniestorecord;
static int printint;
static int recordint;
static int savegenealogyint;
static int swaptries_per_cpu;
static int time_t_size= sizeof(time_t);
static int trenddoublepoint;
static int trendspot = 0;
static long burnduration, chainduration;
static int mcmcrecords = 0;
static int mcmcrecords_old = 0;
static int genealogysamples_old = 0;
static int burninsteps_old = 0;
static int noburn_mcfload = 0;
static int runsteps_old = 0;
static long seed_for_ran1;
static time_t endtime;
static time_t lasttime;
static time_t starttime;

/******* variables that are extern in some other files *********/

/* arrays that are indexed by numpops,  all should have length MAXPOPS_PHYLOGENYESTIMATION + 1 */
int numdistinctpopulationpairs[] = {0,0,1,6,25,90,301,966,3025,9330}; /* total number of possible distinct pairs of populations that could engage in gene flow (don't share any descendant pops) */
int numtreesarray[] = {0,1,1,3,18,180,2700,56700,1587600}; // numer of distinct tree topologies with different #'s of populations
int hashvalmaxes[] = {0,0,2,10,50,100,300,1000,3000,10000}; // hashvalue range

#ifdef MPI_ENABLED
MPI_Datatype MPI_updatescalar;  // used for passing updatescalar structures among processes
MPI_Op myOp_updatescalarsum;   // used to make a function for summing MPI_updatescalar using MPI_Reduce
#endif

/***** extern function prototypes *****/
extern void addoutgroup(char s[]);  // declared in alltreestings.cpp

/***** Local function prototypes  *****/
static int checkrunfile(char *runfilename);
static void free_ima_main_stuff ();
static char *releaseinfostring();
static void scan_commandline (int argc, char *argv[], int currentid);
static void begin_outputfile_info_string (void);
static void set_mode_check_inconsistent_options(void);
static void start (int argc, char *argv[], int currentid);
static void qupdate (int currentid );
static void savegenealogyinfo (int currentid);
static void reset_after_burn (int currentid);
static void writeburntrendfile (int currentid);
static int run (int currentid);
static void inctrend (int m, int t, struct value_record *v, double newval);
static void trend_reset (struct value_record *v, int nv);
static void trendrecord (int loadarrayj, int currentid);
static void recordval (struct value_record *v, double val);
static void record_migrations (int z);
static void checkhighs (int z, int currentid, int reset);
static void record (int currentid);
static void loadgenealogyvalues (void);
static void printoutput (int currentid,int finaloutput);
static void intervaloutput (FILE * outto, int currentid);
static void output_update_scalars(int z,int currentid, char updatestr[]);
static void check_to_record (int currentid);
static void record_migration_names();
static int comparethermos(const void *a, const void *b);
void commit_mpi_updatescalar(void);
#ifdef MPI_ENABLED
void jh_mpi_sum_scalerstruct(struct updatescalarinfo *in,struct updatescalarinfo *inout, int *len, MPI_Datatype *type);
#endif
//void printupdatescalarinfo (FILE * outto, int currentid);  stopped using 6/7/2017
#ifdef XMLOUTPUT
std::stack<TiXmlElement*> xstack;
#endif

int main ( int argc, char *argv[]);

int checkrunfile(char *runfilename)
{
  FILE* runfile;
  char ch;
  runfile = fopen (runfilename, "r");
  if (runfile == NULL)
  {
    return 0;
  }
  else
  {
    ch = (char) getc (runfile);
	   fflush(runfile);
    FCLOSE (runfile);
    if ((char) toupper (ch) == 'Y')
    {
		    return 1;
    }
    else
      return 0;
  }
}  // checkrunfile

void
free_ima_main_stuff ()
{
  int i;

  if (memforgenealogiessaved > 0)
  {
    for (i = 0; i < memforgenealogiessaved; i++)
    {
      XFREE (gsampinf[i]);
    }
    XFREE (gsampinf);


    if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)  // debugging genealogies when sampling trees
    {
      for (i = 0; i < memforgenealogiessaved; i++)
      {
        XFREE (debug_ti_addinfo[i]);
      }
      XFREE(debug_ti_addinfo);
    }
  }

  XFREE (fpstri);
  if (modeloptions[EXPOMIGRATIONPRIOR] || (runoptions[LOADRUN] && calcoptions[FINDJOINTPOSTERIOR]))
     XFREE(eexpsum);

  if (outputoptions[MIGRATEHIST])
  {
    free_2Dint_contiguous(migcount,migcountdata);
    XFREE (migfrom);
    XFREE (migto);
  }


  if (infilename != NULL)
  {
    XFREE (infilename);
  }
  if (outfilename != NULL)
  {
    XFREE (outfilename);
  }

  if (loadfilebase != NULL)
  {
    XFREE (loadfilebase);
  }
  freeanymemory ();
}        //free_ima_main_stuff


static char *releaseinfostring()
{
  static char rs[100];
#ifdef IMA3RELEASE
    sprintf(rs,"IMa3 Release Version: %s\n",IMA3RELEASEVERSION);
#endif
#ifdef INDEVELOPMENT
    sprintf (rs,"IMa3 in Development. Most Recent Release Version: %s",IMA3RELEASEVERSION);
#endif
  rs[sizeof(rs)-1] = '\0';  // just in case it has written past the end
  return rs;
}

/* each process scans the command line */
void
scan_commandline (int argc, char *argv[], int currentid)
{
  static int i, j, k;
  static char ch, ch1;
  static char pstr[256];
  const char *buildtimestring = "Executable compiled on " __DATE__ ", " __TIME__ ".";
  /* option flags that, although values are assigned to the flags, the flag
   * is never used by the code:  Cp, Ep, Wp, Yp.  These variable flags
   * are being left in for completeness
   */
   //JK: Adding Np var for unique IMburn/run files
  int Bp, Cp, Dp, Fp, Gp, Hfp, Hnp,/* Hkp, */ Hap, Hbp, Ip, Jp, Lp, Mp, Np, Op, Pp,
    Qp, Rp, Sp, Tp, Up, Vp, Wp, Yp,Xp, Zp;
  double tempf;

  int samplemcmcinterval;
  struct tm *starttimeinfo;
  int numthingstosave;
  int Ap = 0;  // use -a for debugging things that need a flag
  /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */

  int rc = 0; //AS: MPI C-binding error code handler
  time (&starttime);
  //time (&remained_starttime);
#ifdef MPI_ENABLED
		if (numprocesses > 1)
  {
    rc = MPI_Bcast(&starttime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
//printf("1 %d %llu\n",currentid,starttime);
	 	 if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
    //rc = MPI_Bcast(&remained_starttime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
//printf("2 %d %llu\n",currentid,remained_starttime);
	 	 //if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	 }
#endif


  starttimeinfo = localtime(&starttime);
  strftime(timeinfostring,80,"%x - %X", starttimeinfo);

  numchainstotal  = DEFAULTNUMCHAINS;      /* default value */
  numchainspp =  DEFAULTNUMCHAINS;      /* default value */
  genealogiestosave = -1;
  Bp = 0;                       /* duration of burnin */
  Cp = 0;                       /* calculation options  - flag not used */
  Dp = 0;                       /* number of steps in between genealogy saves */
  Fp = 0;                       /* name of mcf file */
  Gp = 0;                       /* name of prior file */
  Hfp = 0;                      /* heating model */
  Hnp = 0;                      /* # of chains */
  //Hkp = 0;                      /* # of swap attempts */  // jh no longer used as of 6/16/2016  now using setswaptries()
  Hap = 0;                      /* heat term1 */
  Hbp = 0;                      /* heat term2 */
  Ip = 0;                       /*input file */
  Jp = 0;                       /* used for programmer options */
  Lp = 0;                       /* duration of chain */
  Mp = 0;                       /* migration rate max */
  Np = 0;                       /* Expanded filepaths for IMburn/run */
  Op = 0;                       /* output file */
  Pp = 0;                       /* output options */
  Qp = 0;                       /* Theta max scalar */
  Rp = 0;                       /* run options */
  Sp = 0;                       /* random number seed */
  Tp = 0;                       /* Time maximum */
  Up = 0;                       /* generation time in years */
  Vp = 0;                       /* genealogy load file name base */
  Wp = 0;                       /* name of file with nested models - */
  Xp = 0;                       /* set prior probability for population topology for pairs of sister populations */
  Yp = 0;                       /* mutation rate scalar for loci with mutation rates given in input file - for use with LOADRUN mode  - flag not used */
  Zp = 0;                       /* screen printout frequency */
  if (currentid == HEADNODE) {
	  printf ("executing program");
    if (numprocesses > 1)
    {
      printf ("...This is the head node. Screen output appears here. \n");
    }
    else
      printf("\n");
    fflush(stdout);
  }
  if ((argc == 2 && ((char) toupper (argv[1][0]) == 'H' || (char) toupper (argv[1][1]) == 'H') ) || argc == 1)
  {
    if (currentid == HEADNODE)  // print help screen then exit
    {
      printf ("IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman \n");
  #ifdef STDTEST
      printf ("\n\n==================\nSTANDARD TEST MODE\n");
  #ifdef RANDOM_NUMBERS_FROM_FILE
      printf ("USES randNums FILE\n");
  #endif
      printf ("==================\n\n");
  #endif
      printf("%s\n",releaseinfostring());
      printf ("%s\n",buildtimestring);
  #ifdef DEBUG
      printf ("  -Executable compiled with DEBUG defined\n");
  #endif
  #ifdef MPI_ENABLED
      printf ("  -Executable compiled with MPI enabled\n");
  #endif
      printf ("This program is run using a command line interface\n");
      printf ("To execute, type the program name followed by command line flags and options \n");
      printf ("If compiled with MPI_ENABLED,  IMa3 can be run with multiple cpus\n");
      printf ("-b  Duration of burn  (MCMC mode only)\n");
      printf ("    - If integer, the number of burnin steps \n");
      printf ("    - If floating point, the time in hours between writing of burntrend file\n");
      printf ("         run continues until file " "IMburn" " is no longer present\n");
      printf ("         in the directory, or if present, does not begin with 'y'\n");
      printf ("-c  Calculation options: \n");
      printf ("    0 Likelihood of data functions return a constant - posterior should equal prior \n");
      printf ("    1 Include ranges on mutation rates as priors on mutation rate scalars\n");
      printf ("    2 Joint posterior density calculations, for LLR tests of nested models, use with -r0 -w\n");
      printf ("    3 Get prior distribution terms from file (requires filename given with -g )\n");
      /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
      /* as of 8/23/2011 this is still in progress,  do not show this help output in release versions
      printf ("    4 Calculate the marginal likelihood, must specify -hn (odd number > 50 chains),  -hf, -ha, and -hb are ignored\n"); */
      printf ("-d  Number of steps between sampling genealogies (default 100)\n");
      printf ("-f  Name of file with saved Markov chain state generated in previous run (use with -r3)\n");
      printf ("-g  Name of file for parameter priors, default: '%s'\n",defaultpriorfilename);
      printf ("-h  Heating terms (MCMC mode only): \n");
      //printf ("  -hf Heating model: l linear ; g geometric\n");  // 6/15/2017 stopped showing this in menu, geometric is default,
      printf ("  -hn Number of chains \n");
      //printf ("  -hk Number of chain swap attempts per step (default = number of chains)\n"); //  jh no longer used as of 6/16/2016  now using setswaptries()
      printf ("  -ha First heating parameter (less than, but near 1.  The more chains, the closer to 1)\n");
      printf ("  -hb Second heating parameter  (the smallest Beta value) \n");
      printf ("-i  Input file name (no spaces) \n");
      printf ("-j  Model options: \n");
      printf ("    0  Population Topology Updating and Estimation (for 3 or more populations)\n");
      printf ("    1  Add a non-sampled ghost population to the model \n");
      printf ("    2  Migration prior follows exponential distribution with mean given by -m or priorfile \n");
      printf ("    3  Use hyperpriors. -q and -m specify hyperprior distributions (if not -j0, priorfile is created)\n");
      printf ("    4  Migration only between sister populations (do not use with -j0)\n");
      printf ("    5  One migration parameter for each pair of populations (do not use with -j0)\n");
      printf ("    6  Migration only between sampled populations (do not use with -j0)\n");
      printf ("    7  Separate population size and migration parameters in each period (do not use with -j0) \n");
      printf ("    8  No migration in the model (do not use with -j0)\n");
      printf ("    9  One single migration parameter for all pairs of populations (do not use with -j0)\n");
      printf ("    x  One single population size parameter for all populations (do not use with -j0)\n");

      /* old ordering,  updated this on 10/30/2017
      printf ("    1  Migration only between sister populations (do not use with -j0)\n");
      printf ("    2  One migration parameter for each pair of populations (do not use with -j0)\n");
      printf ("    3  Migration only between sampled populations (do not use with -j0)\n");
      printf ("    4  Add a non-sampled ghost population to the model \n");
      printf ("    5  Separate population size and migration parameters in each period (do not use with -j0) \n");
      printf ("    6  No migration in the model (do not use with -j0)\n");
      printf ("    7  Migration prior follows exponential distribution with mean given by -m or priorfile \n");
      printf ("    8  One single migration parameter for all pairs of populations (do not use with -j0)\n");
      printf ("    9  Use hyperpriors.  Treat -q and -m values as specifying hyperprior distributions\n"); */


      printf ("-L  Run duration. Duration depends on integer or floating point value and other settings\n");
      printf ("     if integer value, run duration depends on value and other settings\n");
      printf ("       If estimating population phylogeny (-j0), # of phylogenies to save, 10 steps per save\n");
      printf ("       If fixed phylogeny in MCMC mode, # of genealogies to save\n");
      printf ("         This value times -d value sets the # of steps in the chain after burnin \n");
      printf ("       If fixed phylogeny,LOAD-GENEALOGY mode (-r0, -v), # genealogies to load from file(s)\n");
      printf ("	    If floating point: value = the time in hours between outputs. \n");
      printf ("         Run continues until file " "IMrun" " is absent from directory, or does not begin with 'y'\n");
      printf ("-m  Migration prior value (maximum for uniform,  mean if exponential distribution is used) \n");
      printf ("-o  Output file name (no spaces) default is 'outfile.txt' \n");
      printf ("-p  Output options: \n");
      printf ("    0 Turn off trend plots in outfile (default is to print trend plots)\n");
      printf ("    1 Turn off plots of marginal curves in outfile (default is to print marginal density plots)\n");
      printf ("    2 Print TMRCA histogram for each genealogy (MCMC mode only)\n");
      printf ("    3 Print histogram of splitting times divided by prior (do not use with -j0)\n");
      printf ("    4 Turn off printing estimates and histograms of population migration rate (2NM)\n");
      printf ("    5 Print pairwise probabilities that one parameter is greater than another \n");
      printf ("    6 Print histograms of the number of migration events (do not use with -j0)\n");
      printf ("    7 Print joint estimate for splitting times (not with -j0, models with 3, 4 or 5 populations)\n");
      printf ("-q  Maximum for population size parameters (4Nu) \n");
      printf ("-r  Run options \n");
      printf ("    0 LOAD-GENEALOGY Mode - load from previous run(s); also requires -v (do not use with -j0)\n");
      printf ("    1 Do not save genealogies (default saves sampled genealogies, ignored with -j0)\n");
      printf ("    2 Save the state of the Markov chain in a file - named with extension .mcf (MCMC mode only)\n");
      printf ("    3 Start by loading *.mcf file; requires -f (only the state space is loaded)\n");
      printf ("    4 Write all mutation related updates rates to stdout during the run (default is to suppress)\n");
      printf ("    5 Print burntrend file at end of burnin period\n");
      printf ("    6 When loading mcf files (-r3,-r7) do not load sampled values (i.e. use previous run as burnin)\n");
      printf ("    7 Load *.mcf file if present AND save to *.mcf file.\n");
      printf ("       No burn done if *.mcf is present. Use to continue a run by repeating the same command line.\n");
      printf ("    8 Use uniqpe IMburn/run files with outfilename as a prefix \n");
      printf ("-s  Random number seed (default is taken from current time)\n");
      printf ("-t  Maximum time of population splitting\n");
      printf ("-u  Generation time in years (default is %d) \n", GENERATIONTIMEDEFAULT);
      printf ("      Only needed if demographically scaled histograms are desired in fixed phylogeny runs\n");
      printf ("-v  Base name (no extension) of *.ti files with genealogy data  (requires use of -r0) \n");
      printf ("-w  Name of file with nested models to be tested (requires use of -r0), invokes �c2\n");
      printf ("-x  Set phylogeny priors for pairs of sampled populations (requires -j0)\n");
      printf ("      e.g. -x 0 1 0.5  set relative prior on all trees with 0,1 as sisters to 0.5 (default is 1)\n");
      printf ("      can have multiple -x terms in command line\n");
      printf ("-y  Mutation rate scalar for relevant loci\n");
      printf ("      Only needed if demographically scaled histograms are desired in fixed phylogeny runs\n");
      printf ("      Only needed in LOAD-GENEALOGY Mode (-r0) and if mutation rates not known for all loci\n");
      printf ("-z  Number of steps between screen output (default is %d) (don't use with LOAD-GENEALOGY mode -r0)\n",PRINTINTDEFAULT);
      fflush(stdout);
    }
#ifdef MPI_ENABLED
  MPI_Finalize();
#endif
    exit (0);
  }
  else
  {
/*
command line circumstances:
all flags begin with '-'
-most flags are single letter flags
-some are double letter flags: (e.g. jh)
-some flags are followed by a string or a character, others by an in int or a float
-All flags are followed by at least something,no flag is followed by nothing
it is ok to have spaces between a flag and its values

*/
    strcpy (command_line, "");
    strcpy (heatingterm_str,"");
    strcpy (modeloptions_str,"");
    strcpy (calcoptions_str,"");
    strcpy (outputoptions_str,"");
    strcpy (runoptions_str,"");
    strcpy (priors_str,"");
    for (i=0;i<HIDDENOPTIONNUMBER;i++) hiddenoptions[i] = 0;
    for (i=0;i<MODELOPTIONSNUMBER;i++) modeloptions[i]= 0;
    for (i=0;i<CALOPTIONSNUMBER;i++)calcoptions[i]=0;
    for (i=0;i<OUTPUTOPTIONSNUMBER;i++) outputoptions[i]=0;
    for (i=0;i<PRINTBURNTREND + 1;i++) runoptions[i]=0;

    for (i = 1; i < argc; i++)
    {
      strcat (command_line, argv[i]);
      strcat (command_line, " ");
    }
    for (i = 1; i < argc; i++)
    {
      //printf("cid %d %s\n",currentid,argv[i]);
      strcpy (pstr, argv[i]);
/*      if (toupper(pstr[1]) != 'X')
      {
        strcat (command_line, " ");
        strcat (command_line, pstr);
      } */

      if (strlen (pstr) < 2)
        IM_err (IMERR_COMMANDLINEFORMAT, " one of the command line strings is too short: %s ",pstr);
      if (pstr[0] != '-')
        IM_err (IMERR_COMMANDLINEFORMAT, "command line flag not preceded by '-' : %s", pstr);
      ch = toupper (pstr[1]);
      ch1 = ' ';
      if (ch == 'J')
      {
        strcat (modeloptions_str, " ");
        strcat (modeloptions_str, pstr);
        ch1 = toupper (pstr[2]);
      }
      if (ch == 'C')
      {
        strcat (calcoptions_str, " ");
        strcat (calcoptions_str, pstr);
      }
      if (ch == 'P')
      {
        strcat (outputoptions_str, " ");
        strcat (outputoptions_str, pstr);
      }
      if (ch == 'R')
      {
        strcat (runoptions_str, " ");
        strcat (runoptions_str, pstr);
      }
      if (ch == 'H')
      {
        strcat (heatingterm_str, " ");
        strcat (heatingterm_str, pstr);
        ch1 = toupper (pstr[2]);
      }
      if (strlen(pstr) > 1 && strchr(pstr+1,'-') != NULL)
        IM_err(IMERR_COMMANDLINEFORMAT,"command line format problem: %s",pstr);
      if (strlen (argv[i]) == 2 || (i < argc - 1 && isdigit (argv[i + 1][0])))  // space separates flag from its number
      {
        i++;
        strcpy (pstr, argv[i]);
      }
      else
      {
        if ((ch == 'H'))
          strdelete (pstr, 1, 3);
        else
          strdelete (pstr, 1, 2);
      }
      switch ((char) toupper (ch))
      {
       case 'B':
        tempf = atof (&pstr[0]);
        /* check to see if the value is floating point, in which case treat it as being in fractions of an hour  and convert to seconds */
        if (strchr (pstr, '.'))
        {
          burnduration = (int) (3600 * tempf);
          burndurationmode = TIMEINF;
	///AS: Should time be broadcast or should each processor have its own time?
          time (&lasttime);
#ifdef MPI_ENABLED
		        if (numprocesses > 1)
          {
            rc = MPI_Bcast(&lasttime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	         if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	         }
#endif

          runoptions[PRINTBURNTREND] = 1;
        }
        else
        {
          burnduration = (int) tempf;
          burndurationmode = TIMESTEPS;
          if (burnduration==0)
          {
            burndone = 1;
#ifdef MPI_ENABLED
		        /*rc = MPI_Bcast(&burndone, 1, MPI_INT, currentid, MPI_COMM_WORLD);
		        if (rc !=MPI_SUCCESS)
            {
			        MPI_Abort(MPI_COMM_WORLD,-1);
			        return;
		        } */
#endif
          }
        }
        Bp = 1;
        break;
      case 'C':
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
            IM_err (IMERR_COMMANDLINEFORMAT, "calculation option flag -c should be followed by a digit: %s",pstr);
          if (atoi (&pstr[j]) >= CALOPTIONSNUMBER )
              IM_err (IMERR_OPTIONNUMBERFAIL, "model option flag -c followed by option number out of allowed range : %d", atoi (&pstr[j]));
          calcoptions[atoi (&pstr[j])] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
      case 'D':
        samplemcmcinterval = atoi (&pstr[0]);  //samplemcmcinterval is just a local variable. used befow depending on whether -j0 or not
        Dp = 1;
        break;
      case 'F':
        strcpy (mcfreadfilename, pstr);
        put_spaces_in_filepaths(mcfreadfilename);
        Fp = 1;
        break;
      case 'G':
        strcpy (priorfilename, pstr);
        put_spaces_in_filepaths(priorfilename);
        Gp = 1;
        break;
      case 'H':
        switch ((char) toupper (ch1))
        {
        case 'A':
          hval1 = atof (&pstr[0]);
          Hap = 1;
          break;
        case 'B':
          hval2 = atof (&pstr[0]);
          Hbp = 1;
          break;
        case 'N':
          numchainstotal = atoi (&pstr[0]);
          /* set and check numchainspp */
          numchainspp = numchainstotal/numprocesses;
          Hnp = 1;
          break;
      /*  case 'K':  jh no longer used as of 6/16/2016  now using setswaptries()
          swaptries = atoi (&pstr[0]);   // this is Bcast down further down in this function
          Hkp = 1;
          break; */
        case 'F':
          Hfp = 1;
          switch ((char) toupper (pstr[0]))
          {
          case 'G':
            heatmode = HGEOMETRIC;  // default
            break;
          case 'S':
            heatmode = HFULL;  //JH added to deal with hidden genealogies and topology updating // not really in use as of 6/2017
            break;
          default:
            heatmode = HLINEAR;
            break;
          }
          break;
        default:
          IM_err (IMERR_COMMANDLINEHEATINGTERMS, "mistake in use of -h flag : %s", pstr);
        }
        break;
      case 'I':
        infilename = strdup (pstr);
        put_spaces_in_filepaths(infilename);
	       //strncpy(infile_name, infilename, sizeof(infilename));	 //don't think sizeof makes sense here.
        strcpy(infile_name, infilename);	  // should be safe
        Ip = 1;
        break;
      case 'J':
        Jp = 0;
        if (!(toupper(ch1) == 'H'))
        {
          j = (int) (strlen (pstr) - 1);
          while (j >= 0)
          {
            if (!isdigit (pstr[j]))
            {
              if (toupper(pstr[j]) == 'X')  // 10/30/2017  ran out of digits for modeloptions, awkward
              {
                modeloptions[ONEPOPSIZEPARAMETER] = 1;
              }
              else
                IM_err (IMERR_COMMANDLINEFORMAT, "model option flag -j should be followed by a digit: %s", pstr);
            }
            //if (atoi (&pstr[j] ) >= MODELOPTIONSNUMBER  )  can't use this when we have 10 options and using 'X'
              //IM_err (IMERR_OPTIONNUMBERFAIL, "model option flag -j followed by option number out of allowed range : %d", atoi (&pstr[j] ));
            else
              modeloptions[atoi (&pstr[j])] = 1;
            pstr[j] = '\0';
            j--;
          }
        }
        else // -jh  hiddenoptions[]
        {
          if (isdigit(pstr[1])  && pstr[1] == '1') // hiddenoptions[WRITEMIGRATIONNAME]
          {
            char tempc = pstr[2];
            migrationnamefrom = atoi (&tempc);
            tempc = pstr[3];
            migrationnameto = atoi (&tempc);
            strcpy(migrationnamefilename, &pstr[4]);
            put_spaces_in_filepaths(migrationnamefilename);
            hiddenoptions[WRITEMIGRATIONNAME] = 1;
          }
          else
          {
            j = (int) (strlen (pstr) - 1);  // position of last charater of the string
            while (j >= 0)
            {
              if (isdigit(pstr[j]))
              {
                if (atoi (&pstr[j] ) >= HIDDENOPTIONNUMBER  )
                    IM_err (IMERR_OPTIONNUMBERFAIL, "model option flag -jh followed by option number out of allowed range : %d", atoi (&pstr[j] ));
                hiddenoptions[atoi (&pstr[j])] = 1;
              }
              else
              {
                if (toupper(pstr[j])=='A') // treat 'A' as SKIPMOSTUSCALAROUTPUT, can be necessary if # of digits gets too high for hiddenoptions
                  hiddenoptions[SKIPMOSTUSCALAROUTPUT] = 1;
                if (toupper(pstr[j])=='B') // treat 'B' as STOPMOSTINTERVALOUTPUT
                  hiddenoptions[STOPMOSTINTERVALOUTPUT] = 1;
                if (toupper(pstr[j])=='C') // treat 'C' as READOLDMCFFILE  //jh1_17_2018
                  hiddenoptions[READOLDMCFFILE] = 1;
              }
              pstr[j] = '\0';
              j--;
            }
          }
        };
        break;
      case 'L':
        tempf = atof (&pstr[0]);
        /* check to see if the value is floating point, in which case treat it as being in fractions of an hour  and convert to seconds */
        if (strchr (pstr, '.'))
        {
          chainduration = (int) (3600 * tempf);
          cdurationmode = TIMEINF;
          genealogiestosave = -1;
          phylogeniestorecord = -1;  // only used if modeloptions[POPTREETOPOLOGYUPDATE] == 1
        }
        else
        {
          numthingstosave = (int) tempf;
          genealogiestosave = (int) tempf;
          phylogeniestorecord = genealogiestosave;  // only used modeloptions[POPTREETOPOLOGYUPDATE] == 1,  then genealogiestosave gets reset to -1 below
          cdurationmode = TIMESTEPS;
        }
        Lp = 1;
        break;
      case 'M':
        strcat (priors_str, " -m ");
        strcat (priors_str, pstr);
        mprior = (double) atof (&pstr[0]);
        if (mprior == 0)
          modeloptions[NOMIGRATION] = 1;
        Mp = 1;
        break;
      case 'O':
        outfilename = strdup (pstr);
        put_spaces_in_filepaths(outfilename);
        Op = 1;
        break;
      case 'P':
        Pp = 1;
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
          {
            IM_err (IMERR_COMMANDLINEFORMAT, "print option flag -p should be followed by a digit : %s",pstr);
          }
          k = atoi (&pstr[j]);
          if (k >= OUTPUTOPTIONSNUMBER  )
              IM_err (IMERR_OPTIONNUMBERFAIL, "model option flag -p followed by option number out of allowed range : %d", k);
          outputoptions[k] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
      case 'Q':
        strcat (priors_str, " -q ");
        strcat (priors_str, pstr);
        thetaprior = (double) atof (&pstr[0]);
        Qp = 1;
        break;
      case 'R':
        Rp = 1;
        j = (int) (strlen (pstr) - 1);
        while (j >= 0)
        {
          if (!isdigit (pstr[j]))
            IM_err (IMERR_COMMANDLINEFORMAT, "run option flag -r should be followed by a digit : %s ", pstr);
          k = atoi (&pstr[j]);
          if (k >= RUNOPTIONSNUMBER )
              IM_err (IMERR_OPTIONNUMBERFAIL, "model option flag -r followed by option number out of allowed range : %d", k);
          runoptions[k] = 1;
          pstr[j] = '\0';
          j--;
        }
        break;
      case 'S':
        seed_for_ran1 = atoi (&pstr[0]);
	      seed_for_ran1 = seed_for_ran1 * (1+currentid); //jh changed to 1+currentid ///Just to make sure that there's a different seed on each process
        if (!seed_for_ran1)
          seed_for_ran1 = currentid;
        Sp = 1;
        break;
      case 'T':
        strcat (priors_str, " -t ");
        strcat (priors_str, pstr);
        tprior = (double) atof (&pstr[0]);
        Tp = 1;
        break;
      case 'U':
        generationtime = atof (&pstr[0]);
        Up = 1;
        break;
      case 'V':
        loadfilebase = strdup (pstr);
        put_spaces_in_filepaths(loadfilebase);
        Vp = 1;
        break;
      case 'W':
        strcpy (nestedmodelfilename, pstr);
        put_spaces_in_filepaths(nestedmodelfilename);
        calcoptions[FINDJOINTPOSTERIOR] = 1;
        break;
      case 'Y':
        scaleumeaninput = atof (&pstr[0]);
        Yp = 1;
        break;
      case 'X':
          strcpy (pstr, argv[i]);
          strcat(topologypriorinfostring,pstr);
          strcat(topologypriorinfostring," ");
          i += 1;
          strcat(topologypriorinfostring,argv[i]);
          strcat(topologypriorinfostring," ");
          i += 1;
          strcat(topologypriorinfostring,argv[i]);
          strcat(topologypriorinfostring," ");
          Xp = 1;
        break;
      case 'Z':
        printint = atoi (&pstr[0]);
        Zp = 1;

        break;
      /*case 'A':   // not in use for command line
        tempnumgupdates = atoi (&pstr[0]);
        Ap = 1;
        break; */
      default:
        IM_err (IMERR_COMMANDLINEFORMAT, &ch);
      }
    }
  }
  if (infilename==0)
    IM_err (IMERR_READFILEOPENFAIL,  "pointer to input file not set,  check -i on command line");
  npops = imaInfileNpops (infilename);
  if (npops < 1 || npops > 10)
  {
    IM_err (IMERR_COMMANDLINEFORMAT, "Number of populations must be nonnegative and less than 10");
  }
  if (npops == 1)
  {
    if (modeloptions[NOMIGRATION] == 0)
      modeloptions[NOMIGRATION] = 1;
    if (Mp == 1)
    {
      IM_err (IMERR_COMMANDLINEFORMAT, "model option [single population] should not go with -m");
    }
    if (Tp == 1)
    {
      IM_err (IMERR_COMMANDLINEFORMAT, "model option [single population] should not go with -t");
    }
  }
  if (numchainstotal != numchainspp * numprocesses)
  {
      IM_err (IMERR_CHAINNUM, "Error - user specified number of chains (%d) is not evenly divisible by the number of CPUs (%d)",numchainstotal,numprocesses);
  }
  if (numprocesses > 1 && numchainspp < MINNUMCHAINSPERPROCESSOR)
  {
    IM_err (IMERR_CHAINNUM, "Error - user specified number of chains (%d) divided by the number of CPUs (%d) leaves too few chains per processor (%d) - at least %d required",numchainstotal,numprocesses,numchainspp,MINNUMCHAINSPERPROCESSOR);
  }
  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    hiddenoptions[HIDDENGENEALOGY] = 1;
  if (!Ip)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " No data file given on command line");
  }
  if (!Lp && !runoptions[LOADRUN])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " No run duration (-L) given on command line");
  }
  if (runoptions[LOADRUN] || modeloptions[NOMIGRATION])
    outputoptions[MIGRATEHIST] = 0;
  if (!Op)
    //strcpy (outfilename, "outfile.txt");
    outfilename = strdup ("outfile.txt");
#ifdef DEBUG
  if (currentid == HEADNODE)
    checkoutfileclosed (&outfile, outfilename);   // just make sure that outfilename does not name a file that is already opened
#endif
  if (runoptions[LOADRUN] && !Vp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO, " -r0 invoked without -v information, i.e. no base name for files containing genealogys was given on the command line");
  }
  /* 5/19/2011 JH adding thermodynamic integration  - only the likelihood ratio gets raised to beta,  not the prior ratio */
  /*if (calcoptions[CALCMARGINALLIKELIHOOD])   this section not really relevant as not clear if marginallikelihood calculations are working
  {
    if (!Hnp)  // no multiple chains set on command line
      IM_err (IMERR_MISSINGCOMMANDINFO, " marginal likelihood calculations invoked (-c4) but multiple chains (-hn) not specified on command line\n");
    if (runoptions[LOADRUN])
    {
      IM_err (IMERR_COMMANDLINEINCOMPAT, " Conflicting command line arguments, cannot estimate marginal likelihood in Load mode (-r0)");
    }
    if (!Hfp)  // no heat mode set on command line
      heatmode = HFULL;
    else
    {
      if (heatmode != HFULL)
        IM_err (IMERR_MISSINGCOMMANDINFO, " marginal likelihood calculations invoked (-c4) but heating model used (-hf) not appropriate. Use -hfs or even heating (default)\n");
    }
  } */
  if (numchainspp > 1)
  {
    if ( calcoptions[CALCMARGINALLIKELIHOOD] && heatmode != HFULL)
      IM_err (IMERR_COMMANDLINEHEATINGTERMS, "wrong heating mode.  -hf s required when calculating marginal likelihood ");
    if ( (heatmode == HGEOMETRIC && (numchainspp * numprocesses) < 4) || (heatmode == HFULL && (numchainspp * numprocesses) < 10)  )
    {
      IM_err (IMERR_COMMANDLINEHEATINGTERMS, "too few chains specified in heating model");
    }
    else if (!Hfp && calcoptions[CALCMARGINALLIKELIHOOD]==0)  // no heat mode set on command line
    {
      heatmode = HGEOMETRIC;  //setting default to HGEOMETRIC 6/15/2017
    }
    else  // heat mode given on command line
    {
      if (heatmode > HLINEAR)
      {
        if (heatmode == HGEOMETRIC)
        {
          if (!Hap) {
            hval1 = 0.95;  // default value
	        }
          if (!Hbp) {
            hval2 = 0.8;
	          }
          // if (hval1 > 1.0)  //6/11/2010 JH  stopped this,  it turns out numbers slightly higher than 1 can be useful when the are
            // a large number of chains
            // IM_err (IMERR_COMMANDLINEHEATINGTERMS, "ha commandline term is out of range, should be <= 1");
          if (hval1 > 1.1) //6/11/2010  JH  added this to avoid values much larger than 1
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "for geometric heating it is not useful to have the ha term be greater than 1.1");
          if (hval1 < 0.9)
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "for geometric heating it is not useful to have the ha term be less than 0.9");
          if (hval2 >= 1.0|| hval2 <= 0.0)
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "hb commandline term is out of range, should be < 1 and > 0)");
        }
        if (heatmode == HFULL)   //JH added to deal with hidden genealogies and topology updating
        {
          if (!Hap  ||  hval1 > 1.0 || hval1 <= 0.0)
              IM_err (IMERR_COMMANDLINEHEATINGTERMS, "for full heating model, -ha term not specified or out of range");
          if (Hbp == 1)
            IM_err (IMERR_COMMANDLINEHEATINGTERMS, "for full heating model, -hb should not be specified");
        }
      }
      else if (!Hap)
        hval1 = 0.05;           /* default value */
	    }
  }
  /* setting  mcmc step intervals:
    regardless of sampling phylgoenies or not
      default mcmc checking is RECORDINTERVALDEFAULT 10
      this applies to everything an ESS can be calcualted on
      can't mess with this because it might mess up calculating autocorrelations and ESS values
    if NOT sampling phylogenies:
      default genealogy sampling is SAMPLEGENEALOGYINTERVALDEFAULT  100
        savegenealogyint = SAMPLEGENEALOGYINTERVALDEFAULT
      -d gives the interval for sampling genealogies
      mcmc checking is set to the default RECORDINTERVALDEFAULT 10
    if sampling phylogenies:
      genealogies are not sampled
        savegenealogyint = -1
  */
  recordint = RECORDINTERVALDEFAULT;   // this is fixed and is not a command line option
  if (modeloptions[POPTREETOPOLOGYUPDATE]==0)
  {
    if (!Dp)
      savegenealogyint = SAMPLEGENEALOGYINTERVALDEFAULT;
    else
      savegenealogyint = samplemcmcinterval;
  }
  else
  {
    savegenealogyint = -1;   // do not save genealogies
  }

  if (!Up) {
    generationtime = GENERATIONTIMEDEFAULT;
    usegenerationtimedefault = 1;
  }
  else
    usegenerationtimedefault = 0;
  if (!Zp)
  {
    printint = PRINTINTDEFAULT;
	 }
  if (!Bp && !runoptions[LOADRUN] && !runoptions[SAVELOADSAMEMCFFILE])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            "No burn duration information given on command line (use -b)");
  }
  if (modeloptions[NOMIGRATION] && !Mp)
  {
    mprior = 0;
	 }
  if (modeloptions[NOMIGRATION] && Mp && mprior > 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT, " Conflicting command line arguments, no migration set but migration prior > 0 : %lf",mprior);
  }
  if (Gp && !calcoptions[LOADPRIORSFROMFILE])
  {
    calcoptions[LOADPRIORSFROMFILE] = 1;
  }

  if (!Mp && modeloptions[NOMIGRATION] != 1 && !calcoptions[LOADPRIORSFROMFILE] && npops > 1)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value for migration parameter (-m)");
  }
  if (!Qp && !calcoptions[LOADPRIORSFROMFILE])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value of 4Nu parameters (-q)");
  }
  if (!Tp && !calcoptions[LOADPRIORSFROMFILE])
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            " No information provided for maximum value of t parameters (-t)");
  }
  if (runoptions[LOADMCSTATE] && !Fp)
  {
    IM_err (IMERR_MISSINGCOMMANDINFO,
            "  -r3 invoked without -f information, i.e. no filename given for markov chain state file,  for loading state of markov chain ");
  }
  /* not longer set swaptries on command line.  now using setswaptries()   as of 6/16/2016 */
  swaptries_per_cpu = setswaptries();
	 if (runoptions[PRINTBURNTREND])
  {
    if (runoptions[LOADMCSTATE])
      burntrendstartdelay = 0;
    else
      burntrendstartdelay = BURNTRENDSTARTDELAYDEFAULT;
  }
	 if (!Sp)
  {
    seed_for_ran1 = (long) time (NULL);
  }
  if (strcmp (infilename, outfilename) == 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT, " Input and output file names are identical");
  }

  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
  {
    if (hiddenoptions[GSAMPINFOEXTRA] == 0)  // using GSAMPINFOEXTRA causes a .ti file with added info
        runoptions[DONTSAVEGENEALOGIES] = 1;
    else
        runoptions[DONTSAVEGENEALOGIES] = 0;
    if (npops <= 2)
      IM_err (IMERR_COMMANDLINEINCOMPAT," 3 or more populations required for population topology updating");
    if (modeloptions[NOMIGRATION])
      IM_err (IMERR_COMMANDLINEINCOMPAT," Migration must be included in the model for population topology updating");
    if (runoptions[LOADRUN] ||
        modeloptions[PARAMETERSBYPERIOD] ||
        modeloptions[MIGRATIONBETWEENSAMPLED] ||
        modeloptions[NOMIGBETWEENNONSISTERS] )
        IM_err (IMERR_COMMANDLINEINCOMPAT," Model options incompatible with population topology updating");
    if (runoptions[LOADRUN])
      IM_err (IMERR_COMMANDLINEINCOMPAT," Run options incompatible with population topology updating");
    if (calcoptions[FINDJOINTPOSTERIOR] )//|| calcoptions[LOADPRIORSFROMFILE])  // cannot have complex priors with phylogeny updating, as priors need to change with tree
      IM_err (IMERR_COMMANDLINEINCOMPAT," Calculation options incompatible with population topology updating");
    if ( outputoptions[PRINTTMRCA]||
        outputoptions[THISTDIVIDEBYPRIOR]||
        outputoptions[MIGRATEHIST]||
        outputoptions[PRINTJOINTTEST])
      IM_err (IMERR_COMMANDLINEINCOMPAT," Output options incompatible with population topology updating");
  }
  else if (Xp != 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT,"Population prior terms (-x) used without invoking topology updating (-j0)\n");
  }
  /*if (hiddenoptions[HIDDENGENEALOGY] && calcoptions[CALCMARGINALLIKELIHOOD])
  {
    IM_err(IMERR_COMMANDLINEINCOMPAT," Population updating (or any use of hidden genealogies) is not compatible with marginal likelihood calculation (-c4)");
  }*/
  if (hiddenoptions[WRITEMIGRATIONNAME])
  {
    migrationnamefile = fopen (migrationnamefilename, "w");
  }
  if (( modeloptions[NOMIGBETWEENNONSISTERS] ||
        modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] ||
        modeloptions[MIGRATIONBETWEENSAMPLED] ||
        //modeloptions[ADDGHOSTPOP] ||
        modeloptions[PARAMETERSBYPERIOD] ||
       /* modeloptions[EXPOMIGRATIONPRIOR] ||  */
        npops == 1 ||
        modeloptions[NOMIGRATION]) && calcoptions[LOADPRIORSFROMFILE])
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT,"incompatibility between a model option and the use of a file with parameter priors");
  }
  if (( modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] ||
        npops == 1 ||
        modeloptions[NOMIGRATION]) &&  modeloptions[ONEMIGRATIONPARAMETER])
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT,"incompatibility between a model option and the use of a single migration parameter");
  }
  if ( modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] && outputoptions[NOPOPMIGPARAMHIST] == 0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT,"incompatibility between a output option for population migration rate histogram and model option a single migration parameter");
  }
  if (calcoptions[FINDJOINTPOSTERIOR]==1 && runoptions[LOADRUN]==0)
  {
    IM_err (IMERR_COMMANDLINEINCOMPAT," finding joint posterior and tests of nested models requires L mode");
  }

  if (runoptions[LOADRUN] == 1)
  {
    cdurationmode = TIMESTEPS;
    if (genealogiestosave < 0)
      genealogiestosave = (int) DEFAULTNUMGENEALOGIES;
  }
  else if (cdurationmode == TIMESTEPS)
  {
    if (modeloptions[POPTREETOPOLOGYUPDATE] == 0)
    {
      genealogiestosave = numthingstosave;
      phylogeniestorecord = -1;
      chainduration = genealogiestosave * savegenealogyint;
    }
    else
    {
      if (hiddenoptions[GSAMPINFOEXTRA]==1)
      {
        phylogeniestorecord = numthingstosave;
        genealogiestosave = numthingstosave;
        chainduration = phylogeniestorecord * recordint;
      }
      else
      {
        phylogeniestorecord = numthingstosave;
        genealogiestosave = -1;
        chainduration = phylogeniestorecord * recordint;
      }
    }
	 }
  if (calcoptions[LOADPRIORSFROMFILE] && !Gp)
  {

    strcpy(priorfilename,defaultpriorfilename);
  }

 if (!Gp && (modeloptions[POPTREETOPOLOGYUPDATE]==0 && modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]))
 {
     sprintf(priorfilename,"%s.%s",outfilename,defaultpriorfilename);
 }

if (modeloptions[POPTREETOPOLOGYUPDATE]==0 && modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]==1)
  calcoptions[LOADPRIORSFROMFILE] = 0;
  const int ucbuffer = 200;
  if (burndurationmode == TIMEINF) {
      std::stringstream s("");
      if (runoptions[UNIQUEBURNRUN]) {
          s << outfilename << ".IMburn";
      } else {
          s << "IMburn";
      }
      burnfilename = (char *) malloc(s.str().size()+1);
      strcpy(burnfilename,s.str().c_str());
  }
  if (cdurationmode == TIMEINF) {
      std::stringstream s("");
      if (runoptions[UNIQUEBURNRUN]) {
          s << outfilename << ".IMrun";
      } else {
          s << "IMrun";
      }
      runfilename = (char *) malloc(s.str().size()+1);
      strcpy(runfilename,s.str().c_str());
  }

 return;
}                               // scan_commandline

/*  this prints basic info to a string, fpstr, that later gets printed to the output file */
void
begin_outputfile_info_string (void)
{
  const char *buildtimestring = "IMa3 program compiled on " __DATE__ ", " __TIME__ ".";
  fpstri = static_cast<int *> (malloc (sizeof (int)));
  *fpstri = 0;
  SP "IMa3 - Isolation with Migration Analysis  -  Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Arun Sethuraman, Yujin Chung 2018 \n");

#ifdef STDTEST
  SP "\n\n==================\nSTANDARD TEST MODE\n");
#ifdef RANDOM_NUMBERS_FROM_FILE
  SP "USES randNums FILE\n");
#endif
  SP "==================\n\n");
#endif
  SP "\n%s\n",releaseinfostring());
  SP "%s\n",buildtimestring);
#ifdef TURNONCHECKS
  SP "\n**Compiled and Run with TURNONCHECKS**\n");
#endif //TURNONCHECKS

  SP "\nJob Started: %s\n",timeinfostring);
  SP "%s",outputbanner("input and starting information"));


  if (calcoptions[DONTCALCLIKELIHOODMUTATION])
    SP "\n\n**NO DATA ** - Data likelihoods set to constant  posterior should equal prior \n\n");
  SP "\nCommand line string : %s \n", command_line);
  SP "\nRun Settings:\n");
  SP "  Input filename : %s \n", infilename);
  SP "  Output filename: %s \n", outfilename);
  if (calcoptions[LOADPRIORSFROMFILE])
    SP "  Prior distribution filename: %s \n",priorfilename);



  if (runoptions[LOADMCSTATE])
  {
    SP "  LOAD MCF: Initial Markov chain state space loaded, filename: %s\n",  mcfreadfilename);
  }
  if (runoptions[SAVEMCSTATEFILE])
  {
    SP "  SAVE MCF: Final Markov chain state space saved, filename (1 file per cpu used in the run): %s\n", mcfwritefilename);
  }
  if (runoptions[SAVELOADSAMEMCFFILE])
  {
    SP "  SAVE/LOAD MCF: Markov chain state space loaded from, and saved to %s\n",mcfwritefilename);
  }

  if (runoptions[DONTSAVEGENEALOGIES]==0)
  {
    SP "  Saved genealogy sample filename: %s\n", genealogyinfosavefilename);
  }
  if (outputoptions[MIGRATEHIST])
  {
    SP "  Saved distributions of migration event counts filename: %s%s\n",  outfilename,".mpt");
  }
    if (usegenerationtimedefault)
    SP "\n  Generation Time set to DEFAULT :  1 year\n");
  else
    SP "\n  Generation Time (years): %.2lf\n",generationtime);
  SP "\n");
  if (calcoptions[LOADPRIORSFROMFILE]==0)
    SP "  Prior distribution terms on command line : %s \n\n",priors_str);
  SP "  Calculation options on command line : %s \n",calcoptions_str);
  SP "  Model options on command line : %s \n",modeloptions_str);
  SP "  Run options on command line : %s \n",runoptions_str);
  SP "  Output options on command line : %s \n",outputoptions_str);
  SP "\n");
  if (strlen(topologypriorinfostring) > 0)
    SP "  Population tree topology sister group priors on command line: %s \n",topologypriorinfostring);

  if (modeloptions[EXPOMIGRATIONPRIOR]==1)
  {
    SP"  Exponential priors used for migration rate parameters\n");
  }

  if (hiddenoptions[NOMUTATIONSCALARUPATES]==1)
    SP"  Mutation rate scalars fixed at 1.  No updates.\n");
  if (hiddenoptions[FIXMUTATIONSCALARUPATES]==1)
    SP"  Mutation rate scalars fixed at relative values of mutation rate per year,  given on command line. No updates.\n");

  if (outputoptions[PRINTTMRCA])
    SP "  TMRCA  histograms printed \n");

  if (!runoptions[LOADRUN])
  {
    SP "\n  Run Duration:\n  -------------\n");
    switch (burndurationmode)
    {
    case TIMESTEPS:
      SP "  Burnin period (-B), # steps: %li ", burnduration);
      break;
    case TIMEINF:
      SP "  Burnin period (-B), # seconds: %li (total burn duration depends on IMburn file status)", burnduration);
      break;

    };
    if (runoptions[SAVELOADSAMEMCFFILE])
    {
      SP " (If mcf file(s) loaded using -r%d  no burn will be done)\n",(int) SAVELOADSAMEMCFFILE);
      noburn_mcfload = 1;
    }
    else
      SP "\n");
    if (runoptions[PRINTBURNTREND])
    {
      SP "   User option for printing trendline during, or at end of burnin period, invoked\n");
      SP "   Initial burn duration prior to beginning recording burntrend : %d steps\n", burntrendstartdelay);
    }

    switch (cdurationmode)
    {
    case TIMESTEPS:
      {
        if (genealogiestosave > 0)
          SP "  Record period (-L), #genealogies saved: %d  #steps each: %li   total #steps: %li \n", genealogiestosave, chainduration / genealogiestosave, chainduration);
      if (phylogeniestorecord > 0)
          SP "  Record period (-L), #phylogenies recorded: %d  #steps each: %li   total #steps: %li \n", phylogeniestorecord, chainduration / phylogeniestorecord, chainduration);
      }
      break;
    case TIMEINF:
      SP "  Record period (-L), # seconds: %ld \n",
        chainduration);
      SP "  Run period indefinite in length - determined by user using 'IMrun' file\n");
      break;
    }
    SP "\n  Metropolis Coupling:\n  --------------------\n");
    if (numchainstotal > 1)
    {
      SP "  Heating terms on command line : %s \n",heatingterm_str);
      //SP "     Metropolis Coupling implemented using %d chains (%d processor(s) with %d chains per processor)\n", numchainstotal,numprocesses,numchainspp);
      SP "  Metropolis Coupling implemented using %d chains\n", numchainstotal);
      SP "  Number of processors: %d",numprocesses);
      if (numprocesses > 1)
        SP " (%d chains per processor)\n",numchainspp);
      else
        SP "\n");
      switch (heatmode)
      {
      case HLINEAR:
        SP "  Heating Model: Linear   term: %.3f\n", hval1);
        break;
      case HGEOMETRIC:
        SP "  Heating Model: Geometric   term1: %.3f  term2: %.3f\n",
          hval1, hval2);
        break;
      case HFULL:    //JH added to deal with hidden genealogies and topology updating
        SP "  Heating Model: Full   term: %.3f\n", hval1);
        break;
      }
      SP "  Number of chain swap attempts per step: %d\n",swaptries_per_cpu);
//      if (numprocesses > 1)
  //      SP "  (# within each processor: %d, # between processors: %d\n",withinprocswaptries,betweenprocswaptries);
    }
    else
      SP "  None \n");
#ifndef RANDOM_NUMBERS_FROM_FILE
    SP "\n  Random number seed:\n  -------------------\n  %li \n\n", seed_for_ran1);
#endif
  }

}                               //  begin_outputfile_info_string


void set_mode_check_inconsistent_options(void)
{
 // set runmode between 0 and NUMRUNMODES
 if (modeloptions[POPTREETOPOLOGYUPDATE])
  {
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      runmode = POPTREEHYPERPRIORmode0; //0
    else
      runmode = POPTREEmode1; //1
  }
  else
  {
     if (hiddenoptions[HIDDENGENEALOGY]==1)
     {
       if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
        runmode = HGHYPERPRIORmode5;//5
       else
         runmode = HGmode6; //6
     }
     else
     {
       if (runoptions[LOADRUN])
         runmode = LOADGmode4; //4
       else
       {
         if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
          runmode = GHYPERPRIORmode2; //2
        else
          runmode = Gmode3; //3
       }
     }
  }
 switch (runmode)
 {
 case POPTREEHYPERPRIORmode0:
   {
     if (calcoptions[FINDJOINTPOSTERIOR] || calcoptions[LOADPRIORSFROMFILE])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some calcoptions (-c) incompatible with mode 0");
     if (modeloptions[NOMIGBETWEENNONSISTERS] || modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] || modeloptions[MIGRATIONBETWEENSAMPLED] || modeloptions[PARAMETERSBYPERIOD] || modeloptions[NOMIGRATION] || modeloptions[ONEMIGRATIONPARAMETER])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some modeloptions (-j) incompatible with mode 0");
     if(outputoptions[MIGRATEHIST])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some outputoptions (-p) incompatible with mode 0");
     if (runoptions[LOADRUN])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some runoptions (-r) incompatible with mode 0");
     break;
   }
 case POPTREEmode1:
   {
     if (calcoptions[FINDJOINTPOSTERIOR] || calcoptions[LOADPRIORSFROMFILE])
       IM_err (IMERR_COMMANDLINEINCOMPAT," calcoptions (-c) incompatible with mode 1");
    if (modeloptions[NOMIGBETWEENNONSISTERS] || modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] || modeloptions[MIGRATIONBETWEENSAMPLED] || modeloptions[PARAMETERSBYPERIOD] || modeloptions[NOMIGRATION] || modeloptions[ONEMIGRATIONPARAMETER])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some modeloptions (-j) incompatible with mode 1");
     if(outputoptions[MIGRATEHIST])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some outputoptions (-p) incompatible with mode 1");
     if (runoptions[LOADRUN])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some runoptions (-r) incompatible with mode 1");

     break;
   }
 case GHYPERPRIORmode2:
   {
     if (calcoptions[FINDJOINTPOSTERIOR] || calcoptions[LOADPRIORSFROMFILE])
       IM_err (IMERR_COMMANDLINEINCOMPAT," calcoptions (-c) incompatible with mode 2");
     if (modeloptions[NOMIGBETWEENNONSISTERS] || modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] || modeloptions[MIGRATIONBETWEENSAMPLED] || modeloptions[PARAMETERSBYPERIOD] || modeloptions[POPTREETOPOLOGYUPDATE] || modeloptions[ONEMIGRATIONPARAMETER])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some modeloptions (-j) incompatible with mode 2");
     if (runoptions[LOADRUN])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some runoptions (-r) incompatible with mode 2");

    break;
   }
 case Gmode3:
   {

    if (modeloptions[POPTREETOPOLOGYUPDATE] || modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some modeloptions (-j) incompatible with mode 3");
     if (runoptions[LOADRUN])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some runoptions (-r) incompatible with mode 3");
     break;
   }
 case LOADGmode4:
   {
    if (modeloptions[POPTREETOPOLOGYUPDATE] || modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some modeloptions (-j) incompatible with mode 4");
     if (runoptions[PRINTMUTATIONUPDATESTOSCREEN] || runoptions[PRINTBURNTREND])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some runoptions (-r) incompatible with mode 4");

     break;
   }
 case HGHYPERPRIORmode5:
   {
     if (calcoptions[FINDJOINTPOSTERIOR] || calcoptions[LOADPRIORSFROMFILE])
       IM_err (IMERR_COMMANDLINEINCOMPAT," calcoptions (-c) incompatible with mode 5");
    if (modeloptions[NOMIGBETWEENNONSISTERS] || modeloptions[SINGLEMIGRATIONBOTHDIRECTIONS] || modeloptions[MIGRATIONBETWEENSAMPLED] || modeloptions[PARAMETERSBYPERIOD] || modeloptions[POPTREETOPOLOGYUPDATE] || modeloptions[ONEMIGRATIONPARAMETER])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some modeloptions (-j) incompatible with mode 5");
     if (runoptions[LOADRUN])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some runoptions (-r) incompatible with mode 5");
     break;
   }
 case HGmode6:
   {
    if (modeloptions[POPTREETOPOLOGYUPDATE] || modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some modeloptions (-j) incompatible with mode 6");
     if (runoptions[LOADRUN])
       IM_err (IMERR_COMMANDLINEINCOMPAT," some runoptions (-r) incompatible with mode 6");
     break;
   }
 default:
   break;
 }
} // set_mode_check_inconsistent_options(void)


void add_previous_run_info_to_output(int *fpstri, char fpstr[])
{
  if (numpriormcfruns > 0)
  {
    SP "\nPrevious Run Info (using mcf files) \n");
    SP "-----------------------------------\n");
    SP "  Number of prior runs: %d\n",numpriormcfruns);
    SP "  Total time of prior runs : %s\n",timestring(totaltime));
    /*SP "  Total time of prior runs : %d hours, %d minutes, %d seconds \n",
        (int) totaltime / (int) 3600,
          ((int) totaltime / (int) 60) - ((int) 60 * ((int) totaltime / (int) 3600)),
              (int) totaltime - (int) 60 *((int) totaltime / (int) 60)); */
    SP "\n");
  }
}

void start (int argc, char *argv[], int currentid)
{
  int i;
  FILE *checkfile;
  int rc = 0;

  thetaprior = mprior = tprior = -1.0;
  scan_commandline (argc, argv,currentid);
  int mcf0found,mcffound;
  set_mode_check_inconsistent_options(); // set runmode
  if (runmode == GHYPERPRIORmode2 || runmode == HGHYPERPRIORmode5)
    runoptions[DONTSAVEGENEALOGIES] = 1;


  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
  {
    totaltopolupdates = 0;
    chain0topolupdates = 0;
    chain0topolswaps = 0;
    poptopologiessampled = 0;
    poptopologyproposedlist = static_cast<int *> (calloc ( (size_t) numtreesarray[npops],sizeof (int)));  // initialize at 0
  }

  if (runoptions[DONTSAVEGENEALOGIES]==0)
  {
	  if (currentid == HEADNODE)
   {
       sprintf(genealogyinfosavefilename, "%s.ti", outfilename); // would not hurt to have other nodes set for this, why only on currentid ?
	  }
   genealogysamples = 0; //not sure if this needed here
  }
  if (cdurationmode == TIMEINF)
  {
    sprintf(oldoutfilename, "%s.old", outfilename);
  }
  setseeds (seed_for_ran1 * (currentid+1));  // JH bad bug,  change from currentid  to currentid+1
  setlogfact ();
  setheat (hval1, hval2, heatmode, currentid);
  if (runoptions[SAVEMCSTATEFILE] || runoptions[SAVELOADSAMEMCFFILE])
  {
    sprintf(mcfwritefilename, "%s.mcf.%d", outfilename, currentid);
  }
  begin_outputfile_info_string ();
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    if (modeloptions[EXPOMIGRATIONPRIOR])
    {
      if (mprior < MINPRIORFROMHYPERPRIOR)
        IM_err (IMERR_HYPERPROB," Migration hyperprior value %.4lg  too low,  can't be less than %.4lg",mprior, (double) MINPRIORFROMHYPERPRIOR);
      hyperprior_expo_m_mean = expo_m_mean = mprior;
      hyperprior_uniform_m_max = -1.0; // i.e. not in use
    }
    else
    {
      if (mprior < MINPRIORFROMHYPERPRIOR)
        IM_err (IMERR_HYPERPROB," Migration hyperprior value %.4lg  too low,  can't be less than %.4lg",mprior, (double) MINPRIORFROMHYPERPRIOR);
      hyperprior_uniform_m_max = m_max = mprior;
      hyperprior_expo_m_mean = -1.0; // not in use
    }
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    {
      if (thetaprior < MINPRIORFROMHYPERPRIOR)
        IM_err (IMERR_HYPERPROB," Population size hyperprior value %.4lg  too low,  can't be less than %.4lg",thetaprior, (double) MINPRIORFROMHYPERPRIOR);
      hyperprior_uniform_q_max = q_max = thetaprior;
    }
  }
  else
  {
    if (modeloptions[EXPOMIGRATIONPRIOR])
    {
      expo_m_mean = mprior;
      hyperprior_expo_m_mean = -1.0; // not in use
      hyperprior_uniform_m_max = -1.0; // i.e. not in use
    }
    else
    {
      m_max = mprior;
      hyperprior_expo_m_mean = -1.0; // not in use
      if (modeloptions[EXPOMIGRATIONPRIOR])
      {
        hyperprior_expo_m_mean = mprior;
        hyperprior_uniform_m_max = -1.0; // i.e. not in use
      }
      else
      {
        hyperprior_uniform_m_max = mprior;
        hyperprior_expo_m_mean = -1.0; // not in use
      }
    }
  }
  // set doRYupdate and doNWupdate
  doRYupdate = 0;
#ifdef DO_RY1UPDATE
      doRYupdate = 1;
#endif
    doNWupdate = 0;
#ifdef DO_NWUPDATE
    if (hiddenoptions[HIDDENGENEALOGY]==0 && modeloptions[NOMIGRATION] == 0)
      doNWupdate = 1;
#endif
  setup (infilename, fpstri, fpstr,priorfilename, topologypriorinfostring,currentid); // setup() is in initialize.cpp and does all the big structures like L, C, T
  add_previous_run_info_to_output(fpstri, fpstr);
  if (runoptions[LOADMCSTATE]) // set name of mcf file to read
  {
	  if (strcmp(&mcfreadfilename[strlen(mcfreadfilename)-4],".mcf")==0)
      sprintf(mcfreadfilename, "%s.%d", mcfreadfilename, currentid);
    else
      sprintf(mcfreadfilename, "%s.mcf.%d", mcfreadfilename, currentid);
    readmcf (mcfreadfilename, &mcmcrecords,&hilike,&hiprob, currentid);
    mcf_was_read = 1;
  }
  else
    mcf_was_read =0;

  if (runoptions[SAVELOADSAMEMCFFILE]) // check if mcf file is present and read if it is
  {
    sprintf(mcfwritefilename, "%s.mcf.%d", outfilename, currentid);
    sprintf(mcfreadfilename, "%s", mcfwritefilename);
    runoptions[SAVEMCSTATEFILE] = 1;
    mcffound = 0;
    if (currentid == HEADNODE)
      mcf0found = 0;
    if ((checkfile = fopen (mcfreadfilename, "r")) != NULL)
    {
      fclose (checkfile);
      sprintf(oldoutfilename, "%s.old", outfilename);
      //runoptions[LOADMCSTATE] = 1;
      readmcf (mcfreadfilename, &mcmcrecords,&hilike,&hiprob, currentid);
      mcf_was_read = 1;
      /*burndurationmode = TIMESTEPS;
      burnduration = (int) 0;  // no burn,  regardless of -b value
      if (runoptions[MCFLOADONLYSTATESPACE] == 1)  // ignore saved values in .mcf and .ti file, set burndone ==0 to overright .ti file with reset_after_burn
        burndone = 0;
      else
        burndone = 1;
      checkhighs(whichiscoldchain(),currentid,0);
      //setheat (hval1, hval2, heatmode, currentid); under some circumstances might want to reset heat terms after loading mcf files
      mcffound = 1;
      if (currentid == HEADNODE)
        mcf0found = 1; */
    }
  }
  if (mcf_was_read)
  {
    numpriormcfruns += 1;
    burndurationmode = TIMESTEPS;
    if (runoptions[MCFLOADONLYSTATESPACE] == 1 )
      burnduration = (int) 0;  // no burn,  regardless of -b value
    if (runoptions[MCFLOADONLYSTATESPACE] == 1 || runoptions[LOADMCSTATE])  // make new .ti file, set burndone ==0 to overright .ti file with reset_after_burn
      burndone = 0;
    else
      burndone = 1;
    //checkhighs(whichiscoldchain(),currentid,0); gets done anyway below
    //setheat (hval1, hval2, heatmode, currentid); under some circumstances might want to reset heat terms after loading mcf files
    mcffound = 1;
    if (currentid == HEADNODE)
      mcf0found = 1;
    add_previous_run_info_to_output(fpstri, fpstr);
    if (hiddenoptions[READOLDMCFFILE]==0)
    {
      if (runoptions[MCFLOADONLYSTATESPACE])
      {
        mcmcrecords_old = 0;
        genealogysamples_old = 0;
      }
      else
      {
        mcmcrecords_old = mcmcrecords;
        genealogysamples_old = genealogysamples;
      }
      mcmcrecords = 0;
      genealogysamples = 0;
      burninsteps_old = burninsteps;
      burninsteps = 0;
      runsteps_old = runsteps;
      runsteps = 0;
    }
    else
    {
      runsteps_old = runsteps = burninsteps_old = burninsteps = genealogysamples_old = 0;
      mcmcrecords_old = mcmcrecords;
      mcmcrecords = 0;
    }

#ifdef MPI_ENABLED
		  if (numprocesses > 1)
    {
      rc = MPI_Bcast(&mcf0found, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	   if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	   }
#endif
    if (mcffound != mcf0found)
      IM_err (IMERR_MCFREADFAIL, "not all cpus found an MCF file");
    if (mcf0found)
    {
      time(&lasttime); // set lasttime, same as if the burn had just finished
#ifdef MPI_ENABLED
		    if (numprocesses > 1)
      {
        rc = MPI_Bcast(&lasttime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	     if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	     }
#endif
    }
    //else  file does not exist,  do a regular run
  }
  else
    numpriormcfruns = 0;
  if (mcf_was_read == 0 || runoptions[MCFLOADONLYSTATESPACE] == 1)
    checkautoc (1, 0, 0, currentid);
  gsampinflength = calc_gsampinf_length ();
  checkhighs(whichiscoldchain(),currentid,1);
  // AS: TODO Looks like outputoptions has to be broadcast as well
  if (outputoptions[MIGRATEHIST])
  {
    sprintf(migplotfilename, "%s.mpt", outfilename);
    migcount = alloc2Dint_contiguous(&migcountdata, (nloci + (nloci > 1)),nummigdirs);
    migfrom = static_cast<int *> (malloc (nummigdirs * sizeof (int)));
    migto = static_cast<int *> (malloc (nummigdirs * sizeof (int)));
    for (i = 0; i < nummigdirs; i++)
    {
      migfrom[i]= atoi(strchr(migration_counts[0][i].str,'m')+1);
      migto[i] = atoi(strchr(migration_counts[0][i].str,'>')+1);
    }

  }
  if (runoptions[PRINTBURNTREND] && modeloptions[POPTREETOPOLOGYUPDATE])
    init_burn_phylogeny_counts();
#ifdef DEBUG
/*  code for reading in mcf file written by IMa2
   useful for debugging
  do
  {
    if (GSAMPINFOEXTRA)
    {
      //static char ima2mcffilename[21] = "ima2makemcf.out.mcf";
      static char ima2mcffilename[120] = "E:\\genemod\\ML-MCMC\\SEAI\\IMa2\\JH_branch_8_12_2014\\code\\test\\ima2makemcf.out.mcf";
      readima2mcf(ima2mcffilename);
    }

  }
  while (1); */
#endif
  /*
    8/25/2017  added mini-burn code,
    when topology is updated, it seems like burnin takes forever
    possibe that a lot of this is associated with getting genealogies and poptrees in sync
    so maybe it would help if just burned genealogies for awhile before starting all the other updates
    just made numminiburn numbers up
    this does seem to help
    makes sense because the genealogies start out as pretty arbitrary with respect to the model
    also played around a lot with extra hg updates during the run,  but this does not seem to help
  */
  if ( (runmode == POPTREEHYPERPRIORmode0 || runmode == POPTREEmode1) && runoptions[LOADRUN]==0 && (runoptions[LOADMCSTATE] == 0 || runoptions[SAVELOADSAMEMCFFILE]==1))
    if (runoptions[SAVELOADSAMEMCFFILE]==0 || (runoptions[SAVELOADSAMEMCFFILE]==1 && mcffound == 0))
    {
      int dominiburn = 1;
      int numminiburn;
      float meanss = total_numgenes/ (float) nloci;
      if (meanss < 10.0)
        numminiburn = 100;
      else if (meanss < 30.0)
        numminiburn = 1000;
      else
        numminiburn = 10000;
      if (dominiburn)
      {
        for (int ci = 0; ci < numchainspp; ci++)
        {
          for (int li = 0;li<nloci;li++)
          {
            assert (hiddenoptions[HIDDENGENEALOGY]==1);
            int topolchange;
            int tmrcachange;
            for (int bi = 0; bi < numminiburn;bi++)
              update_hidden_genealogy (ci, li, &topolchange, &tmrcachange);
          }
        }
        //if (currentid == HEADNODE)
        //  printf("mini-burn done\n");
      }
    }
}                               /* start */

/* 4/13/2017   JH redid qupdate with a different scheme for managing the intervals between updates */
/* interval is # steps _between_,  period is interval + 1,  minimum interval is 0  (every step) */
#define GENEALOGYUPDATEINTERVAL 0
#define POPTOPOLUPDATEINTERVAL 0
#define USCALARUPDATEINTERVAL  9
#define CHAINSWAPINTERVAL 0
#define POPRYUPDATEINTERVAL  1
#define POPNWUPDATEINTERVAL  1
#define POPSLIDEUPDATEINTERVAL 1
#define MIGPRIORUPDATEINTERVAL  9
#define QPRIORUPDATEINTERVAL  9
#define MIGPRIORNOTINTREEUPDATEINTERVAL 500
#define QPRIORNOTINTREEUPDATEINTERVAL 500
void qupdate (int currentid)
{
/* updatecounter is a simple stucture to determine when things get updated
  .interval is # steps between updates  (0 means update every step)
  .count is the current # of steps since last update
  By starting .count at a nonzero value,  can stagger updates that  have the same intervals
  can also set .interval+1 values to be prime numbers

  to turn off an update set .interval to INTMAX
  */
  struct updatecounter{
    int interval;
    int count;
  };
  int j, k, li, ci, ui;
  int changed;
  int qswapped = 0;
  int topolchange, tmrcachange;
  int periodpick;
  static struct updatecounter gene_uc,popslide_uc,poptopol_uc,popRY_uc,popNW_uc,uscale_uc,migprior_uc,chainswap_uc;
  static struct updatecounter migpriornit_uc,qprior_uc,qpriornit_uc;
  int mpa;
  int mpanitattempt,mpanitaccept;
  int poptopolchange,poptmrcachange;
  int trypoptopolchange,trypoptmrcachange;
  if (step==0)
  {
    gene_uc.count = 0;
    gene_uc.interval = GENEALOGYUPDATEINTERVAL;
    popslide_uc.count = 0;
    popslide_uc.interval = POPSLIDEUPDATEINTERVAL;
    popRY_uc.count = 0;
    popRY_uc.interval = POPRYUPDATEINTERVAL;
    popNW_uc.count = 1;
    popNW_uc.interval = POPNWUPDATEINTERVAL;
    poptopol_uc.count = 2;
    poptopol_uc.interval = POPTOPOLUPDATEINTERVAL;
    migprior_uc.count = 3;
    migprior_uc.interval = MIGPRIORUPDATEINTERVAL;
    uscale_uc.count = 1; // stagger these updates with respect to others
    uscale_uc.interval = USCALARUPDATEINTERVAL ;
    chainswap_uc.count = 0;
    chainswap_uc.interval = CHAINSWAPINTERVAL;
    migpriornit_uc.count = 0;
    migpriornit_uc.interval =   MIGPRIORNOTINTREEUPDATEINTERVAL;
    qprior_uc.count = 4; // stagger these updates with respect to others
    qprior_uc.interval  = QPRIORUPDATEINTERVAL;
    qpriornit_uc.count = 250; // stagger these updates with respect to others
    qpriornit_uc.interval = QPRIORNOTINTREEUPDATEINTERVAL;
}

#ifdef TURNONCHECKS
  checkgenealogy(0,0,0);
#endif //TURNONCHECKS
  int z = whichiscoldchain();  // if cold chain not on this cpu, then z=-1

  /* update genealogies */
  if (gene_uc.interval >= gene_uc.count)
  {
    for (ci = 0; ci < numchainspp; ci++)
    {
#ifdef TURNONCHECKS
checkgenealogyweights(ci);
  checkpoptree(ci,0);
  check_hgprob_sums(ci);
  checkprobs(ci,-1);
#endif //TURNONCHECKS
      for (li = 0;li<nloci;li++)
      {
#ifdef TURNONCHECKS
        checkprobs(ci,li);
#endif //TURNONCHECKS
        if (ci == z) // only record acceptance and tries from the cold chain
        {
          L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].tries+= 1;
          L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].tries+= 1;
          L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].tries+= 1;
        }
        if (hiddenoptions[HIDDENGENEALOGY]==1)
        {
          {
            if (update_hidden_genealogy (ci, li, &topolchange, &tmrcachange))
            {

              if (ci == z) // only record acceptance and tries from the cold chain
              {
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp++;
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp += (topolchange > 0);
                L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp += (tmrcachange > 0);
              }
              #ifdef TURNONCHECKS
              pcheck(ci,1);
              #endif

            }
          }

#ifdef TURNONCHECKS
        checkgenealogy(ci,li,1);
        checkprobs(ci,-1);
#endif //TURNONCHECKS
        }
        else
        {
#ifdef TURNONCHECKS
          checkgenealogy(0,0,0);
#endif //TURNONCHECKS

          if (updategenealogy (ci, li, &topolchange, &tmrcachange))
          {
            if (ci == z) // only record acceptance and tries from the cold chain
            {
              L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_ANY].accp++;
              L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TOPOLOGY].accp += (topolchange > 0);
              L[li].g_rec->upinf[IM_UPDATE_GENEALOGY_TMRCA].accp += (tmrcachange > 0);
            }
            #ifdef TURNONCHECKS
            pcheck(ci,2);
            #endif
          }
        }
      }

    }
    gene_uc.count = 0;
  }
  else
  {
      gene_uc.count += 1;
  } // geneealogy

  /*  poptree topology updating  */
  if (npops > 1  && modeloptions[POPTREETOPOLOGYUPDATE] == 1)
  {
    if (poptopol_uc.count >= poptopol_uc.interval )
    {
      for (ci = 0; ci < numchainspp; ci++)
      {
        changed = change_poptree (ci,&trypoptopolchange,&poptopolchange,&trypoptmrcachange,&poptmrcachange, 1);
        #ifdef TURNONCHECKS
        pcheck(ci,3);
        #endif
        if (ci==z) // only record acceptance and tries from the cold chain
        {
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].tries++;
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].tries += trypoptopolchange;
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].tries += trypoptmrcachange;
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].accp += changed;
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].accp += (poptopolchange > 0);
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].accp += (poptmrcachange > 0);
        }
      }
      poptopol_uc.count = 0;
    }
    else
    {
      poptopol_uc.count += 1;
    }
  }//poptopol

  /*  poptree population branch slide updating   */
  if (npops > 1  && modeloptions[POPTREETOPOLOGYUPDATE] == 0 && hiddenoptions[HIDDENGENEALOGY]==1)
  {
    if (popslide_uc.count >= popslide_uc.interval )
    {
      for (ci = 0; ci < numchainspp; ci++)
      {
        changed = change_poptree (ci,&trypoptopolchange,&poptopolchange,&trypoptmrcachange,&poptmrcachange, 0);
        #ifdef TURNONCHECKS
        pcheck(ci,4);
        #endif
        if (ci==z) // only record acceptance and tries from the cold chain
        {
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].tries++;
          //poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].tries += trypoptopolchange;  not attempted, only a slide
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].tries += trypoptmrcachange;
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].accp += changed;
          //poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].accp += (poptopolchange > 0);  not attempted, only a slide
          poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].accp += (poptmrcachange > 0);
        }

      }
      popslide_uc.count = 0;
    }
    else
    {
      popslide_uc.count += 1;
    }
  }//popslide

  /* poptree RY update */
  if (doRYupdate)
  {
    if (popRY_uc.count >= popRY_uc.interval)
    {
      for (ci = 0; ci < numchainspp; ci++)
      {
        periodpick = randposint (numsplittimes);
        if (hiddenoptions[HIDDENGENEALOGY]==1)
          changed = changet_RYhg (ci, periodpick);
        else
          changed = changet_RY1 (ci, periodpick);
        #ifdef TURNONCHECKS
        pcheck(ci,5);
        #endif
        if (ci == z) // only record acceptance and tries from the cold chain
        {
          T[periodpick].upinf[update_type_RY].tries++;
          if (changed)
            T[periodpick].upinf[update_type_RY].accp++;
        }
      }
      popRY_uc.count = 0;
    }
    else
    {
      popRY_uc.count += 1;
    }
  } //RY
  /* poptree NW update */
  if (doNWupdate)
  {
    if (popNW_uc.count ==popNW_uc.interval)
    {
      for (ci = 0; ci < numchainspp; ci++)
      {
        periodpick = randposint (numsplittimes);
        changed = changet_NW (ci, periodpick);
        #ifdef TURNONCHECKS
        pcheck(ci,6);
        #endif
        if (ci == z) // only record acceptance and tries from the cold chain
        {
          T[periodpick].upinf[update_type_NW].tries++;
          if (changed)
            T[periodpick].upinf[update_type_NW].accp++;
        }
      }
      popNW_uc.count = 0;
    }
    else
    {
      popNW_uc.count += 1;
    }
  } //NW

  /* mutation scalar updates*/
  if ( (uscale_uc.count >= uscale_uc.interval) && domutationscalarupdate)
  {
    if (nurates > 1)
    {
      for (ci = 0; ci < numchainspp; ci++)
      {
        for (j = 0; j < (nurates - (nurates == 2)); j++)
        {
          ui = changeu (ci, j, &k);
          #ifdef TURNONCHECKS
          pcheck(ci,7);
          #endif
          if (ci == z) // only record acceptance and tries from the cold chain

          {
            L[ul[j].l].u_rec[ul[j].u].upinf->tries++;
            L[ul[k].l].u_rec[ul[k].u].upinf->tries++;
            if (ui == 1)

            {
              L[ul[j].l].u_rec[ul[j].u].upinf->accp++;
              L[ul[k].l].u_rec[ul[k].u].upinf->accp++;
            }
          }
        }
      }
    }
    else
    {
      if (nloci == 1 && L[0].model == HKY)      /* if there is just one HKY locus kappa needs updating on its own */
        for (ci = 0; ci < numchainspp; ci++)
        {
          changekappa (ci);
          #ifdef TURNONCHECKS
          pcheck(ci,8);
          #endif
        }
    }
    uscale_uc.count = 0;
  }
  else
  {
    uscale_uc.count += 1;
  }//uscale

  /* hyperprior updates */
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR]) // individual migration prior updates
    {
      if (migprior_uc.count >= migprior_uc.interval)
      {
        for (ci = 0; ci < numchainspp; ci++)
        {
          for (j = 0; j < nummigrateparams; j++)
          {
            mpa = update_migration_prior_intree(ci,j);
            #ifdef TURNONCHECKS
            pcheck(ci,9);
            #endif
            if (ci == z)
            {
              mh[j].upinf[PRIOR_UPDATE].tries++;
              if (mpa == 1)
                mh[j].upinf[PRIOR_UPDATE].accp++;
            }
          }
        }
        migprior_uc.count = 0;
      }
      else
      {
        migprior_uc.count += 1;
      }
      //now do popsize terms
      if (qprior_uc.count >= qprior_uc.interval)
      {
        for (ci = 0; ci < numchainspp; ci++)
        {
          for (j = 0; j < numpopsizeparams; j++)
          {
            mpa = update_popsize_prior_intree(ci,j);
            #ifdef TURNONCHECKS
            pcheck(ci,10);
            #endif
            if (ci == z)
            {
              qh[j].upinf[PRIOR_UPDATE].tries++;
              if (mpa == 1)
                qh[j].upinf[PRIOR_UPDATE].accp++;
            }
          }
        }
        qprior_uc.count = 0;
      }
      else
      {
        qprior_uc.count += 1;
      }
    }
    if (modeloptions[POPTREETOPOLOGYUPDATE] == 1)  // priors not in the current tree
    {
      if (migpriornit_uc.count >= migpriornit_uc.interval)  // migration priors not in the current tree
      {
        for (ci = 0; ci < numchainspp; ci++)
        {
          if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
          {
            update_migration_prior_not_intree(ci, &mpanitattempt,&mpanitaccept);
            #ifdef TURNONCHECKS
            pcheck(ci,11);
            #endif
          }
          if (ci == z)
          {
            mhnit->upinf[PRIOR_UPDATE].tries += mpanitattempt;
            mhnit->upinf[PRIOR_UPDATE].accp += mpanitaccept;;
          }
        }
      }
      else
      {
        migpriornit_uc.count += 1;
      }
      //now do popsize terms not in the current tree
      if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
      {
        if (qpriornit_uc.count >= qpriornit_uc.interval)
        {
          for (ci = 0; ci < numchainspp; ci++)
          {
            update_popsize_prior_not_intree(ci, &mpanitattempt,&mpanitaccept);
            #ifdef TURNONCHECKS
            pcheck(ci,12);
            #endif
            if (ci == z)
            {
              qhnit->upinf[PRIOR_UPDATE].tries += mpanitattempt;
              qhnit->upinf[PRIOR_UPDATE].accp += mpanitaccept;;
            }
          }
        }
        else
        {
          qpriornit_uc.count += 1;
        }
      }
    }
  } //migprior

  /* chain swapping */
  if (numchainspp > 1)
  {
    if (chainswap_uc.count >= chainswap_uc.interval)
    {
      if (numprocesses >1 )
      {
#ifdef MPI_ENABLED
        // 7/7/2016  used to call both swapchains and swapchains_bwprocesses  but finally figured out this was unnecessary and just let swapchains_bwprocesses do both
        // thse num... vars are really just for debugging
        int numattemptwithin = 0;
        int numattemptbetween = 0;
        int numsuccesswithin = 0;
        int numsuccessbetween = 0;

        swapchains_bwprocesses(currentid, swaptries_per_cpu,  &numattemptwithin,&numattemptbetween,&numsuccesswithin,&numsuccessbetween);

#endif
	     }
      else // numprocesses == 1
      {
        qswapped = swapchains (1,swaptries_per_cpu); // qswapped never used

      }
      chainswap_uc.count = 0;
    }
    else
    {
      chainswap_uc.count += 1;
    }
  }   // chain swapping update


#ifdef TURNONCHECKS
  for (ci=0;ci<numchainspp;ci++)
  {
    checkgenealogyweights(ci);
    pcheck(ci,0);// just in case I missed something
  }
#endif
  return;
}                               /* qupdate */
#undef GENEALOGYUPDATEINTERVAL
#undef POPTOPOLUPDATEINTERVAL
#undef USCALARUPDATEINTERVAL
#undef CHAINSWAPINTERVAL
#undef POPRYUPDATEINTERVAL
#undef POPNWUPDATEINTERVAL
#undef POPSLIDEUPDATEINTERVAL
#undef MIGPRIORUPDATEINTERVAL
#undef QPRIORUPDATEINTERVAL
#undef MIGPRIORNOTINTREEUPDATEINTERVAL
#undef QPRIORNOTINTREEUPDATEINTERVAL
#undef HIDDENGENEALOGYUPDATEEXTRA

#undef BRANCHUPDATEFRAC

/**********  writes ti files and calls function that writes mcf files****************/

void reset_after_burn (int currentid)
{
  int li, ui, i, j;
  int rc;
  time_t timer;
  //time (&chainstarttime);
  checkhighs(whichiscoldchain(),currentid,1);
#ifdef MPI_ENABLED
		/*if (numprocesses > 1)
  {
    //rc = MPI_Bcast(&chainstarttime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
//printf("5 %d %llu\n",currentid,chainstarttime);
	 	 //if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	 } */
#endif

  burninsteps = step - 1;
  runsteps = 0;
  if (cdurationmode == TIMEINF)
  {
    time (&lasttime);
    time (&timer);
  }
#ifdef MPI_ENABLED
		if (numprocesses > 1)
  {
    rc = MPI_Bcast(&timer, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	 if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
    rc = MPI_Bcast(&lasttime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	 if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	 }
#endif
  for (i = 0; i < lastperiodnumber; i++)
    for (j = 0; j < T[i].num_uptypes; j++)
      T[i].upinf[j].accp = T[i].upinf[j].tries = 0;

  for (li = 0; li < nloci; li++)
  {

    for (ui = 0; ui < L[li].nlinked; ui++)
    {
      for (j = 0; j < L[li].u_rec[ui].num_uptypes; j++)
        L[li].u_rec[ui].upinf[j].accp = L[li].u_rec[ui].upinf[j].tries = 0;
      if (L[li].model == HKY)
        for (j = 0; j < L[li].kappa_rec->num_uptypes; j++)
          L[li].kappa_rec->upinf[j].accp = L[li].kappa_rec->upinf[j].tries =
            0;
      // A_rec not used as of sometime in 2010, A gets enough updates when updating genealogy
      if (L[li].model == STEPWISE || L[li].model == JOINT_IS_SW)
      {
        for (j = 0; j < L[li].A_rec[ui].num_uptypes; j++)
          L[li].A_rec[ui].upinf[j].accp = L[li].A_rec[ui].upinf[j].tries = 0;
      }
      for (j = 0; j < L[li].g_rec->num_uptypes; j++)
        L[li].g_rec->upinf[j].accp = L[li].g_rec->upinf[j].tries = 0;
    }
  }
  if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
  {
    totaltopolupdates = 0;
    chain0topolupdates = 0;
    chain0topolswaps = 0;
    poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].tries = 0;
    poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].tries = 0;
    poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].tries = 0;
    poptreeuinfo->upinf[IM_UPDATE_POPTREE_ANY].accp  = 0;
    poptreeuinfo->upinf[IM_UPDATE_POPTREE_TOPOLOGY].accp  = 0 ;
    poptreeuinfo->upinf[IM_UPDATE_POPTREE_TMRCA].accp = 0 ;
    if (runoptions[LOADMCSTATE] == 0)
    {
      for (i=0;i<numpoptopologies;i++)
        poptopologycounts[i] = 0;
    }
  }
  checkautoc (1, burndone, burninsteps, currentid);
  if (runoptions[DONTSAVEGENEALOGIES]==0)
  {
	   if (currentid == HEADNODE)
    if (runoptions[MCFLOADONLYSTATESPACE] && mcf_was_read) // write a .ti file,  possibly overwriting one of the same name (which is ok)
    {
      if ((genealogyinfosavefile = fopen (genealogyinfosavefilename, "w")) == NULL)
      {
        IM_err (IMERR_CREATEFILEFAIL, "Error creating file for holding genealogy information");
      }
      fprintf (genealogyinfosavefile,
               "-------------------------------------------\n\n");
      fprintf (genealogyinfosavefile, "Header for genealogy file:  %s\n\n",
               genealogyinfosavefilename);
      fprintf (genealogyinfosavefile,
               "-------------------------------------------\n\n");
      fprintf (genealogyinfosavefile, "%s\n", fpstr);
      fprintf (genealogyinfosavefile,
               "-------------------------------------------\n\n");
      fprintf (genealogyinfosavefile, "End of header for genealogy file:  %s\n\n",
               genealogyinfosavefilename);
      fprintf (genealogyinfosavefile,
               "-------------------------------------------\n\n");
      if (hiddenoptions[GSAMPINFOEXTRA]==1)
      {
        for (i=0;i<numpopsizeparams;i++)
          fprintf (genealogyinfosavefile,"q_c\t");
        for (i=0;i<numpopsizeparams;i++)
          fprintf (genealogyinfosavefile,"q_f\t");
        for (i=0;i<numpopsizeparams;i++)
          fprintf (genealogyinfosavefile,"h\t");
        for (i=0;i<nummigrateparams;i++)
          fprintf (genealogyinfosavefile,"m_c\t");
        for (i=0;i<nummigrateparams;i++)
          fprintf (genealogyinfosavefile,"m_f\t");
        for (i=0;i<numpopsizeparams;i++)
          fprintf (genealogyinfosavefile,"q_intg\t");
        for (i=0;i<nummigrateparams;i++)
          fprintf (genealogyinfosavefile,"m_intg\t");
        fprintf (genealogyinfosavefile,"P_data\t");
        fprintf (genealogyinfosavefile,"P_G\t");
        for (i=0;i<npops-1;i++)
          fprintf (genealogyinfosavefile,"t\t");
        fprintf (genealogyinfosavefile,"P_hm\t");
        fprintf (genealogyinfosavefile,"c_hm\t");
        if (modeloptions[POPTREETOPOLOGYUPDATE])
          fprintf (genealogyinfosavefile,"treenum\t");
        fprintf (genealogyinfosavefile,"\n");
      }
      fprintf (genealogyinfosavefile, "VALUESSTART\n");
      FCLOSE (genealogyinfosavefile);
	   }
  }
  if (runoptions[SAVEMCSTATEFILE])
  {
    //writemcf (mcfwritefilename,command_line,mcmcrecords,hilike_rec,hiprob_rec,currentid);   // can be useful to save after the burn in case something happens to the rest of the run
    writemcf (mcfwritefilename,command_line,mcmcrecords,mcmcrecords_old,genealogysamples_old,burninsteps_old,runsteps_old,hilike,hiprob,currentid);
  }

}                               /* reset_after_burn() */

/**********  writes burntrend file ****************/
void writeburntrendfile (int currentid)
{
  FILE *burntrendfile;
  char burntrendfilename[FNSIZE];
  char oldburntrendfilename[FNSIZE];
  assert (runoptions[PRINTBURNTREND]);
  static time_t burntrendtime;
  struct tm *burntrendtimeinfo;
  static char burntimeinfostring[80];
  double seconds;
  int rc;
  time (&burntrendtime);
#ifdef MPI_ENABLED
		if (numprocesses > 1)
  {
    rc = MPI_Bcast(&burntrendtime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	 if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	 }
#endif
  burntrendtimeinfo = localtime(&burntrendtime);
  strftime(burntimeinfostring,80,"%x - %X", burntrendtimeinfo);
  seconds = difftime (burntrendtime, starttime);

  if (currentid == HEADNODE)
  {
    printf("%s",outputbanner("printing burn trend file"));
    fflush(stdout);
    strcpy (burntrendfilename, outfilename);
    strcat (burntrendfilename, ".burntrend.out");
    sprintf(oldburntrendfilename, "%s.old", burntrendfilename); // set the name of the old burntrend file
    if (file_exists(oldburntrendfilename))
      remove (oldburntrendfilename);   // delete that file
    rename (burntrendfilename, oldburntrendfilename);  // now rename the current burtrend file to that old name
    if ((burntrendfile = fopen (burntrendfilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL,
              "Error opening burntrend output file for writing");
    }
    fprintf (burntrendfile,"Plots of Runtime Information and Parameter Trends during Burnin \n");
    fprintf (burntrendfile, "\nBurn trend file printed: %s\n",burntimeinfostring);
    fprintf (burntrendfile, "Command line string : %s \n", command_line);
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1)
    {
      char temppoptreestring[POPTREESTRINGLENGTHMAX];
      fprintf (burntrendfile, "Population Tree in Input File: %s\n", startpoptreestring);
      strcpy(temppoptreestring,startpoptreestring);
      rewrite(temppoptreestring);//  this will put the string in standard order
      if (strcmp(temppoptreestring,startpoptreestring) != 0)
        fprintf (burntrendfile, "  Input File Tree rewritten with standard ordering: %s\n",temppoptreestring);
      if (modeloptions[ADDGHOSTPOP])
      {
        strcpy(temppoptreestring,startpoptreestring);
        addoutgroup(temppoptreestring);
        rewrite(temppoptreestring);
        fprintf (burntrendfile,"    Input File Tree rewritten with standard ordering and ghost population outgroup: %s\n",temppoptreestring);
      }
    }

    fprintf (burntrendfile, "========================================\n");
/*#ifdef INDEVELOPMENT
    output_update_scalars(whichiscoldchain(),currentid);
#endif */
    intervaloutput (burntrendfile, currentid);
    if (modeloptions[POPTREETOPOLOGYUPDATE])
      outputburntopologycounts(burntrendfile,currentid);// add topology distribution to the end of the file,  needs to be called regardless of currentid

    fprintf (burntrendfile, "========================================\n\n");
    fprintf (burntrendfile, "Current Step #: %d \n\n", step);
    if (trendspot > 1 && currentid == HEADNODE)
      callasciitrend (burntrendfile,trenddoublepoint,trendspot);
    else if (trendspot <= 1 && currentid == HEADNODE)
    {
      fprintf(burntrendfile, "burn period too short to plot trends,  trend recording begins at step %d \n",burntrendstartdelay);
    }
    fprintf (burntrendfile, "\nTime Elapsed Since Start of Run : %s\n",timestring(seconds));
		  //fprintf (burntrendfile, "\nTime Elapsed Since Start of Run : %d hours, %d minutes, %d seconds \n",
    //  (int) seconds / (int) 3600,((int) seconds / (int) 60) - ((int) 60 * ((int) seconds / (int) 3600)),(int) seconds - (int) 60 *((int) seconds / (int) 60));
		  fprintf (burntrendfile, "\nEnd of Burn trend output file\n");
    fclose (burntrendfile);
  }
  else
  {
    //need to do some operations on all nodes so that MPI works
    //intervaloutput() and outputburntopologycounts()  make MPI_Reduce() calls,  and these must be made on all processors
/*#ifdef INDEVELOPMENT
    output_update_scalars(whichiscoldchain(),currentid);
#endif */
    intervaloutput (NULL, currentid);
    if (modeloptions[POPTREETOPOLOGYUPDATE])
      outputburntopologycounts(NULL,currentid);
  }
}                               // writeburntrendfile

/* run() determines the status of a run in terms of whether it is in the burnin phase or not
  and of whether the run should keep going.
  if the burnin period has just ended,  some work is done
  run()  also  checks to see if it is time to print an output file
  returns 0 (stop) or 1 (continue run)
*/
int run (int currentid)
{
  static int checkinterval = 0;
  static int printburntrendstep = 0, burnrecordi = 1;
  int tempburndone;
  int rc = 0; //AS: return code for MPI C bindings
  //char runfilename[] = "IMrun";
  //char burnfilename[] = "IMburn";
  time_t timer;
  if (burndone)
  {
    switch (cdurationmode)
    {
      case TIMESTEPS:
        {
          return (step < (chainduration + burninsteps)); ///AS: Is this correct?
          break;
	       }
      case TIMEINF:
      {
        if (checkinterval < CHECKINTERVALSTEPS)
        {
          checkinterval++;
          return (1);
        }
        else
        {
          checkinterval = 0;
          if (maxedoutgenealogysave)
	         {
		          return(0);
	         }
          time (&timer);
#ifdef MPI_ENABLED
		        if (numprocesses > 1)
          {
            rc = MPI_Bcast(&timer, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	         if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	         }
#endif
          if ((timer - lasttime) > chainduration)
          {
            time (&lasttime);
#ifdef MPI_ENABLED
		        if (numprocesses > 1)
          {
            rc = MPI_Bcast(&lasttime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	         if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	         }
#endif
            if (currentid == HEADNODE)
            {

              continuerun = checkrunfile(runfilename);
            }
#ifdef MPI_ENABLED
		          if (numprocesses > 1)
            {
              rc = MPI_Bcast(&continuerun, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	           if (rc !=MPI_SUCCESS) {
			             MPI_Abort(MPI_COMM_WORLD,-1);
              }
	           }
#endif
	           if (continuerun == 0)
            {
		            return(0);
	           }
            else
            {
		            if (currentid == HEADNODE)
              {
                if (runoptions[DONTSAVEGENEALOGIES]==0 && genealogysamples > 0)
                {
                  savegenealogyfile (genealogyinfosavefilename, genealogyinfosavefile, &lastgenealogysaved, gsampinflength);
                }
		            }
              printoutput (currentid,0);
#ifdef MPI_ENABLED
		            if (numprocesses > 1)
              {
			             MPI_Barrier(MPI_COMM_WORLD);// is this necessary ??
		            }
#endif

            }
          }
          return (1);
        }
        break;
      }
      default:
        return (0);
        break;
    }
  }
  else
  {
    tempburndone = 0;
    switch (burndurationmode)
    {
    case TIMESTEPS:
      tempburndone = (step > burnduration);
      break;
    case TIMEINF:
      if (checkinterval < CHECKINTERVALSTEPS)
      {
        checkinterval++;
      }
      else
      {
        checkinterval = 0;
        time (&timer);
        tempburndone = (timer - lasttime) > burnduration;
#ifdef MPI_ENABLED
		      if (numprocesses > 1)
        {
          rc = MPI_Bcast(&timer, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	       if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
          rc = MPI_Bcast(&tempburndone, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	       if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	       }
#endif
      }
      break;
    default:
      return (0);
      break;
    }
    // plot burn trend if called for
    if (runoptions[PRINTBURNTREND] && (tempburndone || step >= burntrendstartdelay))
    {
      if (burnrecordi == recordint)
      {
        trendrecord (-1, currentid);       // record values of parameters that are in mcmc
        burnrecordi = 1;
      }
      else
      {
        burnrecordi++;
      }
      if (tempburndone)
      {
        if (noburn_mcfload == 0)
          writeburntrendfile (currentid);
/*#ifdef MPI_ENABLED
		      if (numprocesses > 1)
        {
			       //MPI_Barrier(MPI_COMM_WORLD);// is this necessary ??
		      }
#endif */
        printburntrendstep = 1;
        /* now check to see if IMburn file is present */
        if (burndurationmode == TIMEINF)
        {
          if (currentid == HEADNODE)
          {

            tempburndone = !checkrunfile(burnfilename);
            if (tempburndone==0) // continue burn
            {
              printf("IMburn status: continue burn\n");
              fflush(stdout);

            }
          }

#ifdef MPI_ENABLED
		        if (numprocesses > 1)
          {
            rc = MPI_Bcast(&tempburndone, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	         if (rc !=MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,-1);
	         }
#endif
	         if (tempburndone == 0)
          {
		          time (&lasttime); // start the clock again, continue burn
#ifdef MPI_ENABLED
		          if (numprocesses > 1)
            {
              rc = MPI_Bcast(&lasttime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	           if (rc !=MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD,-1);
	           }
#endif
	         }

        }
      }
      else
      {
        printburntrendstep++;
      }
    }
    if (tempburndone)
    {
      burndone = 1;
#ifdef MPI_ENABLED
      MPI_Barrier(MPI_COMM_WORLD); // not sure if necessary, but there is some bug that crashes things sometimes at the end of the burn 5/5/2017
#endif
      reset_after_burn (currentid);
    }
    /*if (calcoptions[CALCMARGINALLIKELIHOOD])
      stepstone_get_Lmax(); */
    return (1);
  }
}                               /* run */

#define  TRENDLASTPT (TRENDDIM - 1)
void inctrend (int m, int t, struct value_record *v, double newval)
{
  int j;
  for (j = m; j < TRENDLASTPT; j++)
  {
    /* assert (v->trend[j + 1] != 0); This asserts sometimes. */
    v->trend[j] = v->trend[j + 1];
  }
  v->trend[t] = newval;
  //assert(v != 0);
}

void trend_reset (struct value_record *v, int nv)
{
  int i, j;
  for (i = 0; i < nv; i++)
    for (j = 0; j < TRENDDIM; j++)
      v[i].trend[j] = 0;
}                               // trend_reset

/* Using trendrecord()
This function records trend lines
It works on instances of struct value_record
the value_record must first be initialized (probably in initialize.c)
Code for a particular value_record or array of value_records
can be placed in trendrecord at two places:
  in the "if (burndone && reset == 0) " section
  and in the "if (recordinc == recordtrendinc)" section

explanation for how trendline data are accumulated:
 - movespot is the point at which values are deleted from the array,
 - each new value is added to the end of the array (i.e. at trendspot)
 - all values from one position past movespot up to the end of the array are moved down one position
 - this erases the value at movespot and makes room for a new value at the end.
 -each time the replacement point (movespot) reaches the end of the array the time
	period between additions doubles
values are added more slowly as the run proceeds.
- this is because the time period doubles when movespot reaches the end
- the values to the left of movespot have half the density (in time) of those to the right
*/
//AS: adding currentid, so I can start recording values correctly on the head node
//2/7/2014

void trendrecord (int loadarrayj, int currentid)
{
  static int recordtrendinc = 1, recordinc = 0, movespot =  TRENDLASTPT;
  static int init = 1, reset = 0;
  int j, k, li, ui;
  double probrec;
#ifdef MPI_ENABLED
	MPI_Status status;
#endif
  if (burndone && reset == 0)   // reset all trend-related values after burnin
  {
    init = 1;
    reset = 1;
    if (lpgpd_v->do_trend)
    {
      trend_reset (lpgpd_v, 1);
    }
    if (modeloptions[POPTREETOPOLOGYUPDATE]==1 && poptreeuinfo->v->do_trend == 1)
      trend_reset (poptreeuinfo->v, 1);

    for (li = 0; li < nloci; li++)
    {
      if (runoptions[LOADRUN] == 0)
      {
        if (L[li].g_rec->v->do_trend)
        {
          trend_reset (L[li].g_rec->v, 1);
        }
        for (ui = 0; ui < L[li].nlinked; ui++)
        {
          if (L[li].u_rec[ui].v->do_trend)
          {
            trend_reset (L[li].u_rec[ui].v, 1);
          }
        }
        if (L[li].model == HKY)
        {
          if (L[li].kappa_rec->v->do_trend)
          {
            trend_reset (L[li].kappa_rec->v, 1);
          }
        }

      }
    }
      for (k = 0; k < lastperiodnumber; k++)
        if (T[k].v->do_trend)
          trend_reset (T[k].v, 1);

/* ADD ADDITONAL trend_reset() calls here */
  }
  if (init == 1)
  {
    trendspot = 0;
    recordtrendinc = 1;
    recordinc = 0;
    movespot = TRENDLASTPT;
    init = 0;
  }
  recordinc++;
  if (recordinc == recordtrendinc)
  {
    if (runoptions[LOADRUN])
    {

      if (lpgpd_v->do_trend)
        inctrend (movespot, trendspot, lpgpd_v,
                  gsampinf[loadarrayj][gsamp_pdgp] +
                  gsampinf[loadarrayj][gsamp_probgp]);

    }
    else
    {
      if (lpgpd_v->do_trend)
      {
	      int z = whichiscoldchain();
	      if (z >= 0 && currentid == HEADNODE)
        {
          if (hiddenoptions[HIDDENGENEALOGY] == 0)
	          inctrend (movespot, trendspot, lpgpd_v,
                  C[z]->allpcalc.probg + C[z]->allpcalc.pdg);
          else
            inctrend (movespot, trendspot, lpgpd_v,
                  C[z]->allpcalc.probg + C[z]->allpcalc.pdg + C[z]->allpcalc.probhgg);
	      }
	#ifdef MPI_ENABLED
	      if (z < 0 && currentid == HEADNODE)
        {
		      probrec = 0.0;
		      int rc = MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1313, MPI_COMM_WORLD, &status);
		      if (rc != MPI_SUCCESS)
          {
			      MPI_Abort(MPI_COMM_WORLD, rc);
		      }
		      inctrend (movespot, trendspot, lpgpd_v, probrec);
	      }
	      if (z >=0 && currentid != 0)
       {
          if (hiddenoptions[HIDDENGENEALOGY]==0)
		          probrec = C[z]->allpcalc.probg + C[z]->allpcalc.pdg;
          else
            probrec = C[z]->allpcalc.probg + C[z]->allpcalc.pdg + C[z]->allpcalc.probhgg;
		      int rc = MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 1313, MPI_COMM_WORLD);
		      if (rc != MPI_SUCCESS)
          {
			      MPI_Abort(MPI_COMM_WORLD, rc);
		      }
	      }
	#endif
	    }
      if (modeloptions[POPTREETOPOLOGYUPDATE]==1 && poptreeuinfo->v->do_trend)
      {
	      int z = whichiscoldchain();
	      if (z >= 0 && currentid == HEADNODE)
        {
          if (npops - modeloptions[ADDGHOSTPOP] < MINPOPSFORDISTANCE) // for fewer than 6 populations record the pop id #, else record the distance to tree 0
	           inctrend (movespot, trendspot, poptreeuinfo->v, (double) C[z]->poptreenum);
          else
            inctrend (movespot, trendspot, poptreeuinfo->v,(double) RFtreedis[C[z]->poptreenum]);
	      }
	#ifdef MPI_ENABLED
	      if (z < 0 && currentid == HEADNODE)
       {
		      probrec = 0.0;
		      int rc = MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 13131, MPI_COMM_WORLD, &status);
		      if (rc != MPI_SUCCESS)
          {
			      MPI_Abort(MPI_COMM_WORLD, rc);
		      }
		      inctrend (movespot, trendspot, poptreeuinfo->v, probrec);
	      }
	      if (z >=0 && currentid != 0)
       {
          if (npops - modeloptions[ADDGHOSTPOP] < MINPOPSFORDISTANCE)
		          probrec = (double) C[z]->poptreenum;
          else
            probrec = (double) RFtreedis[C[z]->poptreenum];
		        int rc = MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 13131, MPI_COMM_WORLD);
		        if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	      }
	#endif
	    }
     for (j = 0; j < lastperiodnumber; j++)
     {
        if (T[j].v->do_trend)
        {
		      int z = whichiscoldchain();
		      if (z >= 0 && currentid == HEADNODE)
        {
	         inctrend (movespot, trendspot, T[j].v, C[z]->tvals[j]);
		      }
#ifdef MPI_ENABLED
		      if (z < 0 && currentid == HEADNODE) {
			      probrec = 0.0;
			      int rc =  MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1717, MPI_COMM_WORLD, &status);
			      if (rc != MPI_SUCCESS)
         {
				      MPI_Abort(MPI_COMM_WORLD, rc);
			      }
			      inctrend (movespot, trendspot, T[j].v, probrec);
		      }
		      if (z >= 0 && currentid != 0)
        {
			      probrec = C[z]->tvals[j];
			      int rc =  MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 1717, MPI_COMM_WORLD);
			      if (rc != MPI_SUCCESS)
         {
				        MPI_Abort(MPI_COMM_WORLD, rc);
			      }
		      }
#endif
		    }
	     }
     }
    for (li = 0; li < nloci; li++)
    {
      if (runoptions[LOADRUN] == 0)
      {
        for (ui = 0; ui < L[li].nlinked; ui++)if (L[li].u_rec[ui].v->do_trend)
        {
		        int z = whichiscoldchain();
		        if (z >= 0 && currentid == HEADNODE)
          {
            inctrend (movespot, trendspot, L[li].u_rec[ui].v, C[z]->G[li].uvals[ui]);
		        }
#ifdef MPI_ENABLED
		        if (z < 0 && currentid == HEADNODE)
          {
			         probrec = 0.0;
			         int rc = MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1414, MPI_COMM_WORLD, &status);
			         if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);
			         inctrend (movespot, trendspot, L[li].u_rec[ui].v, probrec);
		        }
		        if (z >= 0 && currentid != 0)
          {
			         probrec = C[z]->G[li].uvals[ui];
			         int rc = MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 1414, MPI_COMM_WORLD);
			         if (rc != MPI_SUCCESS)		MPI_Abort(MPI_COMM_WORLD, rc);
		        }
#endif
		      }
	     }
      if (L[li].model == HKY)
      {
		      int z = whichiscoldchain();
		      if (z >= 0 && currentid == HEADNODE)
        {
	         inctrend (movespot, trendspot, L[li].kappa_rec->v, C[z]->G[li].kappaval);
		      }
#ifdef MPI_ENABLED
		      if (z < 0 && currentid == HEADNODE)
        {
			       probrec = 0.0;
			       int rc = MPI_Recv(&probrec, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1515, MPI_COMM_WORLD, &status);
			       if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
			       inctrend (movespot, trendspot, L[li].kappa_rec->v, probrec);
		      }
		      if (z >=0 && currentid != 0)
        {
			       probrec = C[z]->G[li].kappaval;
			       int rc = MPI_Send(&probrec, 1, MPI_DOUBLE, 0, 1515, MPI_COMM_WORLD);
			       if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
		      }
#endif
	     }
    }
/* ADD ADDITONAL inctrend() calls here */
    if (movespot == TRENDLASTPT && trendspot == TRENDLASTPT)
    {
      movespot = 0;
      recordtrendinc *= 2;
    }
    else
    {
      movespot += (movespot < TRENDLASTPT);
    }
    trendspot += (trendspot < TRENDLASTPT);
    recordinc = 0;
  }
  trenddoublepoint = movespot;
}                               /* trendrecord */

/* calculates the bin number of an xy array of a value_record that a value falls in,  increments the count in that bin */
//AS: I need to change this for MPI version - every time I record something on the cold chain process,
//it also has to be updated on the head node, so printing can happen correctly

void recordval (struct value_record *v, double val)
{

  int k;
  double logval;
  if (v->do_logplot)
  {
    assert (!(val < 0.0));
    logval = log (val);
    k =
      (int) (GRIDSIZE * (logval + v->plotrange.max) /
             (2.0 * v->plotrange.max));
  }
  else
  {
    k = (int) (GRIDSIZE * val / v->plotrange.max);
  }

  if (k < 0)
  {
    v->beforemin++;
  }
  else if (k >= GRIDSIZE)
  {
    v->aftermax++;
  }
  else
  {
    assert (!(k < 0));          // FIXME: it's been crashing
    v->xy[k].y++;
  }
  return;
}  /* recordval */

/**********  writes migrationnamefile ****************/
/* this is used to record the names of actual loci and gene copies that migrate, from to to
write these names to a file */
void record_migration_names(void)
 {
  int i, j, li;
  int from, to;
  struct edge *gtree;
  from = migrationnamefrom;
  to = migrationnameto;
  int z = whichiscoldchain();
  if (z >= 0) {
  for (li = 0; li < nloci; li++)
  {
    gtree = C[z]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > 0)
      {
        if (from == nowedgepop (0, &gtree[i], gtree[i].mig[j].mt) && to == C[z]->G[li].gtree[i].mig[j].mp)
        {
          fprintf(migrationnamefile,"%s ",L[li].name);
          if (i< L[li].numgenes)
            fprintf(migrationnamefile,"%s ",L[li].gNames[i]);
          else
            fprintf(migrationnamefile,"internal ");
        }
        j++;
      }
    }
  }
  fprintf(migrationnamefile,"\n");
  } else {
	return;
 }
 } // record_migration_names

/* INSTRUCTIONS to record a numerical value from the markov chain:
----------------------------------------------------
this works on instances of struct value_record

the value_record must be initiatlized (e.g. in initialize.c,  see e.g. init_g_rec)
this includes a call to init_value_record()

insert line(s) code into record()  below,  to make a call to recordval()

*/

void record_migrations (int z)
{
  int i, j, k, li, from, to, foundparam;
  struct edge *gtree;
	if (z >= 0) {

  for (j = 0; j < nloci + (nloci > 1); j++)
    for (i = 0; i < nummigdirs; i++)
    {
      migcount[j][i] = 0;
    }
  for (li = 0; li < nloci; li++)
  {
    gtree = C[z]->G[li].gtree;
    for (i = 0; i < L[li].numlines; i++)
    {
      j = 0;
      while (gtree[i].mig[j].mt > 0)
      {
        from = nowedgepop (z, &gtree[i], gtree[i].mig[j].mt);
        to = C[z]->G[li].gtree[i].mig[j].mp;
        k = 0;
        foundparam = 0;
        while (k < nummigdirs && foundparam == 0)
        {
          if (from == migfrom[k] && to == migto[k])
            foundparam = 1;
          else
            k++;
        }
        assert(k<nummigdirs);
        if (nloci == 1)
        {
          migcount[0][k]++;
        }
        else
        {
          migcount[li + 1][k]++;
          migcount[0][k]++;
        }
        j++;
      }
    }
  }
 /*for (i = 0; i < nloci + (nloci > 1); i++)
    for (k = 0; k < nummigdirs; k++)
    {
      recordval (&migration_counts [i][k], (double) migcount[i][k]);
    }  */
  }
}                               //record_migrations

void checkhighs (int z,int currentid, int reset)
{
  int rc = 0; //AS: return code for MPI C bindings
  double lowprob = -1e100;
  if (z>= 0)
  {
    currlike = C[z]->allpcalc.pdg;
    if (hilike < currlike || reset == 1)
      hilike = currlike;
    if (hiddenoptions[HIDDENGENEALOGY]==0)
      currprob = C[z]->allpcalc.probg;
    else
      currprob = C[z]->allpcalc.probg + C[z]->allpcalc.probhgg;
    if (hiprob < currprob || reset == 1)
        hiprob = currprob;
  }
  else if (reset)
  {
    hilike = currlike = hiprob = currprob = lowprob;
  }

#ifdef MPI_ENABLED
  MPI_Status status;
  double temp[4];
  if (z >=0 && currentid !=HEADNODE) // then send roottime to zero
  {
    temp[0] = hilike;
    temp[1] = currlike;
    temp[2] = hiprob;
    temp[3] = currprob;
    rc = MPI_Send(temp, 4, MPI_DOUBLE, 0, 12321, MPI_COMM_WORLD);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
  }
	 if (z < 0 && currentid == HEADNODE) // receive roottime
  {
		  rc = MPI_Recv(temp, 4, MPI_DOUBLE, MPI_ANY_SOURCE, 12321, MPI_COMM_WORLD, &status);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
    hilike = temp[0];
    currlike = temp[1];
    hiprob = temp[2];
    currprob = temp[3];
  }

  /*if (z >=0 && currentid !=HEADNODE) // then send roottime to zero
  {
    rc = MPI_Send(&hilike, 1, MPI_DOUBLE, 0, 232321, MPI_COMM_WORLD);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		  rc = MPI_Send(&currlike, 1, MPI_DOUBLE, 0, 232341, MPI_COMM_WORLD);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		  rc = MPI_Send(&hiprob, 1, MPI_DOUBLE, 0, 232361, MPI_COMM_WORLD);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		  rc = MPI_Send(&currprob, 1, MPI_DOUBLE, 0, 232381, MPI_COMM_WORLD);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	 }
	 if (z < 0 && currentid == HEADNODE) // receive roottime
  {
		  rc = MPI_Recv(&hilike, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 232321, MPI_COMM_WORLD, &status);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		  rc = MPI_Recv(&currlike, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 232341, MPI_COMM_WORLD, &status);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		  rc = MPI_Recv(&hiprob, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 232361, MPI_COMM_WORLD, &status);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		  rc = MPI_Recv(&currprob, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 232381, MPI_COMM_WORLD, &status);
		  if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	 }
  MPI_Barrier(MPI_COMM_WORLD); // should be a way to do this better that avoids a barrier
  */
#endif
  //printf("currlike: %.2lf hilike: %.2lf currprob: %.2lf hiprob: %.2lf\n",currlike,hilike,currprob,hiprob);
}                               /* checkhighs */

void record (int currentid)
{
  int j, li, ui;
  struct genealogy *G;
  int z = whichiscoldchain();  // if the current cpu hosts the cold chain,  z will be >= 0,  else it will be -1
  int rc = 0; //AS: return code for MPI C bindings
  int mj;
  int ptn;
  double mval;
#ifdef MPI_ENABLED
	MPI_Status status;
#endif

 /*
    need to get chain 0 stuff to head node
    four possibilities:
      1. chain 0 on current node and current node is headnode  (simply record the info)
      2. chain 0 on current node and current node is NOT headnode (send chain 0 info to headnode)
      3. chain 0 NOT on current node and current node is the headnode  (receive chain 0 info and record the info)
      4. chain 0 NOT on current node and current node is NOT the headnode (do nothing)
 */

 /*
  current and highest likelihoods and genealogy priors
 */
 checkhighs(z,currentid,0);

  if (outputoptions[PRINTTMRCA])
  {
    for (li = 0; li < nloci; li++)
    {

	    if (z >= 0 && currentid == HEADNODE)
     {
	      recordval (L[li].g_rec->v, C[z]->G[li].roottime);
	    }
#ifdef MPI_ENABLED
	    if (z >=0 && currentid !=0) // then send roottime to zero
     {
       double rtime = C[z]->G[li].roottime;
		     rc = MPI_Send(&rtime, 1, MPI_DOUBLE, 0, 2323, MPI_COMM_WORLD);
		     if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	    }
	    if (z < 0 && currentid == HEADNODE) // receive roottime
     {
		     double rtime = 0.0;
		     rc = MPI_Recv(&rtime, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2323, MPI_COMM_WORLD, &status);
		     if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
		     recordval (L[li].g_rec->v, rtime);
	    }
#endif
    }
  }

  if (domutationscalarupdate) for (li = 0; li < nloci; li++)
  {
    if (z >= 0)
    {
    G = &(C[z]->G[li]);
    }
    for (ui = 0; ui < L[li].nlinked; ui++)
    {
	    if (z >= 0 && currentid == HEADNODE)
     {
	      recordval (L[li].u_rec[ui].v, G->uvals[ui]);
	    }
#ifdef MPI_ENABLED
	    if (z >=0 && currentid !=0)
      {
		    double uvals = G->uvals[ui];
		    rc = MPI_Send(&uvals, 1, MPI_DOUBLE, 0, 2424, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)
			    MPI_Abort(MPI_COMM_WORLD, rc);
	    }
	    if (z < 0 && currentid == HEADNODE)
      {
		    double uvals = 0.0;
		    rc = MPI_Recv(&uvals, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2424, MPI_COMM_WORLD, &status);
		    if (rc != MPI_SUCCESS)
			    MPI_Abort(MPI_COMM_WORLD, rc);
		    recordval (L[li].u_rec[ui].v, uvals);
	    }
#endif
    }
    if (L[li].model == HKY)
    {
	    if (z >= 0 && currentid == HEADNODE)
      {
	          recordval (L[li].kappa_rec->v, G->kappaval);
	    }
#ifdef MPI_ENABLED
	    if (numprocesses > 1)
      {
		    if (z >= 0 && currentid != 0)
        {
			    double kval = G->kappaval;
			    rc = MPI_Send(&kval, 1, MPI_DOUBLE, 0, 2345, MPI_COMM_WORLD);
			    if (rc != MPI_SUCCESS)
				    MPI_Abort(MPI_COMM_WORLD, rc);
		    }
		    if (z < 0 && currentid == HEADNODE)
        {
			    double kval = 0.0;
			    rc = MPI_Recv(&kval, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2345, MPI_COMM_WORLD, &status);
			    if (rc != MPI_SUCCESS)
				    MPI_Abort(MPI_COMM_WORLD, rc);
			    recordval (L[li].kappa_rec->v, kval);
		    }
	    }
#endif
    }
  }
  for (j = 0; j < lastperiodnumber; j++)
  {
	  if ( z >= 0)
    {
	    assert (C[z]->tvals[j] > T[j].pr.min && C[z]->tvals[j] < T[j].pr.max); // causes a bug when using a priorfile with -c4 thermodynamic integration
	  }
	  if (z >= 0 && currentid == HEADNODE)
    {
	        recordval (T[j].v, C[z]->tvals[j]);
	  }
#ifdef MPI_ENABLED
	  if (z >= 0 && currentid != 0)
    {
		  double tvals = C[z]->tvals[j];
		  rc = MPI_Send(&tvals, 1, MPI_DOUBLE, 0, 2525, MPI_COMM_WORLD);
		  if (rc != MPI_SUCCESS)
			MPI_Abort(MPI_COMM_WORLD, rc);
	  }
	  if (z < 0 && currentid == HEADNODE)
    {
		  double tvals = 0.0;
		  rc = MPI_Recv(&tvals, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 2525, MPI_COMM_WORLD, &status);
		  if (rc != MPI_SUCCESS)
			MPI_Abort(MPI_COMM_WORLD, rc);
		  recordval (T[j].v, tvals);
	  }
#endif
  }
  if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
  {
    for (j = 0; j < nummigrateparams; j++)
    {
      mj = j;
	     if ( z >= 0)
      {
        if (modeloptions[EXPOMIGRATIONPRIOR])
        {
          mval = C[z]->imig[mj].pr.expomean;
          assert (mval < hyperprior_expo_m_mean * EXPOMIGPLOTSCALE && mval > 0.0);
        }
        else
        {
          mval = C[z]->imig[mj].pr.max;
          assert (mval < hyperprior_uniform_m_max && mval > 0.0);
        }
	     }
	     if (z >= 0 && currentid == HEADNODE)
      {
	           recordval (mh[j].v, mval);
	     }
#ifdef MPI_ENABLED
	     if (z >= 0 && currentid != 0)
      {
		    rc = MPI_Send(&mval, 1, MPI_DOUBLE, 0, 25251, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)
			   MPI_Abort(MPI_COMM_WORLD, rc);
	     }
	     if (z < 0 && currentid == HEADNODE)
      {
		    mval = 0.0;
		    rc = MPI_Recv(&mval, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 25251, MPI_COMM_WORLD, &status);
		    if (rc != MPI_SUCCESS)
			   MPI_Abort(MPI_COMM_WORLD, rc);
		    recordval (mh[j].v, mval);
	     }
#endif
    }
    if (modeloptions[POPSIZEANDMIGRATEHYPERPRIOR])
    {
      for (j = 0; j < numpopsizeparams; j++)
      {
	       if ( z >= 0)
        {
          mval = C[z]->itheta[j].pr.max;
	       }
	       if (z >= 0 && currentid == HEADNODE)
        {
         recordval (qh[j].v, mval);
	       }
#ifdef MPI_ENABLED
	       if (z >= 0 && currentid != 0)
        {
		      rc = MPI_Send(&mval, 1, MPI_DOUBLE, 0, 25251, MPI_COMM_WORLD);
		      if (rc != MPI_SUCCESS)
			     MPI_Abort(MPI_COMM_WORLD, rc);
	       }
	       if (z < 0 && currentid == HEADNODE)
        {
		      mval = 0.0;
		      rc = MPI_Recv(&mval, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 25251, MPI_COMM_WORLD, &status);
		      if (rc != MPI_SUCCESS)
			     MPI_Abort(MPI_COMM_WORLD, rc);
		      recordval (qh[j].v, mval);
	       }
#endif
      }
    }
  }
/*
  for (i = 0; i < nloci + (nloci > 1); i++)
    for (k = 0; k < nummigdirs; k++)
    {
      recordval (&migration_counts [i][k], (double) migcount[i][k]);
    }
 */

  if (outputoptions[MIGRATEHIST])
  {
    int i,k;
    int rc = (nloci + (nloci > 1)) * nummigdirs;
    if (z >= 0)
      record_migrations (z);  // puts counts into migcount

    if (z>= 0  && currentid == HEADNODE)
    {
      for (i = 0; i < nloci + (nloci > 1); i++)
        for (k = 0; k < nummigdirs; k++)
        {
          recordval (&migration_counts [i][k], (double) migcount[i][k]);
        }
     }
#ifdef MPI_ENABLED
    if (z >= 0 && currentid != HEADNODE)
    {
      /*can do MPI_Send and MPI_Recv of 2D array because it was malloced
        as a contiguous block using alloc2Dint_contiguous() */
      rc = MPI_Send(&(migcount[0][0]), rc, MPI_INT, 0, 97525, MPI_COMM_WORLD); // must make sure the address is the location where the integers begin,  not the pointer
      if (rc != MPI_SUCCESS)
      MPI_Abort(MPI_COMM_WORLD, rc);
    }
    if (z < 0 && currentid == HEADNODE)
    {
    double tvals = 0.0;
    rc = MPI_Recv(&(migcount[0][0]), rc, MPI_INT, MPI_ANY_SOURCE, 97525, MPI_COMM_WORLD, &status); // must make sure the address is the location where the integers begin,  not the pointer
    if (rc != MPI_SUCCESS)
    MPI_Abort(MPI_COMM_WORLD, rc);
    for (i = 0; i < nloci + (nloci > 1); i++)
      for (k = 0; k < nummigdirs; k++)
      {
        recordval (&migration_counts [i][k], (double) migcount[i][k]);
      }
    }
#endif
  }

  if (npops >= 3 && npops <= 5 && outputoptions[PRINTJOINTTEST])
      setup_multi_t_arrays (z);
  if (!outputoptions[DONTPRINTASCIITREND])
  {
    trendrecord (-1, currentid);
  }
  if (calcoptions[CALCMARGINALLIKELIHOOD])
  {
    summarginlikecalc();
  }
  if (modeloptions[POPTREETOPOLOGYUPDATE] == 1 && currentid == HEADNODE)
    poptopologiessampled +=1;
  if (modeloptions[POPTREETOPOLOGYUPDATE] == 1  && z >= 0) // should record one tree for every time the CPU with the cold chain comes thru here,
  {
    poptopologycounts[C[z]->poptreenum] += 1;
    assert (beta[z] == 1.0);
  }
  if (modeloptions[POPTREETOPOLOGYUPDATE] == 1) // add the current topology number from the cold chain to poptopologysequence
  {
	    if (z >= 0 && currentid == HEADNODE)
     {
       poptopologysequence.vals[poptopologysequence.currentlength] = C[z]->poptreenum;
       poptopologysequence.disvals[poptopologysequence.currentlength] =  (double) RFtreedis[C[z]->poptreenum];
       poptopologysequence.currentlength += 1;
     }
#ifdef MPI_ENABLED
	    if (z >=0 && currentid !=0) // then send roottime to zero
     {
       ptn =C[z]->poptreenum;
       assert (ptn >= 0 && ptn < numtreesarray[npops - modeloptions[ADDGHOSTPOP]]);
		     rc = MPI_Send(&ptn, 1, MPI_INT, 0, 179, MPI_COMM_WORLD);
		     if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
	    }
	    if (z < 0 && currentid == HEADNODE) // receive roottime
     {
		     ptn = 0;
		     rc = MPI_Recv(&ptn, 1, MPI_INT, MPI_ANY_SOURCE, 179, MPI_COMM_WORLD, &status);
		     if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
       assert (ptn >= 0 && ptn < numtreesarray[npops - modeloptions[ADDGHOSTPOP]]);
       poptopologysequence.vals[poptopologysequence.currentlength] = ptn;
       poptopologysequence.disvals[poptopologysequence.currentlength] = (double) RFtreedis[ptn];
       assert(poptopologysequence.vals[poptopologysequence.currentlength] >= 0 && poptopologysequence.vals[poptopologysequence.currentlength] <  numtreesarray[npops - modeloptions[ADDGHOSTPOP]]);
       poptopologysequence.currentlength += 1;
	    }
#endif
     if (currentid == HEADNODE)
     {
       if ((poptopologysequence.maxlength - poptopologysequence.currentlength) < 1)
       {
         poptopologysequence.vals =  static_cast<int *>  (realloc (poptopologysequence.vals,(poptopologysequence.maxlength + POPTOPOLOGYSEQUENCELENGTH) * sizeof(int)));
         poptopologysequence.disvals =  static_cast<double *>  (realloc (poptopologysequence.disvals,(poptopologysequence.maxlength + POPTOPOLOGYSEQUENCELENGTH) * sizeof(double)));
         poptopologysequence.maxlength += POPTOPOLOGYSEQUENCELENGTH;
       }
     }
  }

  return;
}                               /* record */

/*
  whichiscoldchain()  examines the beta values to see which one is 1.0
  and returns the corresponding chain number.

  Important:  if there are multiple cpus,  then the current cpu may not have a chain with beta == 1.0
  in this case this function returns -1.

  This function is called alot, as the return value can be used to see, not only does the current cpu have the cold chain,
  but also, if it does,  which chain # it is.
*/

int
whichiscoldchain (void)
{
	int which = -1;
	for (int i = 0; i < numchainspp; i++) {
		if (beta[i] == 1.0) {
			which = i;
		}
	}
	//else return a -ve flag
	return which;
}

void savegenealogyinfo (int currentid)        // use floats to save space
{
  int i;
#ifdef MPI_ENABLED
  MPI_Status status;
  MPI_Request request[1];
  int rc = 0; //Return code for C Bindings
#endif
  float *gsampinflocal;
   // debugging genealogies when sampling trees
  char debug_ti_addinfolocalstr[1000];
  char *debug_ti_addinfolocal = debug_ti_addinfolocalstr;

//  request = new MPI_Request[1];

///Allocate memory for the full gsampinf only if I am on the head node
//else, I only need a local structure to be filled up, then I'll send that to the head node
  if (genealogysamples == 0)
  {
    if (cdurationmode == TIMESTEPS)
    {
	///AS: This is a problem - each process is going to now have its own set of genealogies
	//that are going to be saved, and the total number of genealogies should be = genealogiestosave
	//so we need to dynamically allocate this. This current method is wasteful.
	//Better way to do is to create a 2D vector gsampinf, keep push_back() into it dynamically
	//But keeping this for now as on Mon Mar 31 16:14:33 EDT 2014
      gsampinf = static_cast<float **> (malloc (genealogiestosave * sizeof (float *)));
      for (i = 0; i < genealogiestosave; i++)
        gsampinf[i] = static_cast<float *> (malloc (gsampinflength * sizeof (float)));
      memforgenealogiessaved = genealogiestosave;

      if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)  // debugging genealogies when sampling trees
      {
        debug_ti_addinfo = static_cast<char **> (malloc (genealogiestosave  * sizeof (char *)));
        for (i = 0; i < genealogiestosave ; i++)
        {
            debug_ti_addinfo[i] = static_cast<char *> (malloc (1000 * sizeof (char)));
        }
      }


    }
    else                        /* cdurationmode == TIMEINF */
    {
      gsampinf = static_cast<float **> (malloc (MAXGENEALOGIESTOSAVE * sizeof (float *)));
      for (i = 0; i < MAXGENEALOGIESTOSAVE; i++)
        gsampinf[i] = static_cast<float *> (malloc (gsampinflength * sizeof (float)));
      memforgenealogiessaved = MAXGENEALOGIESTOSAVE;
      if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)  // debugging genealogies when sampling trees
      {
        debug_ti_addinfo = static_cast<char **> (malloc (MAXGENEALOGIESTOSAVE  * sizeof (char *)));
        for (i = 0; i < MAXGENEALOGIESTOSAVE ; i++)
        {
            debug_ti_addinfo[i] = static_cast<char *> (malloc (1000 * sizeof (char)));
        }
      }

    }
  }
///AS: creating local gsampinf float array instead of length gsampinflength
//AS: ideally should do a check for cduration mode here as well and allocate gsampinflocal accordingly
	if (currentid != 0) {
		gsampinflocal = static_cast<float *> (malloc (gsampinflength * sizeof (float)));
	}

  if (genealogysamples >= (MAXGENEALOGIESTOSAVE - 1) && cdurationmode == TIMEINF)
  {
    printf (" maximum possible genealogies saved \n");
    fflush(stdout);
    maxedoutgenealogysave = 1;
    //bcast bug .  should replace currentid with 0  ??
#ifdef MPI_ENABLED
		//rc = MPI_Bcast(&maxedoutgenealogysave, 1, MPI_INT, currentid, MPI_COMM_WORLD);
    rc = MPI_Bcast(&maxedoutgenealogysave, 1, MPI_INT, 0, MPI_COMM_WORLD);
		  if (rc != MPI_SUCCESS)	MPI_Abort(MPI_COMM_WORLD, rc);
#endif

  }
  else
  {
	int z = whichiscoldchain();
	if (currentid == HEADNODE) {
		if ( z >= 0) {
      if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)
        savegsampinf_debug_ti_addinfo (gsampinf[genealogysamples], z,&debug_ti_addinfo[genealogysamples][0]);
      else
		    savegsampinf (gsampinf[genealogysamples], z);
			return;

		}
	} else if (currentid != 0) {
		if ( z >= 0) {
      if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)
        savegsampinf_debug_ti_addinfo (gsampinflocal, z,debug_ti_addinfolocal);
      else
			  savegsampinf (gsampinflocal, z);
		}
	}
#ifdef MPI_ENABLED
	///if i am not on head node, send this info to head node
	if (currentid != 0 && numprocesses > 1 && z >= 0) {
		//send my id first so head node knows where to receive from
			int tempcurrid = currentid;
			rc = MPI_Isend(&tempcurrid, 1, MPI_INT, 0, 1212, MPI_COMM_WORLD, &request[0]);
			if (rc != MPI_SUCCESS)
				MPI_Abort(MPI_COMM_WORLD, rc);

			rc = MPI_Waitall(1, &request[0], &status);
			if (rc != MPI_SUCCESS)
				MPI_Abort(MPI_COMM_WORLD, rc);
		for (int v = 0; v < gsampinflength; v++) {
				rc = MPI_Isend(&gsampinflocal[v], 1, MPI_FLOAT, 0, v*2, MPI_COMM_WORLD, &request[0]);
				if (rc != MPI_SUCCESS)
					MPI_Abort(MPI_COMM_WORLD, rc);
				rc = MPI_Waitall(1, &request[0], &status);
				if (rc != MPI_SUCCESS)
					MPI_Abort(MPI_COMM_WORLD, rc);
		}
    if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)
    {
   		rc = MPI_Isend(debug_ti_addinfolocal , 1000, MPI_CHAR, 0, 3971, MPI_COMM_WORLD, &request[0]);
	  	if (rc != MPI_SUCCESS)
		  	MPI_Abort(MPI_COMM_WORLD, rc);
    }
	} else if (currentid == HEADNODE && numprocesses > 1 && z < 0) {
		//receive info from the other proces first about where to receive from
		int tempcurrid = 0;
			rc = MPI_Irecv(&tempcurrid, 1, MPI_INT, MPI_ANY_SOURCE, 1212, MPI_COMM_WORLD, &request[0]);
			if (rc != MPI_SUCCESS)
				MPI_Abort(MPI_COMM_WORLD, rc);
			rc = MPI_Waitall(1, &request[0], &status);
			if (rc != MPI_SUCCESS)
				MPI_Abort(MPI_COMM_WORLD, rc);
		for (int v = 0; v < gsampinflength; v++) {
			rc = MPI_Irecv(&gsampinf[genealogysamples][v], 1, MPI_FLOAT, tempcurrid, v*2, MPI_COMM_WORLD, &request[0]);
			if (rc != MPI_SUCCESS)
				MPI_Abort(MPI_COMM_WORLD, rc);
			rc = MPI_Waitall(1, &request[0], &status);
			if (rc != MPI_SUCCESS)
				MPI_Abort(MPI_COMM_WORLD, rc);

		}
      if (hiddenoptions[HIDDENGENEALOGY] && hiddenoptions[GSAMPINFOEXTRA] == 1)
      {
   			rc = MPI_Irecv(debug_ti_addinfo[genealogysamples], 1000, MPI_CHAR, tempcurrid, 3971, MPI_COMM_WORLD, &request[0]);
	  		if (rc != MPI_SUCCESS)
		  		MPI_Abort(MPI_COMM_WORLD, rc);
      }
	}
  if (currentid != 0)
    XFREE(gsampinflocal);
#endif  //MPI_ENABLED

	return;
  }
}                               /* savegenealogyinfo */


void loadgenealogyvalues (void)  // only called when current node is HEADNODE
{
  char filenamewildcard[FNSIZE];
  char *ctp, *textline, *dataline, *c, tempc;
             /* CR 111019.1 increase from 12 to 20 handles larger data values */
  int charspervalue = 20;
  FILE *sfile;
  int i, j, numgenealogies, totalnumgenealogies;
  /* int filefound, nofile; */
  char *defaultdir;
  struct dirent *dir_entry;
  DIR *dp;
  int numfiles = 0;
  int numtoload, loadall, loaded, notloaded;
  size_t len_base;
  size_t l2;
  size_t len_defaultdir;
  char *basename;
  char *genealogysavefilename;

  float load_fraction;
  if (strlen (loadfilebase))
  {
    strcpy (filenamewildcard, loadfilebase);
  }
  else
  {
    strcpy (filenamewildcard, outfilename);
    strtrunc (filenamewildcard, (char) '.');
    strtrunc (filenamewildcard, (char) '-');
  }
  strcat (filenamewildcard, "*.ti");
  SP "\nLOAD TREES (L) MODE INFORMATION\n");
  SP "============================================================================\n");
  SP "  Base filename for loading files with sampled genealogies: %s\n", filenamewildcard);
  SP "  Files loaded with sampled genealogies:\n");
  numgenealogyfiles = 0;
  textline = static_cast<char *> (malloc (300 * sizeof (char)));
  dataline = static_cast<char *> (malloc (gsampinflength * charspervalue * sizeof (char)));
  ctp = &textline[0];
  numgenealogies = totalnumgenealogies = 0;

  /* Find the default directory based on loadfilebase.
   * /this/directory/a.out -> defaultdir is /this/directory/
   * /a.out                -> /
   * a.out                 -> ./
   */
  imaDirBase (loadfilebase, &defaultdir);
  len_defaultdir = strlen (defaultdir);

  if ((dp = opendir (defaultdir)) == NULL)
  {
    IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
  }
  basename = strrchr (loadfilebase, '/');
  if (basename == NULL)
    {
      basename = strrchr (loadfilebase, '\\');
      if (basename == NULL)
        basename = loadfilebase;
      else
        basename++;
    }
  else
    {
      basename++;
    }
  len_base = strlen (basename);
  while ((dir_entry = readdir(dp)) != NULL)
  {
    if (!strncmp(dir_entry->d_name, basename, len_base))
    {
      l2 = strlen (dir_entry->d_name);
      if (!strcmp (&dir_entry->d_name[l2 - 3], ".ti"))
      {
        numgenealogyfiles++;
        /* We found one. */
        genealogysavefilename = static_cast<char *> (malloc ((len_defaultdir + l2 + 1) * sizeof (char)));
        sprintf (genealogysavefilename, "%s%s", defaultdir, dir_entry->d_name);

        /* Count the number of gene genealogies of the found file. */
        if ((sfile = fopen (genealogysavefilename, "r")) == NULL)
        {
          IM_err (IMERR_TIFILE, " cannot open .ti file");
        }
        while (fgets (textline, 300, sfile)
               && strstr (textline, "VALUESSTART") == NULL && !feof (sfile));
        numgenealogies = 0;
        while ((tempc = fgetc (sfile)) != EOF)    // count lines
        {
          numgenealogies += (tempc == '\n');
        }
        if (numgenealogies < 1)
        {
          printf ("  *no genealogies loaded from file %s\n", dir_entry->d_name);
          SP "  *no genealogies loaded from file %s\n", dir_entry->d_name);
        }
        else
        {
          printf ("  loaded %d genealogies from genealogy file  %s\n", numgenealogies,
                  dir_entry->d_name);
          SP "  loaded %d genealogies from genealogy file  %s\n", numgenealogies,
            dir_entry->d_name);
        }
        fclose(sfile);
        totalnumgenealogies += numgenealogies;

        XFREE (genealogysavefilename);
      }
    }
  }
  closedir(dp);
  /* We have counted gene genealogies. */

  if (genealogiestosave > 0)
    numtoload = IMIN (totalnumgenealogies, genealogiestosave);
  else
    numtoload = totalnumgenealogies;
  numtoload = IMIN (numtoload, MAXGENEALOGIESTOSAVE);
  memforgenealogiessaved = numtoload;
  gsampinf = static_cast<float **> (malloc (numtoload * sizeof (float *)));
  loadall = numtoload >= totalnumgenealogies;
  load_fraction = (float) numtoload/ (float) totalnumgenealogies;
  loaded = 0;
  notloaded = 0;
  SP "\n  HEADER INFORMATION FROM FIRST GENEALOGY FILE\n");
  SP "  ---------------------------------------\n");
// now go through again and save the genealogies
  /* closedir (dp); */
  if ((dp = opendir (defaultdir)) == NULL)
  {
    IM_err (IMERR_TIFILE, "cannot open directory: %s", defaultdir);
  }
  numfiles = 0;

  while ((dir_entry = readdir(dp)) != NULL)
  {
    if (!strncmp(dir_entry->d_name, basename, len_base))
    {
      l2 = strlen (dir_entry->d_name);
      if (!strcmp (&dir_entry->d_name[l2 - 3], ".ti"))
      {
        numfiles++;
        /* We found one. */
        genealogysavefilename = static_cast<char *>
                (malloc ((len_defaultdir + l2 + 1) * sizeof (char)));
        sprintf (genealogysavefilename, "%s%s", defaultdir, dir_entry->d_name);

        /* Count the number of gene genealogies of the found file. */
        if ((sfile = fopen (genealogysavefilename, "r")) == NULL)
        {
          IM_err (IMERR_TIFILE, " cannot open .ti file");
        }

        while (fgets (textline, 300, sfile) && strstr (textline, "VALUESSTART") == NULL && !feof (sfile))
          if (numfiles == 1)
          {
            SP "  ||%s", textline);
          }
        while (fgets (dataline, gsampinflength * charspervalue, sfile)
               && !feof (sfile))
        {
          if (loadall || (loaded == 0)
              ||  ((float) loaded/ (float) (loaded + notloaded)) <= load_fraction)
              //JH 7/24/09  change this to using a proportion((float) loaded / (float) notloaded) <= load_notload_ratio)
          {
            gsampinf[loaded] = static_cast<float *>
                    (malloc (gsampinflength * sizeof (float)));
            c = dataline;
            for (i = 0; i < gsampinflength; i++)
            {
              scanfval = sscanf (c, "%f", &gsampinf[loaded][i]);
              j = allwhitespace (c);
              if (j ==1 || j== -1)
                IM_err (IMERR_TIFILE, "Problem in .ti file %s,  too few values per genealogy, .ti file may have been generated with a different program",genealogysavefilename);
              c = nextwhite (c);
            }
            j = allwhitespace (c);
            if (j==0)
               IM_err (IMERR_TIFILE, "Problem in .ti file %s,  too many values per genealogy, .ti file may have been generated with a different program",genealogysavefilename);
            loaded++;
          }
          else
            notloaded++;
          /* gcounter++; */
        }
        fclose (sfile);

        XFREE (genealogysavefilename);
      }
    }
  }
  closedir(dp);

  SP "  END OF HEADER INFORMATION FROM FIRST GENEALOGY FILE\n");
  SP "  ----------------------------------------------\n\n");

  SP "  Number of files loaded : %d  total number of genealogies used: %d out of a total of: %d \n", numgenealogyfiles, loaded, totalnumgenealogies);
  printf("  Number of files loaded : %d  total number of genealogies used: %d out of a total of: %d \n", numgenealogyfiles, loaded, totalnumgenealogies);
  fflush(stdout);
  if (numgenealogies < 1)
    IM_err (IMERR_TIFILE, "  no genealogies loaded from .ti file(s)");
  /* closedir (dp); */
  genealogysamples = loaded;

  for (j = 0; j < genealogysamples; j++)
  {
    // use full range of t, ignore t.pr.min > 0
    if (npops > 1)
    {
      for (i = 0; i < lastperiodnumber; i++)
      {
        recordval (T[i].v, gsampinf[j][gsamp_tp + i]);
      }
    }
    if (!outputoptions[DONTPRINTASCIITREND])
      trendrecord (j, 0);
  }
// }
  XFREE (textline);
  XFREE (dataline);
  XFREE (defaultdir);
}                               /* sang chul's loadgenealogyvalues */

/* reorganized printoutput()  6/13/2017 */
void
printoutput (int currentid, int finaloutput)         // mostly calls functions in output.c   - called only for printing a results file
{
  int i;
  long seconds;
  int p;
  float *holdpeakloc;
  double multitpeak[MAXPOPS - 1];
  struct tm *endtimeinfo;
  static time_t startoutputtime,endoutputtime,totaloutputseconds = 0,previoustime = 0;
  int rc = 0;
  time(&startoutputtime);

  #ifdef XMLOUTPUT
  char s[200];
  int n = sprintf(s,"%s.xml",outfilename);
  TiXmlDocument doc;
  TiXmlDeclaration *decl = new TiXmlDeclaration("1.0","","");
  doc.LinkEndChild(decl);
  TiXmlElement *output = new TiXmlElement("Output");
  doc.LinkEndChild(output);
  xstack.push(output);
  #endif

/*   printoutput()  open outout file,  deal with old outputfile */
  if (runoptions[LOADRUN] == 0)
  {
    if ((cdurationmode == TIMEINF || runoptions[SAVELOADSAMEMCFFILE]) && currentid == HEADNODE)
    {
      if (file_exists(oldoutfilename))
        remove (oldoutfilename);
      rename (outfilename, oldoutfilename);
    }
  }
  if (currentid == HEADNODE)
  {
    if ((outfile = fopen (outfilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL, "Error opening text file for writing");
    }
  }
  else
    outfile = NULL;

/*   printoutput()  call printrunbasics()*/
  runsteps = step - burninsteps;
  printrunbasics (outfile, runoptions[LOADRUN], fpstr, burninsteps, burninsteps_old,runsteps_old,mcmcrecords_old,genealogysamples_old,recordint,mcmcrecords, savegenealogyint, hilike,hiprob);

#ifdef STDTEST
printf("printed run basics\n");
#endif

/*   printoutput()  print acceptance rates and autoc table*/
  if (runmode != LOADGmode4) // any mode running mcmc
  {
    if (currentid == HEADNODE)
    {
      FP "%s",outputbanner("MCMC Information"));
    }

    callprintacceptancerates (outfile, currentid);
#ifdef STDTEST
printf("printed acceptance rates\n");
#endif
    if (numchainstotal > 1)
    {
	     printchaininfo (outfile, heatmode, hval1, hval2, currentid);
    }
    if (currentid == HEADNODE)
    {
      callprintautoctable (outfile/*, step*/);
      /* get joint splittime peak */

      if ((runmode == POPTREEHYPERPRIORmode0 || runmode == POPTREEmode1) && burndone && hiddenoptions[CALCGEWEKEZ])
        printgewekez(outfile);
    }
#ifdef STDTEST
printf("printed autoc table\n");
#endif

/*   printoutput()  share topology counts across chains, print topology posterior*/
    if (runmode == POPTREEHYPERPRIORmode0 || runmode == POPTREEmode1)
    {
/*#ifdef INDEVELOPMENT
    //output_update_scalars(whichiscoldchain(),currentid, updatestr); // gets update scalar info from chain 0,  puts it in updatestr and sends updatestr to cpu 0
    //if (currentid == HEADNODE)
      //FP "%s",updatestr);
    //printupdatescalarinfo (outfile,currentid);  stopped using 6/7/2017
#endif */

      int *totalpoptopologycounts = static_cast<int *> (malloc ((size_t) numpoptopologies * sizeof (int))); // had calloc here and it would not work with mpi_reduce ??
      int *poptopologyproposedlist_rec = static_cast<int *> (malloc ( (size_t) numpoptopologies * sizeof (int)));
      int rc;
      if (numprocesses > 1)
      {
#ifdef MPI_ENABLED
        rc = MPI_Reduce(poptopologycounts, totalpoptopologycounts, numpoptopologies,MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
			     if (rc != MPI_SUCCESS)  MPI_Abort(MPI_COMM_WORLD, rc);
        rc = MPI_Reduce(poptopologyproposedlist, poptopologyproposedlist_rec, numpoptopologies , MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
#endif
      }
      else
      {

        for (i = 0;i<numpoptopologies;i++)
        {
          totalpoptopologycounts[i] = poptopologycounts[i];
          poptopologyproposedlist_rec[i] = poptopologyproposedlist[i];
        }
      }
      if (currentid == HEADNODE)
      {
        poptopologiessampled = 0;
        for (i = 0;i<numpoptopologies;i++)
          poptopologiessampled += totalpoptopologycounts[i];


        sort_and_print_alltreestrings(outfile, totalpoptopologycounts,poptopologyproposedlist_rec,&topologypriorinfostring[0]);
      }
      XFREE(poptopologyproposedlist_rec);
      XFREE(totalpoptopologycounts);  // uncommented this,  why was it commented out ?? 4/20/2017

#ifdef  TURNONCHECKS
      /* check if topology count matches mcmcrecords  // does not apply if loading mcf files
      int tsum = 0;
      for (i=0;i<numpoptopologies;i++)
      {
        tsum +=  totalpoptopologycounts[i];
        assert (totalpoptopologycounts[i] >= 0);
      }
      assert (tsum == mcmcrecords); */
      /* check that no two beta values are identical
      int j;
      for (i=0;i<numchainspp-1;i++)
        for (j=i+1;j<numchainspp;j++)
          assert (beta[i] != beta[j]); */
#endif  //TURNONCHECKS
    }
  }
  if (currentid == HEADNODE &&(runmode == Gmode3 || runmode == LOADGmode4 || runmode == HGmode6))
  {
    if (npops >=3  && npops <= 5 && !runoptions[LOADRUN] && outputoptions[PRINTJOINTTEST])
    {
/*   printoutput()  joint t distribution*/
      return_joint_t (multitpeak);
      FP "\nEstimated joint splitting time from multi-dimensional histogram\n");
      FP "  number of bins per dimension %d\n", NUMTARRAYBINS);
      FP "  Posterior probability of estimated joint value of splitting time: %7.4lf\n", joint_t_prob (&multitpeak[0]));
      FP "---------------------------------------------------------------\n");
      for (i = 0; i < numsplittimes; i++)
        FP "   %s\t%.3lf\n", T[i].str, multitpeak[i]);
      FP "\n\n");
    }
/*   printoutput()  greater thans,  means and variances */
    if (outputoptions[PARAMGREATERTHAN])
    {
      print_greater_than_tests (outfile);
#ifdef STDTEST
printf("printed greater than tests\n");
#endif
    }
    if (!modeloptions[EXPOMIGRATIONPRIOR] && calcoptions[DONTCALCLIKELIHOODMUTATION]==0)        //as of 11/19/09 have not yet done the math for case of migration with exponential prior
      print_means_variances_correlations (outfile);
#ifdef STDTEST
printf("printed means variances correlations\n");
#endif
/*   printoutput()  get marginal peaks*/
    if (!calcoptions[DONTCALCLIKELIHOODMUTATION])
    {
      p = numpopsizeparams + nummigrateparams ;
      holdpeakloc = static_cast<float *> (malloc (p * sizeof (float)));
      printf ("surface calculations  . . .\n");
      fflush(stdout);
      if (modeloptions[EXPOMIGRATIONPRIOR] || (runoptions[LOADRUN] && calcoptions[FINDJOINTPOSTERIOR]))
        eexpsum = (struct extendnum *) malloc ((size_t) ((genealogysamples + 1) * sizeof (struct extendnum)));

      findmarginpeaks (outfile, holdpeakloc);
/*   printoutput()  get joint peaks*/
      if (runmode == LOADGmode4 && calcoptions[FINDJOINTPOSTERIOR])
      {
        closeopenout (&outfile, outfilename);
        /* CR 110921.1  Change outfile parameter type to (FILE **) */
        findjointpeaks(&outfile,outfilename,nestedmodelfilename,p);
#ifdef STDTEST
printf("printed joint peaks\n");
#endif
      }
      XFREE (holdpeakloc);
#ifdef STDTEST
printf("found peaks\n");
#endif
    }
  }
/*   printoutput()  print histograms*/
  if (currentid == HEADNODE)
  {

    printhistograms (outfile, mcmcrecords, generationtime,usegenerationtimedefault, scaleumeaninput,priorfilename);
#ifdef STDTEST
printf("printed histograms\n");
#endif
/*   printoutput()  print ascii trends*/
    FP "\n\n===================================================\n");
    if (!outputoptions[DONTPRINTASCIITREND] && runmode != LOADGmode4)
    {
      if (trendspot <= 1)
      {
        FP "run too short to plot trends\n");
      }
      else
      {
        FP "%s",outputbanner("ASCII Plots of Parameter Trends"));
        FP " - note points to the left of '!' on X axis have twice the density in time relative to points to the right\n\n");
        callasciitrend (outfile,trenddoublepoint,trendspot);
#ifdef STDTEST
printf("printed trends\n");
#endif
      }
    }
    callasciitrend (outfile,trenddoublepoint,trendspot);
#ifdef STDTEST
printf("printed trends\n");
#endif
/*   printoutput()  print ascii curves*/
    FP "%s",outputbanner("ASCII Curves - Approximate Posterior Densities"));

    callasciicurves (outfile,mcmcrecords);
#ifdef STDTEST
printf("printed ascii curves\n");
#endif
  }
/*   printoutput()  print migration histogram file  - almost certainly broken as of 1/17/2018*/
  if (outputoptions[MIGRATEHIST] && nummigrateparams > 0  && currentid == HEADNODE)
  {
    if ((migplotfile = fopen (migplotfilename, "w")) == NULL)
    {
      IM_err (IMERR_CREATEFILEFAIL,
              "Error opening file for plotting migration amounts and times");
    }
    printmigrationhistograms (migplotfile, mcmcrecords);
  #ifdef STDTEST
  printf("printed migration histogram\n");
  #endif
    FCLOSE (migplotfile);
  }
  if (hiddenoptions[WRITEMIGRATIONNAME])
  {
    FCLOSE(migrationnamefile);
    migrationnamefile = fopen (migrationnamefilename, "a"); // file is kept open
  }
/*   printoutput()  print timing information*/
  time (&endtime);
  endoutputtime = endtime;
  totaloutputseconds += difftime(endoutputtime,startoutputtime);
  endtimeinfo = localtime(&endtime);
  seconds = difftime (endtime, starttime);
  if (numpriormcfruns > 0  && finaloutput)
  {
    previoustime = totaltime;
    totaltime += (time_t) seconds;
  }
  else
  {
    totaltime = (time_t) seconds;
  }
  if (currentid == HEADNODE)
  {
    FP "\nTime Elapsed : %s \n",timestring(seconds));
    FP "Time spent running analyses after mcmc finished : %s\n",timestring(totaloutputseconds));
/*    FP "\nTime Elapsed : %d hours, %d minutes, %d seconds \n",
      (int) seconds / (int) 3600,
        ((int) seconds / (int) 60) - ((int) 60 * ((int) seconds / (int) 3600)),
            (int) seconds - (int) 60 *((int) seconds / (int) 60));
    FP "Time spent running analyses after mcmc finished : %d hours, %d minutes, %d seconds \n",
      (int) totaloutputseconds / (int) 3600,
        ((int) totaloutputseconds / (int) 60) - ((int) 60 * ((int) totaloutputseconds / (int) 3600)),
            (int) totaloutputseconds - (int) 60 *((int) totaloutputseconds / (int) 60)); */
    if (numpriormcfruns > 0  && finaloutput)
    {
      FP "Duration of previous runs : %s \n",timestring(previoustime));
      FP "Total Time of all Runs : %s\n",timestring(totaltime));
  /*    FP "Duration of previous runs : %d hours, %d minutes, %d seconds \n",
        (int) previoustime / (int) 3600,
          ((int) previoustime / (int) 60) - ((int) 60 * ((int) previoustime / (int) 3600)),
              (int) previoustime - (int) 60 *((int) previoustime / (int) 60));
      FP "Total Time of all Runs : %d hours, %d minutes, %d seconds \n",
        (int) totaltime / (int) 3600,
          ((int) totaltime / (int) 60) - ((int) 60 * ((int) totaltime / (int) 3600)),
              (int) totaltime - (int) 60 *((int) totaltime / (int) 60)); */
    }
    strftime(timeinfostring,80,"%x - %X", endtimeinfo);
    FP "\nJob Finished: %s\n",timeinfostring);
    FP "\nEND OF OUTPUT\n");
    FCLOSE (outfile);
    fflush(stdout);
  }
#ifdef MPI_ENABLED
		if (numprocesses > 1)
  {
    rc = MPI_Bcast(&totaltime, time_t_size, MPI_BYTE, 0, MPI_COMM_WORLD); // broadcast from 0 to others // all processes must reach this line
	 	 if (rc !=MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,-1);
	 }
#endif
  if (runoptions[SAVEMCSTATEFILE])
  {
    //jh 1_17_2018 netsteps += step;
    //writemcf (mcfwritefilename,command_line,mcmcrecords,hilike,hiprob,currentid);
    writemcf (mcfwritefilename,command_line,mcmcrecords,mcmcrecords_old,genealogysamples_old,burninsteps_old,runsteps_old,hilike,hiprob,currentid);
  }
#ifdef STDTEST
printf("done printing output\n");
#endif
#ifdef XMLOUTPUT
if (currentid == HEADNODE) {
//doc.LinkEndChild(runbasics);
doc.SaveFile(s);
xstack.pop();
}
#endif
  return;
}                               /* printoutput */

/* this is generally just turned on when in development */
void output_update_scalars(int z,int currentid, char updatestr[])
{
  static int updatestri = 0;
  int ti;
  int zid  = -1;
#ifdef MPI_ENABLED
  MPI_Status status;
#endif
  int rc = 0;
#define UP updatestri += sprintf(&updatestr[updatestri],

  if (z >= 0)
  {
    zid = currentid;
    updatestri = 0;
    UP "\nPopulation Tree Update Scalars, chain 0\n---------------------------------------\n");
    if (doRYupdate) for (ti=0;ti<numsplittimes;ti++)
    {
      UP "RY update, time period %d, update scalar %.4f  overall acceptance rate %.4f \n",ti,C[z]->RYwidthinfo[ti].updatescalarval,C[z]->RYwidthinfo[ti].allaccp/(float) C[z]->RYwidthinfo[ti].alltries);
    }
    if (doNWupdate) for (ti=0;ti<numsplittimes;ti++)
    {
      UP "NW update, time period %d, update scalar %.4f  overall acceptance rate %.4f \n",ti,C[z]->NWwidthinfo[ti].updatescalarval,C[z]->NWwidthinfo[ti].allaccp/(float) C[z]->NWwidthinfo[ti].alltries);
    }
    /* 5/3/2017  not using popslideinfo.updatescalarval
    if (hiddenoptions[HIDDENGENEALOGY] == 1)
    {
      UP "Branch slide update, update scalar %.4f  overall acceptance rate %.4f \n",C[z]->branchslideinfo.updatescalarval,C[z]->branchslideinfo.allaccp/(float) C[z]->branchslideinfo.alltries);
    }*/
#ifdef MPI_ENABLED
    if (currentid != 0 ) //current process is not zero but has updatestr, must send to 0
    {
		   //rc = MPI_Send(&updatestr,1000, MPI_CHAR, 0, 12377, MPI_COMM_WORLD); // fixed bug  tag 123 was used again for an unrelated recv so changed to 12377
     rc = MPI_Send(updatestr,1000, MPI_CHAR, 0, 12377, MPI_COMM_WORLD); // fixed bug  tag 123 was used again for an unrelated recv so changed to 12377
		    if (rc !=MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,-1);
    }
#endif
  }
#undef UP
  if (currentid ==0 && z < 0)// current process is 0, must receive updatestr
  {
#ifdef MPI_ENABLED
		  //rc = MPI_Recv(&updatestr,1000, MPI_CHAR, MPI_ANY_SOURCE, 12377, MPI_COMM_WORLD,&status); // fixed bug  tag 123 was used again for an unrelated recv so changed to 12377
    rc = MPI_Recv(updatestr,1000, MPI_CHAR, MPI_ANY_SOURCE, 12377, MPI_COMM_WORLD,&status); // fixed bug  tag 123 was used again for an unrelated recv so changed to 12377
		  if (rc !=MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD,-1);
#endif  //MPI_ENABLED
  }
}  // output_update_scalars


#ifdef MPI_ENABLED
/* reduction operation for updatescalarinfo -  ignore type, just trust that it's our struct type */
void jh_mpi_sum_scalerstruct(struct updatescalarinfo *in,struct updatescalarinfo *inout, int *len, MPI_Datatype *type){
  for (int i=0; i<*len; i++)
  {
    inout[i].updatescalarval += in[i].updatescalarval;
    inout[i].updatescalaradjustval += in[i].updatescalaradjustval;
    inout[i].updatescalarvalmin += in[i].updatescalarvalmin;
    inout[i].updatescalarvalmax += in[i].updatescalarvalmax;
    inout[i].targetupdaterate += in[i].targetupdaterate;
    inout[i].recentaccp += in[i].recentaccp;
    inout[i].recenttries += in[i].recenttries;
    inout[i].allaccp += in[i].allaccp;
    inout[i].alltries += in[i].alltries;
    inout[i].numattemptscheck += in[i].numattemptscheck;
  }
  return;
}
#endif

/*void printupdatescalarinfo (FILE * outto, int currentid) stopped using 6/7/2017
{
  int z = whichiscoldchain();
  int rc = 0;
  struct updatescalarinfo suminfo,tempinfo, *siarray[2*MAXPOPS];
  //static int dobs=0;
  int numtodo=0;
  int i,j;
  char sistrings[2*MAXPOPS][30]; //hold strings
  static struct updatescalarinfo holdinfo0 = {0.0,0.0,0.0,0.0,0.0,0,0,0,0,0};

  i = 0;
  if (doRYupdate)
  {
    numtodo  += numsplittimes;
    for (j=0;j<numsplittimes;j++)
    {
      if (z >= 0)
        siarray[i] = &(C[z]->RYwidthinfo[j]);
      else
        siarray[i] = &holdinfo0;
      sprintf(sistrings[i],"RY period %d :",j);
      i+=1;
    }
  }
  if (doNWupdate)
  {
    numtodo  += numsplittimes;
    for (j=0;j<numsplittimes;j++)
    {
      if (z >= 0)
        siarray[i] = &(C[z]->NWwidthinfo[j]);
      else
        siarray[i] = &holdinfo0;
      sprintf(sistrings[i],"NW period %d :",j);
      i+=1;
    }
  }
  if (currentid == HEADNODE)
  {
    fprintf(outto,"\n=========================\nUpdate Scalar Information\n=========================\n");
    fprintf(outto,"update      \tscalar\tadjust\tmin\tmax\ttarget\t#recAcp\t#recTry\t#Acp\t#Try\t#check\n");

  }
  for(i=0;i<numtodo;i++)
  {
    if (numprocesses == 1)
    {
        suminfo = *siarray[i];
    }
    else
    {
      suminfo = holdinfo0;
      tempinfo = *siarray[i];
#ifdef MPI_ENABLED
      rc = MPI_Reduce(&tempinfo,&suminfo,1, MPI_updatescalar,myOp_updatescalarsum, 0, MPI_COMM_WORLD);
		    if (rc != MPI_SUCCESS)
        MPI_Abort(MPI_COMM_WORLD, rc);
#endif
    }
    if (currentid == HEADNODE)
    {
      fprintf(outto,"%s\t%0.4lf\t%0.4lf\t%0.4lf\t%0.4lf\t%0.4lf\t%d\t%d\t%d\t%d\t%d\n",sistrings[i],suminfo.updatescalarval,suminfo.updatescalaradjustval,suminfo.updatescalarvalmin,suminfo.updatescalarvalmax,suminfo.targetupdaterate,suminfo.recentaccp,suminfo.recenttries,suminfo.allaccp,suminfo.alltries,suminfo.numattemptscheck);
    }
  }
} //printupdatescalarinfo  stopped using 6/7/2017 */

/* this is called from within qupdate with results to be sent to stdout
  and it is called from writeburntrendfile()  with results to be sent to the burntrendfile
  also gathers topologycounts to cpu 0 from other cpus */

/*JH 4/27/2016  moved the main step condition outside of this so it would not get called so often, to save some waiting in mpi */
void intervaloutput (FILE * outto, int currentid)
{
  int rc = 0; //AS: Return code for MPI C bindings
  int poptreenum_rec;
#ifdef XMLOUTPUT
  char s[200];
  int n = sprintf(s,"%s.intervals.xml",outfilename);
  TiXmlDocument doc;
  TiXmlDeclaration *decl = new TiXmlDeclaration("1.0","","");
  doc.LinkEndChild(decl);
  TiXmlElement *intervals = new TiXmlElement("Interval");
  doc.LinkEndChild(intervals);
  xstack.push(intervals);
#endif
#ifdef MPI_ENABLED
	MPI_Status status;
  int totaltopol_rec, chain0topol_rec, chain0topolswaps_rec;
#endif
  int z = whichiscoldchain();
  checkhighs(z,currentid,0);
  runsteps = step - burninsteps;
  if (currentid == HEADNODE)
    //printsteps (outto, like_rec, probg_rec, burndone,burninsteps);
    printsteps (outto, currlike, currprob, burndone,burninsteps);
  if (z >=0 && currentid == HEADNODE)
    poptreenum_rec = C[z]->poptreenum;
#ifdef MPI_ENABLED
  else
  {
    if (z >=0 && currentid != 0)
    {
      rc = MPI_Send(&C[z]->poptreenum, 1, MPI_INT, 0, 13137, MPI_COMM_WORLD);
	     if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
    }
    if (z < 0 && currentid == HEADNODE)
    {
      rc = MPI_Recv(&poptreenum_rec, 1,MPI_INT, MPI_ANY_SOURCE, 13137, MPI_COMM_WORLD, &status);
	     if (rc != MPI_SUCCESS) MPI_Abort(MPI_COMM_WORLD, rc);
    }
  }
  #endif
   if (modeloptions[POPTREETOPOLOGYUPDATE]==1) // do the reduction for all processors
   {
#ifdef MPI_ENABLED
  if (numprocesses > 1)
  {
	  rc = MPI_Reduce(&totaltopolupdates, &totaltopol_rec, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rc != MPI_SUCCESS)
			MPI_Abort(MPI_COMM_WORLD, rc);
  	rc = MPI_Reduce(&chain0topolupdates, &chain0topol_rec, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	  if (rc != MPI_SUCCESS)
				MPI_Abort(MPI_COMM_WORLD, rc);
		rc = MPI_Reduce(&chain0topolswaps, &chain0topolswaps_rec, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		if (rc != MPI_SUCCESS)
				MPI_Abort(MPI_COMM_WORLD, rc);
   }
#endif
    if (numchainstotal > 1 && currentid == HEADNODE)
    {
      if (outto != NULL) fprintf(outto,"\nCurrent Population Topology #: %d   Current Topology String: %s\n",poptreenum_rec,alltreestrings[poptreenum_rec]);
      if (outto != NULL) fprintf(outto,"\nPopulation Topology Updates Across Chains\n-----------------------------------------\n");
      if (numprocesses > 1)
      {
#ifdef MPI_ENABLED
		    if (outto != NULL) fprintf(outto," Total Number of Accepted Topology Updates Across all Chains: %d\n",totaltopol_rec);
	      if (outto != NULL) fprintf(outto," Number of Accepted Topology Updates to the Cold Chain: %d\n",chain0topol_rec);
        if (outto != NULL) fprintf(outto," Number of Accepted Swaps that Changed Topology for the Cold Chain: %d\n",chain0topolswaps_rec);
#endif
      }
      else
      {
        if (outto != NULL) fprintf(outto," Total Number of Topology Updates Across all Chains: %d\n",totaltopolupdates);
        if (outto != NULL) fprintf(outto," Number of Accepted Topology Updates to the Cold Chain: %d\n",chain0topolupdates);
	      if (outto != NULL) fprintf(outto," Number of Accepted Swaps that Changed Topology for the Cold Chain: %d\n",chain0topolswaps);
      }
      if (outto != NULL) fprintf(outto,"\n");
    }
    if (numchainstotal == 1) // 1 chain, currentid must be 0
    {
      if (outto != NULL) fprintf(outto,"\nTotal Number of Topology Updates: %d\n",totaltopolupdates);
    }
/*#ifdef INDEVELOPMENT
      /// causing some MPI crash on MUST ??
   //output_update_scalars(whichiscoldchain(),currentid, updatestr); // gets update scalar info from chain 0,  puts it in updatestr and sends updatestr to cpu 0
   // if (currentid == HEADNODE)
    //  if (outto != NULL) fprintf(outto,"%s",updatestr);
    //printupdatescalarinfo (outto,currentid);   stopped using 6/7/2017
#endif */
   }
   if (hiddenoptions[STOPMOSTINTERVALOUTPUT]==0) // turning off chunk of interval output to reduce screen text when debugging
   {

     if (numprocesses > 1)
        printcurrent_tvals (stdout,currentid);
     else
        printcurrentvals (outto); // only works with 1 processor   4/19/2017
     callprintacceptancerates (outto, currentid);
     if (currentid == HEADNODE)
     {
        callprintautoctable (outto );
        if (burndone && modeloptions[POPTREETOPOLOGYUPDATE] && hiddenoptions[CALCGEWEKEZ])
          printgewekez(outto);
     }
     if (numchainstotal > 1)
     {
        printchaininfo (outto, heatmode, hval1, hval2, currentid);

     }
   } // hiddenoptions[STOPMOSTINTERVALOUTPUT]==0
  if (outto==stdout)
    fflush(stdout);
#ifdef XMLOUTPUT
  else if (currentid == HEADNODE) {
      doc.SaveFile(s);
      xstack.pop();
  }
 #endif
  return;
}                               /* intervaloutput */

// check if it is time to call record() and savegenealogyinfo(), and call if it is
void check_to_record (int currentid)
{
  static int i;
  static int j;
  static int init = 0;

  if (init == 0)
  {
    i = recordint;
    j = savegenealogyint; // start out by saving a genealogy
    init = 1;
  }
  if (i == recordint)
  {
    /*if (calcoptions[CALCMARGINALLIKELIHOOD])
        stepstone_get_Lmax(); */
    record (currentid);                  // record some values that are in mcmc

    mcmcrecords++;
    i = 1;
  }
  else
  {
    i++;
  }

  /* only save genealogies in mode 3 or mode 1 if GSAMPINFOEXTRA */
  if (j == savegenealogyint && (runmode == Gmode3 || runmode == HGmode6 ||(runmode == POPTREEmode1 && hiddenoptions[GSAMPINFOEXTRA]==1)))
  {
    savegenealogyinfo (currentid);            // record values associated with genealogies
    if (hiddenoptions[WRITEMIGRATIONNAME])
      record_migration_names();
    genealogysamples++;
    j = 1;
  }
  else
  {
    j++;
  }
}                               // check_to_record

void commit_mpi_updatescalar(void)
{
#ifdef MPI_ENABLED
  /* makes an MPI holder (MPI_updatescalar) for struct updatescalarinfo
  struct updatescalarinfo
  {
    double updatescalarval;  // the current scalar
    double updatescalaradjustval;  //  amount to adjust the scalar
    double updatescalarvalmin;
    double updatescalarvalmax;
    double targetupdaterate;  // target update rate
    int recentaccp;
    int recenttries;
    int allaccp;
    int alltries;
    int numattemptscheck;   // the number of attempts to use for calculating the recent acceptance rate
  }; */
  //MPI_Datatype MPI_updatescalar;  // declared at top of this file, and with extern in swapchains
  int mpi_updatescalar_blocklengths[10] = {1,1,1,1,1,1,1,1,1,1};
  MPI_Datatype mpi_updatescalar_types[10]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT};
  MPI_Aint mpi_updatescalar_displacements[10];
  MPI_Aint intex,doublex;
  MPI_Type_extent(MPI_INT,&intex);
  MPI_Type_extent(MPI_DOUBLE,&doublex);
  mpi_updatescalar_displacements[0] = (MPI_Aint) 0;
  mpi_updatescalar_displacements[1] = doublex;
  mpi_updatescalar_displacements[2] = 2*doublex;
  mpi_updatescalar_displacements[3] = 3*doublex;
  mpi_updatescalar_displacements[4] = 4*doublex;
  mpi_updatescalar_displacements[5] = 5*doublex;
  mpi_updatescalar_displacements[6] = 5*doublex + intex;
  mpi_updatescalar_displacements[7] = 5*doublex + 2*intex;
  mpi_updatescalar_displacements[8] = 5*doublex + 3*intex;
  mpi_updatescalar_displacements[9] = 5*doublex + 4*intex;
  MPI_Type_struct(10,mpi_updatescalar_blocklengths,mpi_updatescalar_displacements,mpi_updatescalar_types,&MPI_updatescalar);
  MPI_Type_commit(&MPI_updatescalar);

#endif
}

/*
  one cpu is the head node and handles the output  this has currentid value of 0 == HEADNODE
  currentid has nothing to do with the where the cold chain is
  this could be on any node

  the head node (i.e. cpu with currentid == HEADNODE) is used for almost all output operations (exception is mcf files - every node has one of those)
  the head node broadcasts, sends and receives stuff from the others as needed
*/
int main (int argc, char *argv[])
{
  int rc;
  int currentid; // this gets used a lot,  so could be global. but keep passing it around as needed to make things more clear as to when it is needed
#ifdef MPI_ENABLED
  rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS) 	MPI_Abort(MPI_COMM_WORLD, rc);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocesses);
	 MPI_Comm_rank(MPI_COMM_WORLD, &currentid);
  commit_mpi_updatescalar();
  MPI_Op_create((MPI_User_function *) jh_mpi_sum_scalerstruct, true, &myOp_updatescalarsum); // makes the jh_mpi_sum_scalerstruct function work wtih MPI_Reduce
#else
  numprocesses = 1;
  currentid = 0;
#endif

#ifdef TURNONCHECKS
#ifdef MPI_ENABLED
    MPI_Comm_rank(MPI_COMM_WORLD, &currentid_debug); // just to have a global access to the current cpu #  when debugging
#else
    currentid_debug = 1;
#endif
#endif

  start (argc, argv, currentid);
  if (runoptions[LOADRUN])
  {
	   if (currentid == HEADNODE)
    {
	    loadgenealogyvalues ();
	    mcmcrecords = genealogysamples; // why is mcmcrecords needed when runoptions[LOADRUN] ?
	  }
   printoutput (currentid,0);
	  if (currentid == HEADNODE)
	    free_ima_main_stuff ();
  }
  else  // !runoptions[LOADRUN]
  {
    if (currentid == HEADNODE)
    {
  	   printf ("Starting Markov chain.\n");
      fflush(stdout);
    }

    step = 0;
    mcmcrecords = 0;

    while (run (currentid))  // main mcmc loop
    {
      qupdate (currentid);
      if ((step / (int) printint) * (int) printint == step && step > 0)
        intervaloutput (stdout, currentid);
      if (step >= CHECKAUTOCWAIT)  // moved from qupdate()
        checkautoc (0, burndone, burninsteps, currentid);
      if (burndone)
      {
        check_to_record (currentid);
      }
      else
      {
        if (runoptions[PRINTBURNTREND] && modeloptions[POPTREETOPOLOGYUPDATE])
          recordburntopology();
      }
#ifdef TURNONCHECKS
//    if (burndone)
  //    chaininfo_print(currentid,recordint);
#endif //TURNONCHECKS
	     step++;
    }  //  run() mcmc loop

#ifdef STDTEST
printf("done mcmc\n");
#endif
#ifdef MPI_ENABLED
	   MPI_Barrier(MPI_COMM_WORLD);  // is this necessary ?
#endif

 /* save genealogy info in *.ti file */
    if (runoptions[DONTSAVEGENEALOGIES]==0 && genealogysamples > 0 && currentid == HEADNODE)
    {
      savegenealogyfile (genealogyinfosavefilename, genealogyinfosavefile, &lastgenealogysaved, gsampinflength);
#ifdef STDTEST
printf("saved genealogies\n");
#endif
    }
    printoutput (currentid,1);
    if (hiddenoptions[WRITEMIGRATIONNAME])
    {
      FCLOSE(migrationnamefile);
    }
    free_ima_main_stuff ();
#ifdef STDTEST
printf("freed stuff\n");
#endif
  }

#ifdef MPI_ENABLED
  MPI_Op_free(&myOp_updatescalarsum);
  MPI_Type_free(&MPI_updatescalar);  //1_26_2018
  MPI_Finalize();
#ifdef STDTEST
printf("finalized mpi\n");
#endif

#endif

  return 0;
}                               /* main */
