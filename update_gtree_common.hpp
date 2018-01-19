/*IMa3 2018 Jody Hey, Rasmus Nielsen, Sang Chul Choi, Vitor Sousa, Janeen Pisciotta, Yujin Chung and Arun Sethuraman */
#ifndef _UPDATE_GTREE_COMMON_H_
#define _UPDATE_GTREE_COMMON_H_
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS          /* empty */
# define __END_DECLS            /* empty */
#endif

__BEGIN_DECLS
/* updating ancestral alleles at nodes,  under stepwise model */
#define updateAfrac 0.05        // inverse of updateAi, just a useful proportion,  don't need to do all of the nodes every step
#define updateAi  20
void init_update_assignment (void);
void free_update_assignment (void);

/* prototypes */

double likelihoodDG (int ci, int li);
int findperiod_Assignment (int ci, double t);
void treeweight_Assignment (int ci, int li);
void integrate_tree_prob_Assignment (int ci,
                                     struct genealogy_weights *gweight,
                                     struct probcalc *pcalc);
void copyfraclike_Assignment (int ci, int li);
void storescalefactors_Assignment (int ci, int li);
void restorescalefactors_Assignment (int ci, int li);
double finishSWupdateA_Assignment (int ci, int li, int ai, int edge,
                                   int downedge, int sisedge, int newsisedge,
                                   double u, double *Aterm);
double updateA_Assignment (int ci, int li, int ai, double u, int *count);
void setcopyedge(void);
void init_gtreecommon (void);
void free_gtreecommon (void);
void  init_gtreecommonhg (void);         // initialize copyedgehg
void free_gtreecommonhg (void);  
double calcmrate (int mc, double mt);
int joinsisdown (int ci, int li, int sis, int *tmrcachange);
void splitsisdown (int ci, int li, int slidingedge, int down, int newsis);
int getm (int ci, struct edgemiginfo *edgem,struct edgemiginfo *sisem, struct edgemiginfo *oldedgem,struct edgemiginfo *oldsisem);
void slider (int ci, int li, int slidingedge, int *sis, double *timepoint,
             double *slidedist);
void slider_nomigration (int ci, int li, int slidingedge, int *sis,
                         double *timepoint, double *slidedist);
void IMA_reset_edgemiginfo (struct edgemiginfo *em);
void storeoldedges (int ci, int li, int edge, int sisedge, int downedge);
void restoreedges (int ci, int li, int edge, int sisedge, int downedge,
                   int newsisedge);
double getmprob(int ci, struct edgemiginfo *edgem,
          struct edgemiginfo *sisem,struct edgemiginfo *oldedgem,
          struct edgemiginfo *oldsisem);
void storeAinfo (int li, struct edge *gtree, int edge, int sisedge,
                 int downedge);
void fillmiginfoperiods (int ci, struct edgemiginfo *em);
void fillmiginfo (int ci, int li, struct edge *gtree, int edge, int sisedge);
void copynewmig_to_gtree (int ci, int li);
void copynewmighg_to_gtree (int ci, int li); //moved into update_hg.cpp
void storegenealogystats (int ci, int li, int mode);
int picktopop (int nowpop, int plist[], int numpops);  
int picktopop2 (int nowpop, int plist[], int numpops, int notother); 
int mwork_single_edge (int ci, struct edgemiginfo *edgem,struct edgemiginfo *oldedgem, int lastmigperiod); 
int mwork_two_edges(int ci, struct edgemiginfo *edgem, struct edgemiginfo *sisem, 
      struct edgemiginfo *oldedgem, struct edgemiginfo *oldsisem,int lastmigperiod, int* empall, int* smpall);
double getnewt (double t_u_prior, double t_d_prior, double oldt, double wadjust);

__END_DECLS
#endif /* _UPDATE_GTREE_COMMON_H_ */
